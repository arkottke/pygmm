"""Bray and Travasarou (2007) displacement model for shallow crustal events."""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .._utils import dist_lognorm

__author__ = "Albert Kottke"


@dist_lognorm
def calc_disp_bt07(yield_coef: float, period_slide: float, psa_dts: float, **kwds):
    """Bray and Travasarou (2007) median non-zero displacement.

    Parameters
    ----------
    yield_coef : float
        Slope yield coefficient [g].
    period_slide : float
        Fundamental period of the sliding mass, :math:`T_s` [s].
    psa_dts : float
        5%-damped spectral acceleration at the degraded period,
        :math:`1.5 T_s` [g].
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    ln_yield_coef = np.log(yield_coef)
    ln_psa_dts = np.log(psa_dts)

    # Equation (4)
    ln_mean = -0.22 if period_slide < 0.05 else -1.10
    ln_mean += (
        -2.83 * ln_yield_coef
        - 0.333 * ln_yield_coef**2
        + 0.566 * ln_yield_coef * ln_psa_dts
        + 3.04 * ln_psa_dts
        - 0.244 * ln_psa_dts**2
        + 1.5 * period_slide
    )
    ln_std = 0.66
    return ln_mean, ln_std


def calc_prob_disp_bt07(
    yield_coef: float, period_slide: float, psa_dts: float, **kwds
) -> float:
    """Bray and Travasarou (2007) probability of non-zero displacement.

    Parameters are as for :func:`calc_disp_bt07`.

    Returns
    -------
    float
        Probability of a non-zero displacement.
    """
    # Modified from Equation (3)
    ln_yield_coef = np.log(yield_coef)

    return norm.cdf(
        -1.76
        - 3.22 * ln_yield_coef
        - 0.484 * period_slide * ln_yield_coef
        + 3.52 * np.log(psa_dts)
    )
