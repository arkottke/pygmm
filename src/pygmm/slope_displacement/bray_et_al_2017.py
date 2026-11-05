"""Bray et al. (2017) displacement model for subduction-zone events."""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from .._utils import dist_lognorm

__author__ = "Albert Kottke"


@dist_lognorm
def calc_disp_bea17(
    yield_coef: float, period_slide: float, psa_dts: float, mag: float, **kwds
):
    """Bray et al. (2017) median non-zero displacement for subduction events.

    Parameters
    ----------
    yield_coef : float
        Slope yield coefficient [g].
    period_slide : float
        Fundamental period of the sliding mass, :math:`T_s` [s].
    psa_dts : float
        5%-damped spectral acceleration at the degraded period,
        :math:`1.5 T_s` [g].
    mag : float
        Moment magnitude.
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    ln_yield_coef = np.log(yield_coef)
    ln_psa_dts = np.log(psa_dts)

    ln_mean = (
        -3.535 * ln_yield_coef
        - 0.390 * ln_yield_coef**2
        + 0.538 * ln_yield_coef * ln_psa_dts
        + 3.060 * ln_psa_dts
        - 0.225 * ln_psa_dts**2
        + 0.550 * mag
    )
    if np.isclose(period_slide, 0.0):
        ln_mean -= 5.864
    else:
        if period_slide < 0.1:
            a1, a2, a3 = -5.864, -9.421, 0.0
        else:
            a1, a2, a3 = -6.896, 3.081, -0.803
        ln_mean += a1 + a2 * period_slide + a3 * period_slide**2

    ln_std = 0.73
    return ln_mean, ln_std


def calc_prob_disp_bea17(
    yield_coef: float, period_slide: float, psa_dts: float, **kwds
) -> float:
    """Bray et al. (2017) probability of non-zero displacement.

    Parameters are as for :func:`calc_disp_bea17`.

    Returns
    -------
    float
        Probability of a non-zero displacement.
    """
    ln_yield_coef = np.log(yield_coef)

    if period_slide <= 0.7:
        z = (
            -2.64
            - 3.20 * ln_yield_coef
            - 0.17 * ln_yield_coef**2
            - 0.49 * period_slide * ln_yield_coef
            + 2.09 * period_slide
            + 2.91 * np.log(psa_dts)
        )
    else:
        z = (
            -3.53
            - 4.78 * ln_yield_coef
            - 0.34 * ln_yield_coef**2
            - 0.30 * period_slide * ln_yield_coef
            - 0.67 * period_slide
            + 2.66 * np.log(psa_dts)
        )

    return norm.cdf(z)
