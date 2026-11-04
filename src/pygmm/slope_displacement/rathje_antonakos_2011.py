"""Rathje and Antonakos (2011) flexible sliding-mass displacement model."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm
from .rathje_saygili_2008 import calc_disp_rs08

__author__ = "Albert Kottke"


@dist_lognorm
def calc_disp_ra11(
    yield_coef: float,
    pga: float | np.ndarray,
    period_slide: float,
    period_mean: float,
    mag: float | None = None,
    pgv: float | None = None,
    **kwargs,
):
    """Rathje and Antonakos (2011) slope displacement model.

    Extends :func:`~pygmm.slope_displacement.calc_disp_rs08` to flexible
    sliding masses by scaling PGA (and PGV) to the maximum seismic
    coefficient (and its velocity) of the sliding mass.

    Parameters
    ----------
    yield_coef : float
        Slope yield coefficient [g].
    pga : float or array_like
        Peak ground acceleration [g].
    period_slide : float
        Fundamental period of the sliding mass [s].
    period_mean : float
        Mean period of the ground motion [s].
    mag : float, optional
        Moment magnitude. Selects the scalar (PGA, M) model.
    pgv : float, optional
        Peak ground velocity [cm/s]. Selects the vector (PGA, PGV) model.
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    if (mag is None) == (pgv is None):
        raise NotImplementedError("Provide only one of `mag` or `pgv`.")

    period_ratio = period_slide / period_mean
    if period_ratio < 0.1:
        k_max = pga
    else:
        ln_pr = np.log(period_ratio / 0.1)
        ln_ratio = (0.459 - 0.702 * pga) * ln_pr + (-0.228 + 0.076 * pga) * ln_pr**2
        k_max = pga * np.exp(ln_ratio)

    # Compute median
    if mag is not None:
        ln_disp = calc_disp_rs08(yield_coef, k_max, mag=mag, stats=True)[0]
        f_1 = np.minimum(3.69 * period_slide - 1.22 * period_slide**2, 2.78)
        ln_disp += f_1
    else:
        assert pgv is not None
        if period_ratio < 0.2:
            k_velmax = pgv
        else:
            ln_pr = np.log(period_ratio / 0.2)
            ln_ratio = 0.240 * ln_pr + (-0.091 - 0.171 * pga) * ln_pr**2
            k_velmax = pgv * np.exp(ln_ratio)

        ln_disp = calc_disp_rs08(yield_coef, k_max, pgv=k_velmax, stats=True)[0]

        f_2 = np.minimum(1.42 * period_slide, 0.71)
        ln_disp += f_2

    # Compute the standard deviation
    if mag is not None:
        ln_std = 0.694 + 0.322 * (yield_coef / k_max)
    else:
        ln_std = 0.400 + 0.284 * (yield_coef / k_max)

    return ln_disp, ln_std
