"""Cho and Rathje (2022) sliding-mass displacement model."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm
from ..model import Coefficients

__author__ = "Albert Kottke"

_COEFFS_PGV = Coefficients(
    b_0=-1.01,
    b_1=1.57,
    b_2=-0.25,
    c_0=0.81,
    c_1=-1.05,
    c_2=-0.60,
    b_3=-4.50,
    b_4=-1.37,
    c_3=1.51,
    c_4=0.10,
    p_0=0.26,
    p_1=-0.14,
    p_2=0.40,
    p_3=0.16,
)
_COEFFS_PGA = Coefficients(
    b_0=3.94,
    b_1=1.28,
    b_2=-0.27,
    c_0=1.21,
    c_1=0.29,
    c_2=0.25,
    b_3=1.96,
    b_4=-0.89,
    c_3=1.47,
    c_4=0.14,
    e_0=0.56,
    e_1=0.66,
)


@dist_lognorm
def calc_disp_cr22(
    yield_coef: float,
    period_slide: float,
    height_ratio: float,
    pgv: float | None = None,
    pga: float | None = None,
    mag: float | None = None,
    **kwds,
):
    """Displacement model from Cho and Rathje (2022).

    The PGV model is used unless `pga` (and `mag`) are provided.

    Parameters
    ----------
    yield_coef : float
        Yield coefficient of the slope [g].
    period_slide : float
        Period of the sliding mass [s].
    height_ratio : float
        Height ratio of the sliding mass.
    pgv : float, optional
        Peak ground velocity [cm/s].
    pga : float, optional
        Peak ground acceleration [g].
    mag : float, optional
        Moment magnitude. Required with `pga`.
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    ln_yield_coef = np.log(yield_coef)
    ln_period_slide = np.log(period_slide)

    pgv_model = pga is None
    C = _COEFFS_PGV if pgv_model else _COEFFS_PGA

    # Mean model
    if height_ratio <= 0.6:
        a_0 = C.b_0 + C.b_1 * ln_period_slide + C.b_2 * ln_yield_coef
        a_1 = C.c_0 + C.c_1 * ln_period_slide + C.c_2 * ln_period_slide**2
    else:
        a_0 = C.b_3 + C.b_4 * ln_yield_coef
        a_1 = C.c_3 + C.c_4 * ln_yield_coef

    if pgv_model:
        if pgv is None:
            raise ValueError("Provide `pgv`, or `pga` and `mag`.")
        ln_mean = a_0 + a_1 * np.log(pgv)
    else:
        if pga is None or mag is None:
            raise ValueError("`mag` is required with `pga`.")
        d_0 = C.e_0 if height_ratio <= 0.6 else C.e_1
        ln_mean = a_0 + a_1 * np.log(pga) + d_0 * (mag - 6.5)

    # Uncertainty
    if pgv_model:
        if height_ratio <= 0.6:
            phi = C.p_0 + C.p_1 * ln_period_slide
        else:
            flag_ky = 0 if yield_coef < 0.2 else 1
            phi = C.p_2 + C.p_3 * np.log(yield_coef / 0.2) * flag_ky
        tau = 0.139
    else:
        phi = 0.72
        tau = 0.143

    ln_std = np.sqrt(phi**2 + tau**2)

    return ln_mean, ln_std
