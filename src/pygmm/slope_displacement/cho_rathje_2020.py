"""Cho and Rathje (2020) sliding-mass displacement model."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm

__author__ = "Albert Kottke"

_COEFFS_PGV = np.rec.fromrecords(
    [
        (
            -3.1355,
            -0.4253,
            0.8802,
            -3.2971,
            1.3502,
            0.0313,
            0.3719,
            -0.1137,
            0.0433,
            0.1356,
        )
    ],
    names="b0,b1,b2,b3,c0,c1,e0,e1,e2,p1",
)[0]
_COEFFS_PGA = np.rec.fromrecords(
    [(2.6482, -0.2530, 0.8802, -3.2971, 1.3058, 0.0577, 0.6002)],
    names="b0,b1,b2,b3,c0,c1,d0",
)[0]


@dist_lognorm
def calc_disp_cr20(
    yield_coef: float,
    period_slide: float,
    pga: float | None = None,
    mag: float | None = None,
    pgv: float | None = None,
    **kwds,
):
    """Calculate the displacement predicted by Cho and Rathje (2020).

    Two models are provided, selected by the inputs:

    - PGV model: provide `pgv` only
    - PGA and magnitude model: provide `pga` and `mag`

    The standard deviation is only provided for the PGV model; it is NaN for
    the PGA model.

    Parameters
    ----------
    yield_coef : float
        Yield coefficient of the slope [g].
    period_slide : float
        Period of the sliding mass [s].
    pga : float, optional
        Peak ground acceleration [g].
    mag : float, optional
        Moment magnitude.
    pgv : float, optional
        Peak ground velocity [cm/s].
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    if pgv is not None and pga is None and mag is None:
        C = _COEFFS_PGV
    elif pgv is None and pga is not None and mag is not None:
        C = _COEFFS_PGA
    else:
        raise NotImplementedError("Provide either `pgv`, or `pga` and `mag`.")

    # Equation 5
    a0 = C.b0 + C.b1 * np.log(yield_coef)
    if period_slide > 0.1:
        a0 += (C.b2 + C.b3 * yield_coef) * np.log(period_slide / 0.1)

    a1 = C.c0 + C.c1 * np.log(yield_coef)

    if pgv is not None:
        ln_disp = a0 + a1 * np.log(pgv)

        p0 = C.e0
        if pgv < 2:
            p0 -= C.e1 * np.log(pgv / 2)
        elif pgv > 50:
            p0 += C.e2 * np.log(pgv / 50)

        phi = p0
        if yield_coef > 0.2:
            phi += C.p1 * np.log(yield_coef / 0.2)

        tau = 0.177
        ln_std = np.sqrt(phi**2 + tau**2)
    else:
        assert pga is not None and mag is not None
        ln_disp = a0 + a1 * np.log(pga) + C.d0 * (mag - 6.5)
        ln_std = np.nan

    return ln_disp, ln_std
