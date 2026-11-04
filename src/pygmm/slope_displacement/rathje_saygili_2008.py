"""Rathje and Saygili (2008) rigid sliding-block displacement model."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm

__author__ = "Albert Kottke"

_COEFFS = {
    method: np.rec.fromrecords([values], names="a1,a2,a3,a4,a5,a6,a7")[0]
    for method, values in zip(
        ["mag", "pgv"],
        [
            (4.89, -4.85, -19.64, 42.49, -29.06, 0.72, 0.89),
            (-1.56, -4.58, -20.84, 44.75, -30.50, -0.64, 1.55),
        ],
        strict=True,
    )
}


@dist_lognorm
def calc_disp_rs08(
    yield_coef: float,
    pga: float | np.ndarray,
    mag: float | None = None,
    pgv: float | None = None,
    **kwargs,
):
    """Rathje and Saygili (2008) slope displacement model.

    Two models are provided, selected by which of `mag` and `pgv` is given:

    - scalar model: PGA and magnitude
    - vector model: PGA and PGV

    Parameters
    ----------
    yield_coef : float
        Slope yield coefficient [g].
    pga : float or array_like
        Peak ground acceleration [g].
    mag : float, optional
        Moment magnitude.
    pgv : float, optional
        Peak ground velocity [cm/s].
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    if (mag is None) == (pgv is None):
        raise NotImplementedError("Provide only one of `mag` or `pgv`.")
    P = _COEFFS["pgv" if mag is None else "mag"]

    yp_ratio = yield_coef / np.asarray(pga)
    ln_disp = (
        P.a1
        + P.a2 * yp_ratio
        + P.a3 * yp_ratio**2
        + P.a4 * yp_ratio**3
        + P.a5 * yp_ratio**4
        + P.a6 * np.log(pga)
    )
    if mag is not None:
        ln_disp += P.a7 * (mag - 6)
        ln_std = 0.73 + 0.79 * yp_ratio - 0.54 * yp_ratio**2
    elif pgv is not None:
        ln_disp += P.a7 * np.log(pgv)
        ln_std = 0.41 + 0.52 * yp_ratio

    return ln_disp, ln_std
