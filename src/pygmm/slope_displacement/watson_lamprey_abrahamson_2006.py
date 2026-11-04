"""Watson-Lamprey and Abrahamson (2006) sliding-block displacement model."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm

__author__ = "Albert Kottke"


def _calc_wla06_ln_a_rms(pga):
    """Calculate :math:`\\ln a_{rms}`, Equation (2)."""
    return -1.167 + 1.02 * np.log(pga)


def _calc_wla06_ln_dur_key(yield_coef: float, pga, psa_1s, mag: float):
    """Calculate the natural log of the duration the yield coefficient is exceeded."""
    ln_pga = np.log(pga)
    ln_pga_ky = np.log(pga / yield_coef)
    ln_psa_1s = np.log(psa_1s)

    return (
        -2.775
        + 0.956 * ln_pga_ky
        - 1.554 / (ln_pga_ky + 0.390)
        - 0.597 * ln_pga
        + 0.381 * ln_psa_1s
        + 0.334 * mag
    )


@dist_lognorm
def calc_disp_wla06(yield_coef: float, pga, psa_1s, mag: float, **kwds):
    """Watson-Lamprey and Abrahamson (2006) slope displacement model.

    The model is not valid for PGA less than the yield coefficient, where NaN
    is returned.

    Parameters
    ----------
    yield_coef : float
        Slope yield coefficient [g].
    pga : float or array_like
        Peak ground acceleration [g].
    psa_1s : float or array_like
        5%-damped spectral acceleration at 1 s [g].
    mag : float
        Moment magnitude.
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    # Constants from Table 1, with updated values provided by the authors.
    a1 = 5.463
    b1 = 0.451
    b2 = 0.0191
    c1 = 0.591
    d1 = 0.205
    d2 = 0.0892
    e1 = 1.039
    e2 = 0.0427
    f1 = -1.409
    f2 = 0.1

    ln_psa_1s = np.log(psa_1s)

    ln_dur_ky = _calc_wla06_ln_dur_key(yield_coef, pga, psa_1s, mag)
    ln_a_rms = _calc_wla06_ln_a_rms(pga)

    pga = np.asarray(pga)
    ln_pga_ky = np.log(pga / yield_coef)
    ln_psa_pga = np.log(psa_1s / pga)

    ln_mean = (
        a1
        + b1 * (ln_psa_1s + 0.45)
        + b2 * (ln_psa_1s + 0.45) ** 2
        + c1 * (ln_a_rms - 1)
        + d1 * ln_psa_pga
        + d2 * ln_psa_pga**2
        + e1 * (ln_dur_ky - 0.74)
        + e2 * (ln_dur_ky - 0.74) ** 2
        + 1 / (f1 * (ln_pga_ky + f2))
    )

    # Equation not valid for PGA values less than the yield coefficient
    ln_mean = np.where(pga < yield_coef, np.nan, ln_mean)
    if ln_mean.ndim == 0:
        ln_mean = float(ln_mean)

    # From Figure 5
    ln_std = 0.53

    return ln_mean, ln_std
