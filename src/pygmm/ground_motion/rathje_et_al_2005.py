"""Rathje et al. (2005) frequency-content period parameters."""

from __future__ import annotations

import numpy as np

from .._utils import check_bounds, check_options, dist_lognorm

__author__ = "Albert Kottke"

KINDS = ("period_mean", "period_avg", "period_pred")

# Model coefficients from Table 2
_COEFFS = {
    kind: np.rec.fromrecords([values], names="c1,c2,c3,c4,c5,c6")[0]
    for kind, values in zip(
        KINDS,
        [
            (-1.00, 0.18, 0.0038, 0.078, 0.27, 0.40),
            (-0.89, 0.29, 0.0030, 0.070, 0.25, 0.37),
            (-1.78, 0.30, 0.0045, 0.150, 0.33, 0.24),
        ],
        strict=True,
    )
}

_INTRA = {"b": 0.42, "c": 0.38, "d": 0.31}
_INTER = {"period_mean": 0.17, "period_avg": 0.13, "period_pred": 0.22}


@dist_lognorm
def calc_period_rea05(
    kind: str,
    mag: float,
    dist_rup: float,
    site_class: str = "c",
    directivity: bool = False,
    **kwargs,
):
    """Rathje et al. (2005) frequency-content period parameters.

    Parameters
    ----------
    kind : str
        One of ``"period_mean"`` (mean period, :math:`T_m`),
        ``"period_avg"`` (average spectral period, :math:`T_{avg}`), or
        ``"period_pred"`` (predominant period, :math:`T_p`).
    mag : float
        Moment magnitude.
    dist_rup : float
        Closest distance to the rupture [km].
    site_class : str, default "c"
        Site class: ``"b"``, ``"c"``, or ``"d"``.
    directivity : bool, default False
        If forward-directivity effects are included.
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median period [s], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    check_options(kind, KINDS, "kind")
    check_options(site_class, "bcd", "site_class")
    C = _COEFFS[kind]

    if kind in ["period_avg", "period_pred"]:
        check_bounds(mag, 4.7, 7.6, "mag")
    else:
        mag = np.minimum(mag, 7.25)
        check_bounds(mag, 5.0, 7.25, "mag")

    s_c = int(site_class == "c")
    s_d = int(site_class == "d")
    f_d = int(directivity)

    ln_period = (
        C.c1
        + C.c2 * (mag - 6)
        + C.c3 * dist_rup
        + C.c4 * s_c
        + C.c5 * s_d
        + C.c6 * (1 - dist_rup / 20) * f_d
    )

    ln_std = np.sqrt(_INTRA[site_class] ** 2 + _INTER[kind] ** 2)

    return ln_period, ln_std
