"""Park and Kishida (2019) fundamental period of earth dams."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm

__author__ = "Albert Kottke"

_COEFFS = {
    "longitudinal": (-2.629, 0.377, 0.464),
    "transverse": (-2.685, 0.430, 0.375),
    "vertical": (-2.793, 0.283, 0.456),
}
_ALIASES = {"lon": "longitudinal", "trans": "transverse", "vert": "vertical"}


@dist_lognorm
def calc_dam_period_pk19(height: float, direction: str = "transverse", **kwds):
    """Fundamental period of an earth dam from Park and Kishida (2019).

    Parameters
    ----------
    height : float
        Height of the dam [m].
    direction : str, default "transverse"
        Direction of the mode: ``"longitudinal"``, ``"transverse"``, or
        ``"vertical"`` (or the abbreviations ``"lon"``, ``"trans"``,
        ``"vert"``).
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median period [s], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    direction = _ALIASES.get(direction, direction)
    try:
        a, b, ln_std = _COEFFS[direction]
    except KeyError:
        raise NotImplementedError(f"Unknown direction: {direction}") from None

    ln_period = a + b * np.log(height)

    return ln_period, ln_std
