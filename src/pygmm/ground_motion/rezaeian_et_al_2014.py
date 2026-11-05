"""Rezaeian et al. (2014) damping scaling factors."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .. import model
from .._utils import check_bounds, check_options

__author__ = "Albert Kottke"

COMPONENTS = ("rotd50", "roti50", "vertical")

_COEFFS: dict[str, np.recarray] = {}


def calc_damping_scaling_rea15(
    damping: float,
    mag: float,
    dist_rup: float,
    comp: str = "rotd50",
    periods: npt.ArrayLike | None = None,
) -> np.recarray:
    """Compute the damping scaling factor of Rezaeian et al. (2014).

    Parameters
    ----------
    damping : float
        Damping ratio of the oscillator in percent (2 for 2%). Valid from 0.5
        to 30.
    mag : float
        Moment magnitude.
    dist_rup : float
        Closest distance to the rupture [km].
    comp : str, default "rotd50"
        Component: one of ``"rotd50"``, ``"roti50"``, or ``"vertical"``.
    periods : array_like, optional
        Periods [s] at which to interpolate the values. If *None*, the periods
        of the model coefficients are used.

    Returns
    -------
    :class:`numpy.recarray`
        Record array with fields ``period``, ``dsf`` (damping scaling factor
        relative to 5%), and ``ln_std`` (natural log standard deviation).
    """
    check_bounds(damping, 0.5, 30, "damping")
    check_options(comp, COMPONENTS, "comp")

    # Only load the coefficients once
    if comp not in _COEFFS:
        _COEFFS[comp] = model.load_data_file(
            f"rezaeian_et_al_2014-{comp}.csv", skip_header=1
        )
    C = _COEFFS[comp]

    ln_damp = np.log(damping)
    ln_dsf = (
        C.b0
        + C.b1 * ln_damp
        + C.b2 * ln_damp**2
        + (C.b3 + C.b4 * ln_damp + C.b5 * ln_damp**2) * mag
        + (C.b6 + C.b7 * ln_damp + C.b8 * ln_damp**2) * np.log(dist_rup + 1)
    )
    ln_damp_5 = np.log(damping / 5)
    ln_std = np.abs(C.a0 * ln_damp_5 + C.a1 * ln_damp_5**2)

    if periods is not None:
        # Interpolate over the provided periods
        periods = np.asarray(periods, dtype=float)
        ln_dsf = np.interp(np.log(periods), np.log(C.period), ln_dsf)
        ln_std = np.interp(np.log(periods), np.log(C.period), ln_std)
    else:
        periods = C.period

    return np.rec.fromarrays(
        [periods, np.exp(ln_dsf), ln_std], names="period,dsf,ln_std"
    )
