"""Boore, Thompson, and Cadet (2011) generic Vs(z)-Vs30 relation."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .. import model
from ..contracts import VelocityProfile

# Depth-dependent quadratic-in-ln(Vs30) coefficients (Boore, Thompson, and
# Cadet, 2011), tabulated for depths of 5-29 m in 1 m increments. By
# construction, the relation converges to Vz = Vs30 (c0 = 0, c1 = 1, c2 = 0)
# as depth approaches 30 m.
_COEFF = model.load_data_file("boore_thompson_cadet_2011.csv")

_MIN_DEPTH_M = _COEFF.depth_m.min()
_MAX_DEPTH_M = _COEFF.depth_m.max()


def btc11_profile(depth: npt.ArrayLike, vs30: float) -> VelocityProfile:
    """Generic time-averaged velocity profile from Boore et al. (2011).

    Gives the time-averaged shear-wave velocity :math:`\\bar V_S(z)` to
    depth ``z`` as a quadratic function of :math:`\\ln(V_{S30})`::

        ln(Vz) = c0 + c1 * ln(Vs30) + c2 * ln(Vs30) ** 2

    with depth-dependent coefficients regressed from a global database of
    velocity profiles.

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the time-averaged velocity. Must be
        between 5 and 29 m, the range over which coefficients were
        published.
    vs30 : float
        Time-averaged shear-wave velocity in the upper 30 m [m/s].

    Returns
    -------
    VelocityProfile
        Dataclass with ``depth``, ``vs_median`` (the time-averaged
        velocity to each depth), and ``std_vs_ln`` (the reported
        log-standard deviation).

    References
    ----------
    Boore, D. M., E. M. Thompson, and H. Cadet (2011). Regional correlations
    of VS30 and velocities averaged over depths less than and greater than
    30 meters, *BSSA*, 101(6), 3046-3059.
    """
    depth = np.asarray(depth, dtype=float)
    if np.any(depth < _MIN_DEPTH_M) or np.any(depth > _MAX_DEPTH_M):
        raise ValueError(
            f"depth must be between {_MIN_DEPTH_M:g} and {_MAX_DEPTH_M:g} m."
        )

    ln_vs30 = np.log(vs30)
    c0 = np.interp(depth, _COEFF.depth_m, _COEFF.c0)
    c1 = np.interp(depth, _COEFF.depth_m, _COEFF.c1)
    c2 = np.interp(depth, _COEFF.depth_m, _COEFF.c2)
    std_vs_ln = np.interp(depth, _COEFF.depth_m, _COEFF["std"])

    vs_median = np.exp(c0 + c1 * ln_vs30 + c2 * ln_vs30**2)

    return VelocityProfile(depth=depth, vs_median=vs_median, std_vs_ln=std_vs_ln)
