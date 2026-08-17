"""Boore & Joyner (1997) / Boore (2016) generic rock velocity profiles.

Boore and Joyner (1997; BJ97) defined two reference shear-wave velocity
profiles used widely for computing generic crustal amplifications:

- ``BJ97gr``, the "generic rock" profile, with a 30 m time-averaged
  velocity :math:`\\bar V_S(30\\ \\mathrm{m})` of 618 m/s.
- ``BJ97gvhr``, the "generic very hard rock" profile, with
  :math:`\\bar V_S(30\\ \\mathrm{m})` of 2780 m/s.

Boore (2016) interpolates these two profiles (in slowness space) to derive a
third profile, ``BJ97gr760``, with the more commonly used reference velocity
:math:`\\bar V_S(30\\ \\mathrm{m})` = 760 m/s, and tabulates it (with a
companion density model) in Table 1 of that article. ``BJ97gr`` itself is
reconstructed here algebraically by inverting Boore's (2016) interpolation
equations using the tabulated ``BJ97gr760`` and ``BJ97gvhr`` profiles; this
recovers :math:`\\bar V_S(30\\ \\mathrm{m})` of about 616 m/s (versus the
618 m/s target), the small difference being due to the decimation of the
tabulated profiles.

References
----------
Boore, D. M., and W. B. Joyner (1997). Site amplifications for generic rock
sites, *BSSA*, 87(2), 327-341.

Boore, D. M. (2016). Determining generic velocity and density models for
crustal amplification calculations, with an update of the Boore and Joyner
(1997) generic site amplification for :math:`\\bar V_S(Z)` = 760 m/s,
*BSSA*, 106(1), 316-320.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .. import model
from ..contracts import VelocityProfile

# Table 1 of Boore (2016): a decimated (36-point) BJ97gr760 profile with
# V-bar_S(30 m) = 759 m/s (760 m/s for the undecimated, 273-point model).
_BJ97GR760 = model.load_data_file("boore_2016_bj97gr760.csv")

# Boore & Joyner (1997) generic very hard rock (BJ97gvhr) profile, with
# V-bar_S(30 m) = 2780 m/s, sampled every 50 m to a depth of 750 m.
_BJ97GVHR = model.load_data_file("boore_joyner_1997_bj97gvhr.csv")

# Reference 30 m time-averaged velocities [km/s] used to derive the
# interpolation coefficient beta (Boore, 2016, eq. 7).
_VBAR_30_BJ97GR = 0.618
_VBAR_30_BJ97GVHR = 2.780
_VBAR_30_BJ97GR760 = 0.759  # matches the decimated table used here
_BETA = (1 / _VBAR_30_BJ97GR760 - 1 / _VBAR_30_BJ97GR) / (
    1 / _VBAR_30_BJ97GVHR - 1 / _VBAR_30_BJ97GR
)

_MAX_DEPTH_BJ97GVHR_M = _BJ97GVHR.depth_km.max() * 1000


def bj97gr760_profile(depth: npt.ArrayLike) -> VelocityProfile:
    """Generic rock velocity profile BJ97gr760 from Boore (2016).

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile.

    Returns
    -------
    VelocityProfile
        Dataclass with ``depth``, ``vs_median``, and ``std_vs_ln``.
        ``std_vs_ln`` is zero because this is a single deterministic
        reference profile (Boore, 2016, Table 1) with no reported
        log-standard-deviation for velocity.

    References
    ----------
    Boore, D. M. (2016). *BSSA*, 106(1), 316-320.
    """
    return _tabulated_profile(
        depth, _BJ97GR760.depth_km, _BJ97GR760.vel_shear_kps, site_class="BJ97gr760"
    )


def bj97gr760_density(depth: npt.ArrayLike) -> npt.NDArray[np.floating]:
    """Density profile associated with BJ97gr760 (Boore, 2016, Table 1).

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile.

    Returns
    -------
    np.ndarray
        Mass density [g/cm^3] at each depth.
    """
    depth = np.asarray(depth, dtype=float)
    depth_km = depth / 1000.0
    return np.interp(depth_km, _BJ97GR760.depth_km, _BJ97GR760.density_gcc)


def bj97gvhr_profile(depth: npt.ArrayLike) -> VelocityProfile:
    """Generic very hard rock velocity profile BJ97gvhr (Boore & Joyner, 1997).

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile. Must not exceed 750 m,
        the maximum depth of the tabulated profile.

    Returns
    -------
    VelocityProfile
        Dataclass with ``depth``, ``vs_median``, and ``std_vs_ln`` (zero;
        this is a single deterministic reference profile).

    References
    ----------
    Boore, D. M., and W. B. Joyner (1997). *BSSA*, 87(2), 327-341.
    """
    depth = np.asarray(depth, dtype=float)
    if np.any(depth > _MAX_DEPTH_BJ97GVHR_M):
        raise ValueError(
            f"depth must not exceed {_MAX_DEPTH_BJ97GVHR_M:g} m, the maximum "
            "depth of the tabulated BJ97gvhr profile."
        )
    return _tabulated_profile(
        depth, _BJ97GVHR.depth_km, _BJ97GVHR.vel_shear_kps, site_class="BJ97gvhr"
    )


def bj97gr_profile(depth: npt.ArrayLike) -> VelocityProfile:
    """Generic rock velocity profile BJ97gr (Boore & Joyner, 1997).

    This profile is reconstructed by inverting the slowness-interpolation
    procedure of Boore (2016, eqs. 5-7), using the tabulated ``BJ97gr760``
    and ``BJ97gvhr`` profiles::

        S760(z) = (1 - beta) * Sgr(z) + beta * Sgvhr(z)

    solved for the BJ97gr slowness ``Sgr(z) = 1 / Vgr(z)``, with ``beta``
    fixed by matching the 30 m time-averaged velocities of the three
    profiles (618, 2780, and 759 m/s, respectively). Because the
    tabulated ``BJ97gvhr`` profile only extends to 750 m, this
    reconstruction is limited to that depth range.

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile. Must not exceed 750 m.

    Returns
    -------
    VelocityProfile
        Dataclass with ``depth``, ``vs_median``, and ``std_vs_ln`` (zero;
        this is a single deterministic reference profile).

    References
    ----------
    Boore, D. M., and W. B. Joyner (1997). *BSSA*, 87(2), 327-341.

    Boore, D. M. (2016). *BSSA*, 106(1), 316-320.
    """
    depth = np.asarray(depth, dtype=float)
    if np.any(depth > _MAX_DEPTH_BJ97GVHR_M):
        raise ValueError(
            f"depth must not exceed {_MAX_DEPTH_BJ97GVHR_M:g} m; the BJ97gr "
            "profile can only be reconstructed where the tabulated BJ97gvhr "
            "profile is available."
        )
    depth_km = depth / 1000.0

    s_combined = 1 / np.interp(depth_km, _BJ97GR760.depth_km, _BJ97GR760.vel_shear_kps)
    s_gvhr = 1 / np.interp(depth_km, _BJ97GVHR.depth_km, _BJ97GVHR.vel_shear_kps)
    s_gr = (s_combined - _BETA * s_gvhr) / (1 - _BETA)

    vs_median = 1000 / s_gr
    std_vs_ln = np.zeros_like(vs_median)

    return VelocityProfile(
        depth=depth, vs_median=vs_median, std_vs_ln=std_vs_ln, site_class="BJ97gr"
    )


def _tabulated_profile(
    depth: npt.ArrayLike,
    depth_km_table: npt.NDArray[np.floating],
    vs_kps_table: npt.NDArray[np.floating],
    site_class: str,
) -> VelocityProfile:
    depth = np.asarray(depth, dtype=float)
    depth_km = depth / 1000.0
    vs_median = 1000 * np.interp(depth_km, depth_km_table, vs_kps_table)
    std_vs_ln = np.zeros_like(vs_median)
    return VelocityProfile(
        depth=depth, vs_median=vs_median, std_vs_ln=std_vs_ln, site_class=site_class
    )
