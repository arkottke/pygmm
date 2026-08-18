"""Shi and Asimaki (2018) sediment velocity model (SVM)."""

from __future__ import annotations

import warnings

import numpy as np
import numpy.typing as npt

from ..contracts import VelocityProfile

# Depth [m] separating the constant near-surface layer from the power-law
# increase with depth.
_DEPTH_0 = 2.5

# Coefficients for V_S0 (eq. 4)
_P1, _P2, _P3 = -2.1688e-4, 0.5182, 69.452
# Coefficients for k (eq. 4)
_R1, _R2, _R3 = -59.67, -0.2722, 11.132
# Coefficients for n (eq. 4)
_S1, _S2, _S3, _S4 = 4.110, -1.0521e-4, -10.827, -7.6187e-3

# Coefficients for the standard deviation of V_S (eq. 9)
_SD_C0, _SD_C1, _SD_C2 = -89.7085, 1.6434, 0.5204

# Recommended and maximum valid range of Vs30 [m/s] (see Validation of the
# SVM section).
_VS30_RECOMMENDED = (173.1, 1000.0)
_VS30_MAX = 1500.0


def sa18_profile(depth: npt.ArrayLike, vs30: float) -> VelocityProfile:
    """Median shear-wave velocity profile from Shi and Asimaki (2018).

    The sediment velocity model (SVM) predicts a 1D shear-wave velocity
    profile as a function of :math:`V_{S30}` for basin sediments in
    California.

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile.
    vs30 : float
        Time-averaged shear-wave velocity in the upper 30 m [m/s].

    Returns
    -------
    VelocityProfile
        Dataclass with ``depth``, ``vs_median``, and ``std_vs_ln``.

    References
    ----------
    Shi, J., & Asimaki, D. (2018). A generic velocity profile for basin
    sediments in California conditioned on Vs30. Seismological Research
    Letters, 89(4), 1397-1409.
    """
    if vs30 > _VS30_MAX:
        raise ValueError(f"vs30 must not exceed {_VS30_MAX} m/s for sa18_profile.")
    if vs30 < _VS30_RECOMMENDED[0] or vs30 > _VS30_RECOMMENDED[1]:
        warnings.warn(
            "sa18_profile is recommended for Vs30 between "
            f"{_VS30_RECOMMENDED[0]} and {_VS30_RECOMMENDED[1]} m/s.",
            UserWarning,
            stacklevel=2,
        )

    depth = np.asarray(depth, dtype=float)

    # Eq. (4): SVM parameters as a function of Vs30.
    vs_0 = _P1 * vs30**2 + _P2 * vs30 + _P3
    k = np.exp(_R1 * vs30**_R2 + _R3)
    n = _S1 * np.exp(_S2 * vs30) + _S3 * np.exp(_S4 * vs30)

    # Eq. (3): median Vs profile, constant for depths shallower than
    # _DEPTH_0 and a power-law increase below it.
    dz = np.maximum(depth - _DEPTH_0, 0.0)
    vs_median = vs_0 * (1 + k * dz) ** (1 / n)

    # Eq. (9): standard deviation of Vs [m/s] as a function of depth and
    # Vs30, clipped at 0 to avoid unphysical negative values.
    sd_vs = np.maximum(_SD_C0 + _SD_C1 * depth + _SD_C2 * vs30, 0.0)

    # Convert the linear-space standard deviation to a log-space standard
    # deviation using the lognormal moment relationships of Toro (1995)
    # (eqs. 1 and 2), noting that `vs_median` is the arithmetic mean.
    std_vs_ln = np.sqrt(np.log1p((sd_vs / vs_median) ** 2))

    return VelocityProfile(
        depth=depth,
        vs_median=vs_median,
        std_vs_ln=std_vs_ln,
        region="california",
    )
