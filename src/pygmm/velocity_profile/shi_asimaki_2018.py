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


def _svm_params(vs30: float) -> tuple[float, float, float]:
    """Eq. (4): SVM parameters ``(vs_0, k, n)`` as a function of Vs30."""
    vs_0 = _P1 * vs30**2 + _P2 * vs30 + _P3
    k = np.exp(_R1 * vs30**_R2 + _R3)
    n = _S1 * np.exp(_S2 * vs30) + _S3 * np.exp(_S4 * vs30)
    return float(vs_0), float(k), float(n)


def sa18_mean_vs(depth: npt.ArrayLike, vs30: float) -> npt.NDArray[np.floating]:
    """Arithmetic-mean shear-wave velocity from Shi and Asimaki (2018).

    This is eq. (3) as fit by the authors to the *arithmetic mean* of the
    (lognormal) measured velocities, see eq. (1). Use :func:`sa18_profile` for
    the median.

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile.
    vs30 : float
        Time-averaged shear-wave velocity in the upper 30 m [m/s].

    Returns
    -------
    numpy.ndarray
        Mean shear-wave velocity [m/s] at each depth.
    """
    vs_0, k, n = _svm_params(vs30)
    dz = np.maximum(np.asarray(depth, dtype=float) - _DEPTH_0, 0.0)
    return vs_0 * (1 + k * dz) ** (1 / n)


def sa18_depth_to_vs(vs: npt.ArrayLike, vs30: float) -> npt.NDArray[np.floating]:
    """Depth at which the mean SVM profile reaches a shear-wave velocity.

    Inverts eq. (3) of Shi and Asimaki (2018). For example, the authors
    recommend using the profile only down to a velocity of 1000 m/s.

    Parameters
    ----------
    vs : array_like
        Shear-wave velocity [m/s]. Must be at least the near-surface velocity
        :math:`V_{S0}` of the model.
    vs30 : float
        Time-averaged shear-wave velocity in the upper 30 m [m/s].

    Returns
    -------
    numpy.ndarray
        Depth [m] of the mean profile at each velocity.
    """
    vs_0, k, n = _svm_params(vs30)
    vs = np.asarray(vs, dtype=float)
    if np.any(vs < vs_0):
        raise ValueError(
            f"vs must be at least the near-surface velocity ({vs_0:.1f} m/s)."
        )
    return _DEPTH_0 + ((vs / vs_0) ** n - 1) / k


def sa18_profile(depth: npt.ArrayLike, vs30: float) -> VelocityProfile:
    """Median shear-wave velocity profile from Shi and Asimaki (2018).

    The sediment velocity model (SVM) predicts a 1D shear-wave velocity
    profile as a function of :math:`V_{S30}` for basin sediments in
    California.

    The authors fit eq. (3) to the arithmetic mean of lognormally distributed
    velocities (:func:`sa18_mean_vs`). The returned ``vs_median`` is the
    corresponding lognormal median, ``mean * exp(-std_vs_ln**2 / 2)``. The
    standard deviation of eq. (9), which the authors define at the layer
    midpoint, is evaluated here at each ``depth``.

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

    # Eqs. (3) and (4): arithmetic-mean Vs profile.
    vs_mean = sa18_mean_vs(depth, vs30)

    # Eq. (9): standard deviation of Vs [m/s] as a function of depth and
    # Vs30, clipped at 0 to avoid unphysical negative values.
    sd_vs = np.maximum(_SD_C0 + _SD_C1 * depth + _SD_C2 * vs30, 0.0)

    # Convert the linear-space standard deviation to a log-space standard
    # deviation using the lognormal moment relationships of Toro (1995)
    # (eqs. 1 and 2), noting that the SVM profile is the arithmetic mean.
    std_vs_ln = np.sqrt(np.log1p((sd_vs / vs_mean) ** 2))

    # Median of the lognormal distribution (eq. 1 inverted).
    vs_median = vs_mean * np.exp(-(std_vs_ln**2) / 2)

    return VelocityProfile(
        depth=depth,
        vs_median=vs_median,
        std_vs_ln=std_vs_ln,
        region="california",
    )
