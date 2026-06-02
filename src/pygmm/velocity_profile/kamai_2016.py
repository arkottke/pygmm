"""Kamai et al. (2016) shear-wave velocity profile model."""

from __future__ import annotations

import warnings

import numpy as np
import numpy.typing as npt

from ..contracts import VelocityProfile


def kea16_profile(depth: npt.ArrayLike, vs30: float, region: str) -> VelocityProfile:
    """Median shear-wave velocity profile from Kamai et al. (2016).

    Parameters
    ----------
    depth : array_like
        Depths [m] at which to evaluate the profile.
    vs30 : float
        Time-averaged shear-wave velocity in the upper 30 m [m/s].
    region : str
        ``'california'`` or ``'japan'``.

    Returns
    -------
    VelocityProfile
        Dataclass with ``depth``, ``vs_median``, and ``std_vs_ln``.

    References
    ----------
    Kamai, R., et al. (2016). Nonlinear site response from the KiK-net database.
    *BSSA*, 106(4), 1710–1723.
    """
    region = region.lower()
    if region == "california":
        b1, b2, b3, b4, b5 = 0.25, 104, 0.22, 166, 53
        b6, b7, b8 = -0.00388, 0.0002, 0.195
    elif region == "japan":
        b1, b2, b3, b4, b5 = 0.2, 40, 0.25, 260, 28
        b6, b7, b8 = -0.00296, 0, 0.4
    else:
        raise ValueError("region must be 'california' or 'japan'")

    if not (250 <= vs30 <= 850):
        warnings.warn(
            "kea16_profile is recommended for Vs30 between 250 and 850 m/s.",
            UserWarning,
            stacklevel=2,
        )

    depth = np.asarray(depth, dtype=float)
    a0 = b1 * vs30 + b2
    a1 = b3 * vs30 + b4
    a2 = b5 * np.exp(b6 * vs30)

    vs_median = a0 + a1 * np.log(np.maximum((depth + a2) / a2, 1e-9))
    std_vs_ln = np.full_like(vs_median, b7 * vs30 + b8)

    return VelocityProfile(depth=depth, vs_median=vs_median, std_vs_ln=std_vs_ln)
