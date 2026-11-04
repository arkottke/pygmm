"""Wair, DeJong, and Shantz (2012) shear-wave velocity from SPT blow count."""

from __future__ import annotations

from typing import Any

import numpy as np
import numpy.typing as npt

from .._utils import KPA_TO_ATM, check_options

__author__ = "Albert Kottke"

# Section 4.5: coefficient, blow count exponent, and stress exponent, plus the
# age scaling factors.
_COEFFS: dict[str, dict[str, Any]] = {
    "fine_grained": {
        "coeffs": (26, 0.17, 0.32),
        "asf": {"holocene": 0.88, "pleistocene": 1.12},
    },
    "sand": {
        "coeffs": (30, 0.23, 0.23),
        "asf": {"holocene": 0.90, "pleistocene": 1.17},
    },
    "gravel": {
        "coeffs": (78, 0.19, 0.18),
        "asf": {"holocene": 53 / 78, "pleistocene": 115 / 78},
    },
    "all": {
        "coeffs": (30, 0.215, 0.275),
        "asf": {"holocene": 0.87, "pleistocene": 1.13},
    },
}


def calc_vel_shear_spt_wds12(
    blows: npt.ArrayLike,
    stress_vert_eff: npt.ArrayLike,
    soil_type: str = "all",
    age: str | None = None,
) -> npt.NDArray[np.floating]:
    """Estimate the shear-wave velocity from the SPT blow count.

    Equations from Section 4.5 of Wair, DeJong, and Shantz (2012).

    Parameters
    ----------
    blows : array_like
        Blow counts, :math:`N_{60}`.
    stress_vert_eff : array_like
        Vertical effective stress [atm].
    soil_type : str, default "all"
        Soil type: ``"fine_grained"``, ``"sand"``, ``"gravel"``, or ``"all"``.
    age : str, optional
        Age of the soil: ``"holocene"`` or ``"pleistocene"``. If *None*, no
        age scaling is applied.

    Returns
    -------
    np.ndarray
        Shear-wave velocity [m/s].
    """
    check_options(soil_type, _COEFFS, "soil_type")
    # Convert to kPa
    stress_kpa = np.asarray(stress_vert_eff, dtype=float) / KPA_TO_ATM

    coeff, pow_blows, pow_stress = _COEFFS[soil_type]["coeffs"]
    vel_shear = coeff * np.asarray(blows) ** pow_blows * stress_kpa**pow_stress

    if age is not None:
        asf = _COEFFS[soil_type]["asf"]
        check_options(age, asf, "age")
        vel_shear = vel_shear * asf[age]

    return vel_shear
