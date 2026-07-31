"""Rollins et al. (2020) nonlinear model for gravels."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from ._hyperbolic import ModifiedHyperbolicBase
from ._units import convert_units

_KPA_TO_ATM = 1.0 / 101.325


class RollinsEtAlSoilType(ModifiedHyperbolicBase):
    """Rollins et al. (2020) model for gravels.

    Parameters
    ----------
    unit_wt : float, default 0
        Unit weight [kN/m³].
    name : str, optional
        Identification label.
    stress_mean : float, default 101.3
        Mean effective confining pressure [kN/m²].
    coef_unif : float or None, default None
        Uniformity coefficient *Cu* = D60/D10.  Uses Eq. 8 when provided,
        Eq. 5 otherwise.
    num_cycles : float, default 10
        Number of loading cycles.
    damping_min : float or None
        Minimum damping [decimal]; defaults to 0.01 when *None*.
    strains : array_like or None
        Shear strain levels [decimal]; defaults to ``np.logspace(-6, -1.5, 20)``.
    """

    @convert_units(
        unit_wt="kilonewton / meter ** 3",
        stress_mean="kilopascal",
        strains="dimensionless",
    )
    def __init__(
        self,
        unit_wt: float = 0.0,
        name: str = "",
        stress_mean: float = 101.3,
        coef_unif: float | None = None,
        num_cycles: float = 10,
        damping_min: float | None = None,
        strains: npt.ArrayLike | None = None,
    ) -> None:
        self._stress_mean = stress_mean
        self._coef_unif = coef_unif
        self._num_cycles = num_cycles

        if damping_min is None:
            damping_min = 0.01

        if not name:
            if coef_unif is not None:
                name = "Rollins et al. (2020) - Cu={:.1f}, σ'₀={:.0f} kPa".format(
                    coef_unif, stress_mean
                )
            else:
                name = f"Rollins et al. (2020) - σ'₀={stress_mean:.0f} kPa"

        super().__init__(name, unit_wt, damping_min, strains)

    @property
    def strain_ref(self) -> float:
        if self._coef_unif is not None:
            return 0.0046 * self._coef_unif**-0.197 * self._stress_mean**0.52 / 100
        return 0.0039 * self._stress_mean**0.42 / 100

    @property
    def curvature(self) -> float:
        return 0.84

    @property
    def masing_scaling(self) -> float:
        return 0.53 - 0.0057 * np.log(self._num_cycles)
