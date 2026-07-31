"""Menq (2003) nonlinear soil model for gravelly soils."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from ._hyperbolic import ModifiedHyperbolicBase
from ._units import convert_units

_KPA_TO_ATM = 1.0 / 101.325


class MenqSoilType(ModifiedHyperbolicBase):
    """Menq (2003) model for gravelly soils.

    Parameters
    ----------
    name : str, optional
        Identification label.
    unit_wt : float, default 0
        Unit weight [kN/m³].
    coef_unif : float, default 10
        Uniformity coefficient (Cᵤ).
    diam_mean : float, default 5
        Mean diameter D₅₀ [mm].
    stress_mean : float, default 101.3
        Mean effective stress [kN/m²].
    num_cycles : float, default 10
        Number of loading cycles.
    damping_min : float or None
        Minimum damping [decimal]; computed when *None*.
    strains : array_like or None
        Shear strain levels [decimal]; defaults to ``np.logspace(-6, -1.5, 20)``.
    """

    @convert_units(
        unit_wt="kilonewton / meter ** 3",
        stress_mean="kilopascal",
        diam_mean="millimeter",
        strains="dimensionless",
    )
    def __init__(
        self,
        name: str = "",
        unit_wt: float = 0.0,
        coef_unif: float = 10,
        diam_mean: float = 5,
        stress_mean: float = 101.3,
        num_cycles: float = 10,
        damping_min: float | None = None,
        strains: npt.ArrayLike | None = None,
    ) -> None:
        self._coef_unif = coef_unif
        self._diam_mean = diam_mean
        self._stress_mean = stress_mean
        self._num_cycles = num_cycles

        if damping_min is None:
            damping_min = self._calc_damping_min()

        if not name:
            name = self._create_name()

        super().__init__(name, unit_wt, damping_min, strains)

    def _calc_damping_min(self) -> float:
        return (
            0.55
            * self._coef_unif**0.1
            * self._diam_mean**-0.3
            * (self._stress_mean * _KPA_TO_ATM) ** -0.08
        ) / 100

    @property
    def masing_scaling(self) -> float:
        return 0.6329 - 0.00566 * np.log(self._num_cycles)

    @property
    def strain_ref(self) -> float:
        return (
            0.12
            * self._coef_unif**-0.6
            * (self._stress_mean * _KPA_TO_ATM) ** (0.5 * self._coef_unif**-0.15)
        ) / 100

    @property
    def curvature(self) -> float:
        return 0.86 + 0.1 * np.log10(self._stress_mean * _KPA_TO_ATM)

    def _create_name(self) -> str:
        return (
            f"Menq (Cᵤ={self._coef_unif:.1f}, D₅₀={self._diam_mean:.1f} mm, "
            f"σₘ'={self._stress_mean:.1f} kN/m²)"
        )
