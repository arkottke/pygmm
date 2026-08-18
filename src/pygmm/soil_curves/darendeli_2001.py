"""Darendeli (2001) nonlinear soil model."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from ..registry import register
from ._hyperbolic import ModifiedHyperbolicBase
from ._units import convert_units

# 1 kPa = 1/101.325 atm
_KPA_TO_ATM = 1.0 / 101.325


@register(provides=("soil_curves",), input="kwargs")
class DarendeliSoilType(ModifiedHyperbolicBase):
    """Darendeli (2001) model for fine-grained soils.

    Parameters
    ----------
    unit_wt : float
        Unit weight [kN/m³].
    name : str, optional
        Identification label; auto-generated when empty.
    plas_index : float, default 0
        Plasticity index [percent].
    ocr : float, default 1
        Over-consolidation ratio.
    stress_mean : float, default 101.3
        Mean effective stress [kN/m²].
    freq : float, default 1
        Excitation frequency [Hz].
    num_cycles : float, default 10
        Number of loading cycles.
    damping_min : float or None
        Minimum damping at low strains [decimal]; computed when *None*.
    strains : array_like or None
        Shear strain levels [decimal]; defaults to ``np.logspace(-6, -1.5, 20)``.
    """

    @convert_units(
        unit_wt="kilonewton / meter ** 3",
        stress_mean="kilopascal",
        freq="hertz",
        strains="dimensionless",
    )
    def __init__(
        self,
        unit_wt: float = 0.0,
        name: str = "",
        plas_index: float = 0,
        ocr: float = 1,
        stress_mean: float = 101.3,
        freq: float = 1,
        num_cycles: float = 10,
        damping_min: float | None = None,
        strains: npt.ArrayLike | None = None,
    ) -> None:
        self._plas_index = plas_index
        self._ocr = ocr
        self._stress_mean = stress_mean
        self._freq = freq
        self._num_cycles = num_cycles

        if damping_min is None:
            damping_min = self._calc_damping_min()

        if not name:
            name = self._create_name()

        super().__init__(name, unit_wt, damping_min, strains)

    def _calc_damping_min(self) -> float:
        return (
            (0.8005 + 0.0129 * self._plas_index * self._ocr**-0.1069)
            * (self._stress_mean * _KPA_TO_ATM) ** -0.2889
            * (1 + 0.2919 * np.log(self._freq))
        ) / 100

    @property
    def masing_scaling(self) -> float:
        return 0.6329 - 0.00566 * np.log(self._num_cycles)

    @property
    def strain_ref(self) -> float:
        return (
            (0.0352 + 0.0010 * self._plas_index * self._ocr**0.3246)
            * (self._stress_mean * _KPA_TO_ATM) ** 0.3483
        ) / 100

    @property
    def curvature(self) -> float:
        return 0.9190

    def _create_name(self) -> str:
        return (
            f"Darendeli (PI={self._plas_index:.0f}, OCR={self._ocr:.1f}, "
            f"σₘ'={self._stress_mean:.1f} kN/m²)"
        )
