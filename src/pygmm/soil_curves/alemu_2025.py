"""Alemu et al. (2025) nonlinear model for transitional silts."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from ..contracts import NonlinearSoilCurves
from ..registry import register
from ._units import convert_units

_KPA_TO_ATM = 1.0 / 101.325


@register(provides=("soil_curves",), input="kwargs")
class AlemuEtAlSoilType:
    """Alemu et al. (2025) model for transitional silts.

    Based on: Alemu et al. (2025), J. Geotech. Geoenviron. Eng., 151(9): 04025091.
    Valid for: 0 ≤ PI ≤ 32, 10 ≤ p′ ≤ 125 kPa, 1 ≤ OCR ≤ 9.1.

    Parameters
    ----------
    unit_wt : float, default 0
        Unit weight [kN/m³].
    name : str, optional
        Identification label.
    plas_index : float, default 0
        Plasticity index [percent].
    ocr : float, default 1
        Over-consolidation ratio.
    stress_mean : float, default 101.3
        Mean effective stress [kN/m²].
    fines_cont : float, default 1.0
        Fines content [decimal]; typical range 0.5–1.0.
    strains : array_like or None
        Shear strain levels [decimal]; defaults to ``np.logspace(-6, -1.5, 20)``.
    """

    # Table 2 — G/Gmax backbone
    _B1, _B2, _B3, _B4, _B5 = 0.166, -1.678, 0.196, 0.025, 0.675
    _A, _B = 0.905, 1.373
    # Table 3 — minimum damping
    _D1, _D2 = 3.64e-4, 0.012
    # Table 4 — strain-dependent damping
    _E1, _E2, _E3, _E4, _E5, _E6 = 29.805, -0.381, 0.853, 0.021, 1.249, 0.680

    @convert_units(
        unit_wt="kilonewton / meter ** 3",
        stress_mean="kilopascal",
        strains="dimensionless",
    )
    def __init__(
        self,
        unit_wt: float = 0.0,
        name: str = "",
        plas_index: float = 0.0,
        ocr: float = 1.0,
        stress_mean: float = 101.3,
        fines_cont: float = 1.0,
        strains: npt.ArrayLike | None = None,
    ) -> None:
        self._unit_wt = unit_wt
        self._plas_index = plas_index
        self._ocr = ocr
        self._stress_mean = stress_mean
        self._fines_cont = fines_cont

        if strains is None:
            strains = np.logspace(-6, -1.5, num=20)
        else:
            strains = np.asarray(strains, dtype=float)

        if not name:
            name = (
                f"Alemu et al. (2025) - "
                f"PI={plas_index:.0f}, OCR={ocr:.1f}, p'={stress_mean:.0f} kPa"
            )
        self.name = name

        stress_ratio = stress_mean * _KPA_TO_ATM

        # G/Gmax backbone (Eqs. 2, 11)
        gamma_mr = (
            self._B1 * ocr**self._B2 * stress_ratio**self._B3
            + self._B4 * (plas_index + 1) ** self._B5
        ) / 100
        mod_reduc = 1.0 / (1.0 + (strains / gamma_mr) ** self._A) ** self._B

        # Minimum damping (Eq. 16)
        d_min = (
            self._D1 * ocr * (plas_index + 1) + self._D2 * fines_cont
        ) / stress_ratio

        # Strain-dependent damping (Eq. 17)
        gamma_d = (
            ocr**self._E2 * stress_ratio**self._E3 + self._E4 * plas_index
        ) ** self._E5 / 100
        ratio = strains / gamma_d
        damping = d_min + (self._E1 / 100) * ratio / (1.0 + ratio) ** self._E6

        self._strains = strains
        self._mod_reduc = mod_reduc
        self._damping = damping
        self._damping_min = float(d_min)

    def curves(self) -> NonlinearSoilCurves:
        return NonlinearSoilCurves(
            strains=self._strains,
            mod_reduc=self._mod_reduc,
            damping=self._damping,
            damping_min=self._damping_min,
            unit_wt=self._unit_wt,
            name=self.name,
        )
