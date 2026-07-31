"""Kishida (2017) nonlinear model for highly organic soils."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from ..contracts import NonlinearSoilCurves
from ._base import SoilCurveModel
from ._units import convert_units

# Standard gravity [m/s²] — used to convert density to unit weight
_GRAVITY = 9.80665


class KishidaSoilType(SoilCurveModel):
    """Kishida (2017) empirical nonlinear model for highly organic soils.

    Parameters
    ----------
    name : str, optional
        Identification label.
    unit_wt : float or None
        Unit weight [kN/m³].  Computed from the empirical model when *None*.
    stress_vert : float, default 101.3
        Vertical effective stress [kN/m²].
    organic_content : float, default 10
        Organic content [percent].
    lab_consol_ratio : float, default 1
        Laboratory consolidation ratio (use 1 for field applications).
    strains : array_like or None
        Shear strain levels [decimal]; defaults to ``np.logspace(-6, -1.5, 20)``.
    """

    @convert_units(
        unit_wt="kilonewton / meter ** 3",
        stress_vert="kilopascal",
        strains="dimensionless",
    )
    def __init__(
        self,
        name: str = "",
        unit_wt: float | None = None,
        stress_vert: float = 101.3,
        organic_content: float = 10,
        lab_consol_ratio: float = 1,
        strains: npt.ArrayLike | None = None,
    ) -> None:
        self._stress_vert = float(stress_vert)
        self._organic_content = float(organic_content)
        self._lab_consol_ratio = float(lab_consol_ratio)

        if strains is None:
            strains = np.logspace(-6, -1.5, num=20)
        else:
            strains = np.asarray(strains, dtype=float)

        strains_pct = strains * 100

        x_1_mean = -2.5
        x_2_mean = 4.0
        x_3_mean = 0.5

        x_3 = 2.0 / (1 + np.exp(self._organic_content / 23))
        strain_ref = self._calc_strain_ref(x_3, x_3_mean)
        strain_ref_pct = strain_ref * 100
        x_1 = np.log(strains_pct + strain_ref_pct)
        x_2 = np.log(self._stress_vert)

        if unit_wt is None:
            self._unit_wt = self._calc_unit_wt(x_2, x_3)
        else:
            self._unit_wt = float(unit_wt)

        ones = np.ones_like(strains)
        x_2_arr = x_2 * ones
        x_3_arr = x_3 * ones

        mod_reduc = self._calc_mod_reduc(
            strains_pct,
            strain_ref_pct,
            x_1,
            x_1_mean,
            x_2_arr,
            x_2_mean,
            x_3_arr,
            x_3_mean,
        )
        damping = self._calc_damping(mod_reduc, x_2_arr, x_2_mean, x_3_arr, x_3_mean)

        if not name:
            name = (
                f"Kishida (σᵥ'={self._stress_vert:.1f} kN/m², "
                f"OC={self._organic_content:.0f} %)"
            )
        self.name = name

        self._strains = strains
        self._mod_reduc = mod_reduc
        self._damping = damping
        self._damping_min = float(damping[0])

    @staticmethod
    def _calc_strain_ref(x_3: float, x_3_mean: float) -> float:
        return np.exp(-1.41 + -0.950 * (x_3 - x_3_mean)) / 100

    def _calc_mod_reduc(
        self, strains, strain_ref, x_1, x_1_mean, x_2, x_2_mean, x_3, x_3_mean
    ) -> np.ndarray:
        ones = np.ones_like(strains)
        x_4 = np.log(self._lab_consol_ratio) * ones
        x = np.c_[
            ones,
            x_1,
            x_2,
            x_3,
            x_4,
            (x_1 - x_1_mean) * (x_2 - x_2_mean),
            (x_1 - x_1_mean) * (x_3 - x_3_mean),
            (x_2 - x_2_mean) * (x_3 - x_3_mean),
            (x_1 - x_1_mean) * (x_2 - x_2_mean) * (x_3 - x_3_mean),
        ]
        denom = np.log(1 / strain_ref + strains / strain_ref)
        b = np.c_[
            5.11 * ones,
            -0.729 * ones,
            (1 - 0.37 * x_3_mean * (1 + ((np.log(strain_ref) - x_1_mean) / denom))),
            -0.693 * ones,
            0.8 - 0.4 * x_3,
            0.37 * x_3_mean / denom,
            0.0 * ones,
            -0.37 * (1 + (np.log(strain_ref) - x_1_mean) / denom),
            0.37 / denom,
        ]
        shear_mod = np.exp((b * x).sum(axis=1))
        return shear_mod / shear_mod[0]

    @staticmethod
    def _calc_damping(mod_reduc, x_2, x_2_mean, x_3, x_3_mean) -> np.ndarray:
        x_1_mean = -1.0
        x_1 = np.log(np.log(1 / mod_reduc) + 0.103)
        ones = np.ones_like(mod_reduc)
        x = np.c_[
            ones,
            x_1,
            x_2,
            x_3,
            (x_1 - x_1_mean) * (x_2 - x_2_mean),
            (x_2 - x_2_mean) * (x_3 - x_3_mean),
        ]
        c = np.c_[2.86, 0.571, -0.103, -0.141, 0.0419, -0.240]
        return np.exp((c * x).sum(axis=1)) / 100.0

    @staticmethod
    def _calc_unit_wt(x_1: float, x_2: float) -> float:
        x = np.r_[1, x_1, x_2]
        d = np.r_[-0.112, 0.038, 0.360]
        return np.exp(d @ x) * _GRAVITY

    def curves(self) -> NonlinearSoilCurves:
        return NonlinearSoilCurves(
            strains=self._strains,
            mod_reduc=self._mod_reduc,
            damping=self._damping,
            damping_min=self._damping_min,
            unit_wt=self._unit_wt,
            name=self.name,
        )
