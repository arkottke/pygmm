"""Wang and Stokoe (2022) empirical nonlinear model."""

from __future__ import annotations

from functools import wraps

import numpy as np
import numpy.typing as npt

from ..contracts import NonlinearSoilCurves
from ._base import SoilCurveModel
from ._units import convert_kwds_units, convert_units

_KPA_TO_ATM = 1.0 / 101.325


def _to_decimal(*keys):
    """Divide the named kwargs by 100 (percent to decimal) before calling.

    Applied to the keyword arguments named in *keys*.
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(
                *args, **{k: v / 100 if k in keys else v for k, v in kwargs.items()}
            )

        return wrapper

    return decorator


class WangSoilType(SoilCurveModel):
    """Wang and Stokoe (2022) empirical nonlinear model for soils.

    Parameters
    ----------
    soil_group : str
        One of ``'clean_sand_and_gravel'``, ``'nonplastic_silty_sand'``,
        or ``'clayey_soil'``.
    name : str, optional
        Identification label.
    unit_wt : float, default 0
        Unit weight [kN/m³].
    damping_min : float or None
        Minimum damping [decimal]; computed when *None*.
    strains : array_like or None
        Shear strain levels [decimal]; defaults to ``np.logspace(-6, -1.5, 20)``.
    **kwds
        Index properties (``stress_mean`` [kN/m²], ``void_ratio``, ``coef_unif``,
        ``diam_50`` [mm], ``fines_cont`` [decimal], ``plas_index`` [decimal],
        ``water_cont`` [decimal], ``ocr``).
    """

    FACTORS = {
        "clean_sand_and_gravel": [
            "stress_mean",
            "fines_cont",
            "coef_unif",
            "diam_50",
            "water_cont",
            "void_ratio",
            "diam_50",
        ],
        "nonplastic_silty_sand": [
            "stress_mean",
            "void_ratio",
            "fines_cont",
            "water_cont",
        ],
        "clayey_soil": [
            "stress_mean",
            "void_ratio",
            "plas_index",
            "fines_cont",
            "ocr",
            "water_cont",
        ],
    }

    LEVELS = {
        "clean_sand_and_gravel": {
            "gmax_model": [
                "stress_mean",
                "void_ratio",
                "diam_50",
                "coef_unif",
                "fines_cont",
            ],
            "ggmax_model": ["stress_mean", "void_ratio", "coef_unif", "fines_cont"],
            "dmin_model": [
                "stress_mean",
                "fines_cont",
                "water_cont",
                "void_ratio",
                "diam_50",
            ],
            "damping_model": ["stress_mean", "void_ratio", "fines_cont", "coef_unif"],
        },
        "nonplastic_silty_sand": {
            "gmax_model": ["stress_mean", "void_ratio", "water_cont"],
            "ggmax_model": ["stress_mean", "void_ratio", "fines_cont"],
            "dmin_model": ["stress_mean", "void_ratio", "fines_cont"],
            "damping_model": ["stress_mean", "void_ratio", "fines_cont"],
        },
        "clayey_soil": {
            "gmax_model": [
                "stress_mean",
                "void_ratio",
                "ocr",
                "fines_cont",
                "plas_index",
            ],
            "ggmax_model": [
                "stress_mean",
                "void_ratio",
                "fines_cont",
                "ocr",
                "plas_index",
            ],
            "dmin_model": ["stress_mean", "void_ratio", "plas_index", "fines_cont"],
            "damping_model": ["stress_mean", "water_cont", "plas_index", "fines_cont"],
        },
    }

    @convert_units(unit_wt="kilonewton / meter ** 3", strains="dimensionless")
    @convert_kwds_units(stress_mean="kilopascal")
    def __init__(
        self,
        soil_group: str,
        name: str = "",
        unit_wt: float = 0.0,
        damping_min: float | None = None,
        strains: npt.ArrayLike | None = None,
        **kwds,
    ) -> None:
        self._unit_wt = unit_wt
        self._soil_group = soil_group
        if soil_group not in self.LEVELS:
            raise ValueError(
                f"Unknown soil_group {soil_group!r}. Choose from {list(self.LEVELS)}"
            )
        self._index_params = {
            p: kwds[p]
            for p in self.params(soil_group, damping_min is not None)
            if p in kwds
        }

        if strains is None:
            strains = np.logspace(-6, -1.5, num=20)
        else:
            strains = np.asarray(strains, dtype=float)

        if not name:
            name = self._create_name()
        self.name = name

        if "stress_mean" not in self._index_params:
            raise ValueError("`stress_mean` is required for all calculations")

        mod_reduc = self.calc_mod_reduc(strains, soil_group, **self._index_params)

        if damping_min is None:
            damping_min = self.calc_damping_min(soil_group, **self._index_params)

        damping = self.calc_damping(
            strains, soil_group, damping_min, **self._index_params
        )

        self._strains = strains
        self._mod_reduc = mod_reduc
        self._damping = damping
        self._damping_min = float(damping_min)

    @classmethod
    def params(cls, soil_group: str, specified_dmin: bool) -> list[str]:
        models = ["ggmax_model", "damping_model"]
        if not specified_dmin:
            models += ["dmin_model"]
        return list({p for m in models for p in cls.LEVELS[soil_group][m]})

    def _create_name(self) -> str:
        fmt = {
            "coef_unif": "Cᵤ={:.1f}",
            "diam_50": "D₅₀={:.1f} mm",
            "fines_cont": "FC={:.0f} %",
            "ocr": "OCR={:.1f}",
            "plas_index": "PI={:.0f}",
            "stress_mean": "σₘ'={:.1f} kN/m²",
            "void_ratio": "e={:0.2f}",
            "water_cont": "w_c={:.1f}%",
        }
        parts = [fmt[k].format(v) for k, v in self._index_params.items() if k in fmt]
        return "Wang & Stokoe ({})".format(", ".join(parts))

    @property
    def soil_group(self) -> str:
        return self._soil_group

    @property
    def index_params(self) -> dict:
        return self._index_params

    @classmethod
    def get_level(cls, model: str, soil_group: str, **kwds) -> int:
        required = cls.LEVELS[soil_group][model]
        provided = list(kwds.keys())
        lvl = -1
        for req in required:
            if req in provided:
                lvl += 1
            else:
                break
        return lvl

    @classmethod
    @_to_decimal("fines_cont", "plas_index", "water_cont")
    def calc_mod_reduc(
        cls, strains: npt.ArrayLike, soil_group: str, **kwds
    ) -> np.ndarray:
        level = cls.get_level("ggmax_model", soil_group, **kwds)
        if soil_group == "clean_sand_and_gravel":
            if level == 0:
                a, b = 0.729, 0.985
                gamma_mr = 0.068 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.4
            elif level == 1:
                a, b = 0.804, 0.882
                gamma_mr = (0.13 * kwds["void_ratio"] ** 0.545 - 0.043) * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                ) ** 0.45
            elif level == 2:
                a, b = 0.834, 0.839
                gamma_mr = (
                    0.05 * kwds["void_ratio"] ** (0.1 * kwds["coef_unif"]) + 0.011
                ) * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.45
            else:
                a = 0.834 + kwds["fines_cont"]
                b = 0.844 - 1.897 * kwds["fines_cont"]
                gamma_mr = (
                    0.048 * kwds["void_ratio"] ** (0.089 * kwds["coef_unif"]) + 0.008
                ) * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.4
        elif soil_group == "nonplastic_silty_sand":
            if level == 0:
                a = 1.04
                b = 0.438 - 0.007 * kwds["stress_mean"] * _KPA_TO_ATM
                gamma_mr = 0.011 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.318
            elif level == 1:
                a = 1.139 * np.exp(0.093 * kwds["void_ratio"])
                b = 0.475 - 0.007 * kwds["stress_mean"] * _KPA_TO_ATM
                gamma_mr = (0.029 * kwds["void_ratio"] - 0.003) * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                ) ** 0.335
            else:
                a = (1.495 * kwds["void_ratio"] + 3.079 * kwds["fines_cont"]) ** 0.121
                b = 0.486 - 0.006 * kwds["stress_mean"] * _KPA_TO_ATM
                gamma_mr = (0.031 * kwds["void_ratio"] - 0.003) * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                ) ** (0.405 - 0.193 * kwds["fines_cont"])
        elif soil_group == "clayey_soil":
            if level == 0:
                a, b = 1.364, 0.28
                gamma_mr = 0.015 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.205
            elif level == 1:
                a, b = 1.185, 0.475
                gamma_mr = (
                    0.035
                    * kwds["void_ratio"]
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.276
                )
            elif level == 2:
                a = 0.966 + 0.378 * kwds["fines_cont"]
                b = 0.596 - 0.207 * kwds["fines_cont"]
                gamma_mr = (0.031 * kwds["void_ratio"] + 0.004 * kwds["fines_cont"]) * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                ) ** 0.25
            elif level == 3:
                a = 0.972 + 0.419 * kwds["fines_cont"]
                b = 0.571 - 0.2 * kwds["fines_cont"]
                gamma_mr = (
                    0.025 * kwds["void_ratio"] + 0.0015 * kwds["fines_cont"]
                ) * (kwds["stress_mean"] * _KPA_TO_ATM + 0.375 * kwds["ocr"]) ** 0.358
            else:
                a = 0.896 + 0.412 * kwds["fines_cont"] + 0.534 * kwds["plas_index"]
                b = 0.586 - 0.098 * kwds["void_ratio"] - 0.135 * kwds["fines_cont"]
                gamma_mr = (0.02 * kwds["void_ratio"] + 0.004 * kwds["fines_cont"]) * (
                    kwds["stress_mean"] * _KPA_TO_ATM + 0.42 * kwds["ocr"]
                ) ** (0.447 - 0.27 * kwds["plas_index"])
        else:
            raise ValueError("Invalid soil group")
        return 1 / (1 + (100 * np.asarray(strains) / gamma_mr) ** a) ** b

    @classmethod
    @_to_decimal("fines_cont", "plas_index", "water_cont")
    def calc_damping_min(cls, soil_group: str, **kwds) -> float:
        level = cls.get_level("dmin_model", soil_group, **kwds)
        if soil_group == "clean_sand_and_gravel":
            if level == 0:
                d_min = 0.77 * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.03
            elif level == 1:
                d_min = (
                    0.55
                    * (1 + 29.03 * kwds["fines_cont"])
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.13
                )
            elif level == 2:
                d_min = (
                    0.64
                    * (0.26 - kwds["water_cont"]) ** 0.11
                    * (1 + 32.8 * kwds["fines_cont"])
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.12
                )
            elif level == 3:
                d_min = (
                    0.55
                    * (1 - kwds["water_cont"]) ** (-12.49 + 19.45 * kwds["void_ratio"])
                    * (1 + 23.44 * kwds["fines_cont"])
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.14
                )
            else:
                d_min = (
                    0.6
                    * (0.99 + kwds["water_cont"])
                    ** (7.45 - 15.23 * kwds["void_ratio"] + 4.29 * kwds["diam_50"])
                    * (1 + 21.17 * kwds["fines_cont"])
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.14
                )
        elif soil_group == "nonplastic_silty_sand":
            if level == 0:
                d_min = 1.47 * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.2
            elif level == 1:
                d_min = (
                    39.11
                    * (0.44 * kwds["void_ratio"]) ** (4.32 * kwds["void_ratio"])
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.19
                )
            else:
                d_min = (
                    52.16
                    * (0.41 * kwds["void_ratio"])
                    ** (0.81 * kwds["fines_cont"] + 5.2 * kwds["void_ratio"])
                    * (1 + 5.35 * kwds["fines_cont"])
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.19
                )
        elif soil_group == "clayey_soil":
            if level == 0:
                d_min = 2.55 * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.11
            elif level == 1:
                d_min = (
                    7.62
                    * 13.48 ** -kwds["void_ratio"]
                    * (kwds["stress_mean"] * _KPA_TO_ATM) ** -0.29
                    + 0.72 ** -kwds["void_ratio"]
                )
            elif level == 2:
                d_min = 7.29 * 8 ** (
                    -kwds["void_ratio"] - 3.31 * kwds["plas_index"]
                ) * (1 + 148 * kwds["plas_index"] ** 1.95) * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                ) ** -0.2 + (0.5 * kwds["plas_index"]) ** (
                    2.54 - 1.8 * kwds["void_ratio"]
                )
            else:
                d_min = 4.86 * (1.99 + kwds["fines_cont"]) ** (
                    -1.91 * kwds["void_ratio"] - 6.5 * kwds["plas_index"]
                ) * (1 + 106.75 * kwds["plas_index"] ** 1.64) * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                ) ** -0.19 + (0.46 * kwds["plas_index"]) ** (
                    1.73 - 1.34 * kwds["void_ratio"]
                )
        else:
            raise ValueError("Invalid soil group")
        return d_min / 100

    @classmethod
    @_to_decimal("fines_cont", "plas_index", "water_cont")
    def calc_damping(cls, strains, soil_group, damping_min=None, **kwds) -> np.ndarray:
        level = cls.get_level("damping_model", soil_group, **kwds)
        if soil_group == "clean_sand_and_gravel":
            if level == 0:
                c, d = 0.93, 15.64
                gamma_d = 0.09 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.32
            elif level == 1:
                c = 1.08 * np.exp(0.62 - 0.73 * kwds["void_ratio"])
                d = 16.39
                gamma_d = 0.09 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.39
            elif level == 2:
                c = 1.02 * np.exp(0.56 - 0.72 * kwds["void_ratio"])
                d = 21.17
                gamma_d = 0.13 * (
                    kwds["stress_mean"] * _KPA_TO_ATM + 17.94 * kwds["fines_cont"]
                ) ** (0.45 - kwds["fines_cont"])
            else:
                c = 0.93 * np.exp(0.34 - 0.8 * kwds["void_ratio"])
                d = 18.13
                gamma_d = (
                    0.13
                    * kwds["coef_unif"] ** -0.31
                    * (kwds["stress_mean"] * _KPA_TO_ATM + 22.04 * kwds["fines_cont"])
                    ** (0.47 - kwds["fines_cont"])
                )
        elif soil_group == "nonplastic_silty_sand":
            if level == 0:
                c, d = 1.187, 13.125
                gamma_d = 0.045 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.293
            elif level == 1:
                c = 1.38 * np.exp(0.25 * kwds["void_ratio"])
                d = 12.09
                gamma_d = (
                    0.0066
                    * (kwds["stress_mean"] * _KPA_TO_ATM + 5.79 * kwds["void_ratio"])
                    ** 1.01
                )
            else:
                c = 1.39 * np.exp(0.27 * kwds["void_ratio"])
                d = 12.13
                gamma_d = 0.0025 * (
                    kwds["stress_mean"] * _KPA_TO_ATM
                    + 5.73 * kwds["void_ratio"]
                    + 9.17 * kwds["fines_cont"]
                ) ** (1.47 - 0.52 * kwds["fines_cont"])
        elif soil_group == "clayey_soil":
            if level == 0:
                c, d = 1.12, 19.47
                gamma_d = 0.11 * (kwds["stress_mean"] * _KPA_TO_ATM) ** 0.23
            elif level == 1:
                c, d = 1.36, 15.16
                gamma_d = 0.29 * (
                    0.017 * kwds["stress_mean"] * _KPA_TO_ATM + kwds["water_cont"]
                ) ** (1.15 + kwds["water_cont"])
            elif level == 2:
                c = 1.48 ** (0.53 + kwds["plas_index"])
                d = 15.61
                gamma_d = 0.07 * (
                    0.06 * kwds["stress_mean"] * _KPA_TO_ATM + 2.69 * kwds["water_cont"]
                ) ** (1.06 + kwds["water_cont"] - kwds["plas_index"])
            else:
                c = (1.91 * kwds["fines_cont"]) ** (1.62 * kwds["plas_index"])
                d = 21.7
                gamma_d = 0.11 * (
                    0.12 * kwds["stress_mean"] * _KPA_TO_ATM
                    + 5.29 * kwds["water_cont"]
                    - kwds["fines_cont"]
                ) ** (
                    1.45
                    - kwds["plas_index"]
                    + kwds["water_cont"]
                    - 1.09 * kwds["fines_cont"]
                )
        else:
            raise ValueError("Invalid soil group")

        if damping_min is None:
            damping_min = cls.calc_damping_min(soil_group, **kwds)

        gamma_ratio = (100 * strains / gamma_d) ** c
        return (d * gamma_ratio + 100 * damping_min) / (gamma_ratio + 1) / 100

    def curves(self) -> NonlinearSoilCurves:
        return NonlinearSoilCurves(
            strains=self._strains,
            mod_reduc=self._mod_reduc,
            damping=self._damping,
            damping_min=self._damping_min,
            unit_wt=self._unit_wt,
            name=self.name,
        )
