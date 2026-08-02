"""Scalar ground-motion intensity measures conditioned on other measures.

These models predict a single intensity measure -- Arias intensity, CAV,
PGV -- either conditioned on another measure such as PGA or Sa(1s), or in
scenario mode from a backbone ground motion model.
"""

from .abrahamson_bhasin_2020 import AbrahamsonBhasin2020
from .abrahamson_shi_yang_2016 import AbrahamsonShiYang2016
from .macedo_abrahamson_liu_2021 import MacedoAbrahamsonLiu2021

__all__ = [
    "AbrahamsonBhasin2020",
    "AbrahamsonShiYang2016",
    "MacedoAbrahamsonLiu2021",
]
