"""Models whose independent variable is **frequency**.

Fourier-amplitude-spectrum predictions and the inter-frequency correlation
models that describe their variability.
"""

from .bayless_abrahamson_2018 import BaylessAbrahamson2018
from .bayless_abrahamson_2019 import BaylessAbrahamson2019
from .source_theory import SourceTheoryModel
from .stafford_2017 import Stafford2017
from .stafford_2022 import StaffordEtAl2022

__all__ = [
    "BaylessAbrahamson2018",
    "BaylessAbrahamson2019",
    "SourceTheoryModel",
    "Stafford2017",
    "StaffordEtAl2022",
]
