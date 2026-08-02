"""pyGMM: Ground motion models implemented in Python."""

import logging

try:
    from ._version import __version__
except ImportError:
    # For development installs
    __version__ = "unknown"

# Expose the category subpackages as attributes, so that `import pygmm`
# followed by `pygmm.fourier_spectrum.SourceTheoryModel` works. Without this
# they are reachable only via an explicit `import pygmm.fourier_spectrum`, and
# which ones happened to be bound depended on whether some name was imported
# from them below.
from . import (
    contracts,
    correlation,
    cpt,
    duration,
    fault_displacement,
    fourier_spectrum,
    gm_intensity,
    response_spectrum,
    soil_curves,
    tools,
    velocity_profile,
)
from .bayless_abrahamson_2018 import BaylessAbrahamson2018
from .boore_stewart_seyhan_atkinson_2014 import BooreStewartSeyhanAtkinson2014
from .campbell_bozorgnia_2014 import CampbellBozorgnia2014
from .chiou_youngs_2014 import ChiouYoungs2014
from .derras_bard_cotton_2014 import DerrasBardCotton2014
from .duration import (
    AbrahamsonSilva1996,
    AfshariStewart2016,
    KemptonStewart2006,
    PinillaRamosEtAl2023,
    PinillaRamosEtAl2024,
)
from .fourier_spectrum import BaylessAbrahamson2019
from .gm_intensity import (
    AbrahamsonBhasin2020,
    AbrahamsonShiYang2016,
    MacedoAbrahamsonLiu2021,
)
from .gulerce_abrahamson_2011 import GulerceAbrahamson2011
from .model import Scenario
from .response_spectrum import (
    AbrahamsonGregorAddo2016,
    AbrahamsonSilvaKamai2014,
    AkkarSandikkayaBommer2014,
    AtkinsonBoore2006,
    Campbell2003,
    CoppersmithBommer2014,
    Idriss2014,
    PezeshkZandiehTavakoli2011,
    TavakoliPezeshk05,
)
from .soil_curves import (
    AlemuEtAlSoilType,
    DarendeliSoilType,
    KishidaSoilType,
    MenqSoilType,
    RollinsEtAlSoilType,
    WangSoilType,
)
from .stafford_2017 import Stafford2017
from .velocity_profile import kea16_profile

__all__ = [
    # Category subpackages
    "contracts",
    "correlation",
    "cpt",
    "duration",
    "fault_displacement",
    "fourier_spectrum",
    "gm_intensity",
    "response_spectrum",
    "soil_curves",
    "tools",
    "velocity_profile",
    # Models and helpers
    "Scenario",
    "AbrahamsonBhasin2020",
    "AbrahamsonShiYang2016",
    "AbrahamsonSilvaKamai2014",
    "AbrahamsonGregorAddo2016",
    "AfshariStewart2016",
    "AkkarSandikkayaBommer2014",
    "AtkinsonBoore2006",
    "BaylessAbrahamson2018",
    "BaylessAbrahamson2019",
    "BooreStewartSeyhanAtkinson2014",
    "Campbell2003",
    "CampbellBozorgnia2014",
    "ChiouYoungs2014",
    "CoppersmithBommer2014",
    "DerrasBardCotton2014",
    "GulerceAbrahamson2011",
    "KemptonStewart2006",
    "Idriss2014",
    "MacedoAbrahamsonLiu2021",
    "PezeshkZandiehTavakoli2011",
    "PinillaRamosEtAl2023",
    "PinillaRamosEtAl2024",
    "TavakoliPezeshk05",
    "Stafford2017",
    "AlemuEtAlSoilType",
    "DarendeliSoilType",
    "KishidaSoilType",
    "MenqSoilType",
    "RollinsEtAlSoilType",
    "WangSoilType",
    "kea16_profile",
]

__author__ = "Albert Kottke"
__copyright__ = "Copyright 2016 Albert Kottke"
__license__ = "MIT"
__title__ = "pyGMM"

# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:

    class NullHandler(logging.Handler):
        def emit(self, record):
            pass


logging.getLogger(__name__).addHandler(NullHandler())

models = [
    AbrahamsonSilva1996,
    AbrahamsonSilvaKamai2014,
    AfshariStewart2016,
    AkkarSandikkayaBommer2014,
    AtkinsonBoore2006,
    BaylessAbrahamson2019,
    BooreStewartSeyhanAtkinson2014,
    Campbell2003,
    CampbellBozorgnia2014,
    ChiouYoungs2014,
    CoppersmithBommer2014,
    DerrasBardCotton2014,
    GulerceAbrahamson2011,
    KemptonStewart2006,
    Idriss2014,
    PezeshkZandiehTavakoli2011,
    PinillaRamosEtAl2023,
    PinillaRamosEtAl2024,
    TavakoliPezeshk05,
    Stafford2017,
]
