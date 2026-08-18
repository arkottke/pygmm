"""pyGMM: Ground motion models implemented in Python."""

import importlib
import logging
import warnings

try:
    from ._version import __version__
except ImportError:
    # For development installs
    __version__ = "unknown"

# Expose the model subpackages as attributes, so that `import pygmm` followed
# by `pygmm.fourier_spectrum.SourceTheoryModel` works. Without this they are
# reachable only via an explicit `import pygmm.fourier_spectrum`, and which
# ones happened to be bound depended on whether some name was imported from
# them below.
#
# There are three governed by the independent variable: shear strain ->
# soil_curves, frequency -> fourier_spectrum, period (or a scalar conditioned
# on a Scenario) -> ground_motion. Anything else stays at the top level until
# a second producer of the same quantity earns it a package -- velocity_profile
# is that case: multiple producers of contracts.VelocityProfile.
#
# `soil_curves` is the exception: it is imported lazily (see `__getattr__`
# below) because it is the only package pulling a heavy third-party
# dependency, `pint`. Importing it eagerly is how `import pygmm` came to fail
# outright on any install without pint.
from . import (
    contracts,
    fourier_spectrum,
    ground_motion,
    registry,
    velocity_profile,
)
from .fourier_spectrum import (
    BaylessAbrahamson2018,
    BaylessAbrahamson2019,
    Stafford2017,
)
from .ground_motion import (
    AbrahamsonBhasin2020,
    AbrahamsonGregorAddo2016,
    AbrahamsonShiYang2016,
    AbrahamsonSilva1996,
    AbrahamsonSilvaKamai2014,
    AfshariStewart2016,
    AkkarSandikkayaBommer2014,
    AtkinsonBoore2006,
    BooreStewartSeyhanAtkinson2014,
    Campbell2003,
    CampbellBozorgnia2014,
    ChiouYoungs2014,
    CoppersmithBommer2014,
    DerrasBardCotton2014,
    GulerceAbrahamson2011,
    Idriss2014,
    KemptonStewart2006,
    MacedoAbrahamsonLiu2021,
    PezeshkZandiehTavakoli2011,
    PinillaRamosEtAl2023,
    PinillaRamosEtAl2024,
    TavakoliPezeshk05,
)
from .model import Scenario
from .registry import ModelInfo, get_model, register
from .velocity_profile import (
    bj97gr760_profile,
    bj97gr_profile,
    bj97gvhr_profile,
    btc11_profile,
    kea16_profile,
    sa18_profile,
)

__all__ = [
    # Model subpackages
    "contracts",
    "fourier_spectrum",
    "ground_motion",
    "soil_curves",
    "velocity_profile",
    # Registry
    "ModelInfo",
    "find_models",
    "get_model",
    "register",
    # Models and helpers
    "Scenario",
    "AbrahamsonBhasin2020",
    "AbrahamsonShiYang2016",
    "AbrahamsonSilva1996",
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
    "sa18_profile",
    "bj97gr760_profile",
    "bj97gr_profile",
    "bj97gvhr_profile",
    "btc11_profile",
]

__author__ = "Albert Kottke"
__copyright__ = "Copyright 2016 Albert Kottke"
__license__ = "MIT"
__title__ = "pyGMM"

# Set default logging handler to avoid "No handler found" warnings.
logging.getLogger(__name__).addHandler(logging.NullHandler())


#: Names served from `soil_curves`, which is imported on first use rather than
#: at `import pygmm` -- see the note on `pint` above.
_LAZY_SOIL_CURVES = frozenset(
    {
        "AlemuEtAlSoilType",
        "DarendeliSoilType",
        "KishidaSoilType",
        "MenqSoilType",
        "RollinsEtAlSoilType",
        "WangSoilType",
    }
)


def _load_soil_curves():
    # `from . import soil_curves` would re-enter `__getattr__` below and
    # recurse; `import_module` binds the submodule without an attribute
    # lookup on this package.
    return importlib.import_module(".soil_curves", __name__)


def find_models(**kwds):
    """Return registered models matching the filters.

    Wraps :func:`pygmm.registry.find_models`, importing the lazily-loaded
    ``soil_curves`` package first so the registry is always complete no matter
    what the caller has imported.

    See :func:`pygmm.registry.find_models` for the accepted filters.
    """
    _load_soil_curves()
    return registry.find_models(**kwds)


def __getattr__(name):
    if name in _LAZY_SOIL_CURVES:
        return getattr(_load_soil_curves(), name)

    if name == "soil_curves":
        return _load_soil_curves()

    if name == "models":
        warnings.warn(
            "pygmm.models is deprecated; use pygmm.find_models(provides='psa') "
            "or another capability filter. The list mixed response-spectrum, "
            "duration, FAS and correlation models with no way to tell them "
            "apart -- 7 of its 20 entries had no .spec_accels.",
            DeprecationWarning,
            stacklevel=2,
        )
        return [info.cls for info in find_models()]

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
