#!/usr/bin/python
# -*- coding: utf-8 -*-

import logging
import pkg_resources

from .abrahamson_silva_kamai_2014 import AbrahamsonSilvaKamai2014
from .atkinson_boore_2006 import AtkinsonBoore2006
from .akkar_sandikkaya_bommer_2014 import AkkarSandikkayaBommer2014
from .boore_stewart_seyhan_atkinson_2014 import BooreStewartSeyhanAtkinson2014
from .campbell_2003 import Campbell2003
from .campbell_bozorgnia_2014 import CampbellBozorgnia2014
from .chiou_youngs_2014 import ChiouYoungs2014
from .derras_bard_cotton_2014 import DerrasBardCotton2014
from .idriss_2014 import Idriss2014
from .pezeshk_zandieh_tavakoli_2011 import PezeshkZandiehTavakoli2011
from .tavakoli_pezeshk_2005 import TavakoliPezeshk05

__all__ = [
    'AbrahamsonSilvaKamai2014',
    'AkkarSandikkayaBommer2014',
    'AtkinsonBoore2006',
    'BooreStewartSeyhanAtkinson2014',
    'Campbell2003',
    'CampbellBozorgnia2014',
    'ChiouYoungs2014',
    'DerrasBardCotton2014',
    'Idriss2014',
    'PezeshkZandiehTavakoli2011',
    'TavakoliPezeshk05',
]

__author__ = 'Albert Kottke'
__copyright__ = 'Copyright 2016 Albert Kottke'
__license__ = 'MIT'
__title__ = 'pyGMM'
__version__ = pkg_resources.get_distribution('pygmm').version


# Set default logging handler to avoid "No handler found" warnings.
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):
        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

models = [
    AbrahamsonSilvaKamai2014,
    AkkarSandikkayaBommer2014,
    AtkinsonBoore2006,
    BooreStewartSeyhanAtkinson2014,
    Campbell2003,
    CampbellBozorgnia2014,
    ChiouYoungs2014,
    DerrasBardCotton2014,
    Idriss2014,
    PezeshkZandiehTavakoli2011,
    TavakoliPezeshk05,
]
