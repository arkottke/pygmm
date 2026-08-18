"""Models predicting ground-motion intensity from an earthquake scenario.

Everything whose independent variable is spectral **period**, or whose output
is a **scalar** conditioned on a :class:`~pygmm.model.Scenario`: response
spectra, peak parameters, significant duration, Arias intensity, CAV,
vertical-to-horizontal ratios, and inter-period correlations.

Which intensity measures a given model provides is registry metadata, not a
directory: a single :class:`~pygmm.model.GroundMotionModel` commonly emits
PSA, PGA, PGV and PGD from one underlying response vector.

``HermkesKuehnRiggelsen2014`` lives here but is deliberately *not* imported
below: its module attempts to download a 430 kB data file at import time when
that file is absent. Import it directly, which also registers it::

    from pygmm.ground_motion.hermkes_kuehn_riggelsen_2014 import (
        HermkesKuehnRiggelsen2014,
    )
"""

from .abrahamson_bhasin_2020 import AbrahamsonBhasin2020
from .abrahamson_gregor_addo_2016 import AbrahamsonGregorAddo2016
from .abrahamson_shi_yang_2016 import AbrahamsonShiYang2016
from .abrahamson_silva_1996 import AbrahamsonSilva1996
from .abrahamson_silva_kamai_2014 import AbrahamsonSilvaKamai2014
from .afshari_stewart_2016 import AfshariStewart2016
from .akkar_sandikkaya_bommer_2014 import AkkarSandikkayaBommer2014
from .atkinson_boore_2006 import AtkinsonBoore2006
from .baker_jayaram_2008 import calc_cond_mean_spectrum, calc_correls
from .boore_stewart_seyhan_atkinson_2014 import BooreStewartSeyhanAtkinson2014
from .campbell_2003 import Campbell2003
from .campbell_bozorgnia_2014 import CampbellBozorgnia2014
from .chiou_youngs_2014 import ChiouYoungs2014
from .coppersmith_bommer_2014 import CoppersmithBommer2014
from .derras_bard_cotton_2014 import DerrasBardCotton2014
from .gulerce_abrahamson_2011 import GulerceAbrahamson2011
from .idriss_2014 import Idriss2014
from .kempton_stewart_2006 import KemptonStewart2006
from .kishida_2017 import calc_cond_mean_spectrum_vector
from .macedo_abrahamson_liu_2021 import MacedoAbrahamsonLiu2021
from .pezeshk_zandieh_tavakoli_2011 import PezeshkZandiehTavakoli2011
from .pinilla_ramos_et_al_2023 import PinillaRamosEtAl2023
from .pinilla_ramos_et_al_2024 import PinillaRamosEtAl2024
from .tavakoli_pezeshk_2005 import TavakoliPezeshk05

__all__ = [
    "AbrahamsonBhasin2020",
    "AbrahamsonGregorAddo2016",
    "AbrahamsonShiYang2016",
    "AbrahamsonSilva1996",
    "AbrahamsonSilvaKamai2014",
    "AfshariStewart2016",
    "AkkarSandikkayaBommer2014",
    "AtkinsonBoore2006",
    "BooreStewartSeyhanAtkinson2014",
    "Campbell2003",
    "CampbellBozorgnia2014",
    "ChiouYoungs2014",
    "CoppersmithBommer2014",
    "DerrasBardCotton2014",
    "GulerceAbrahamson2011",
    "Idriss2014",
    "KemptonStewart2006",
    "MacedoAbrahamsonLiu2021",
    "PezeshkZandiehTavakoli2011",
    "PinillaRamosEtAl2023",
    "PinillaRamosEtAl2024",
    "TavakoliPezeshk05",
    "calc_cond_mean_spectrum",
    "calc_cond_mean_spectrum_vector",
    "calc_correls",
]
