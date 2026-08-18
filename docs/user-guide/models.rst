=======================
Ground Motion Models
=======================

pyGMM provides a unified interface to numerous ground motion prediction equations (GMPEs)
developed by researchers worldwide. All models follow a consistent API for easy comparison
and analysis.

.. grid:: 1 2 2 2
    :gutter: 3

    .. grid-item-card:: 🏗️ Interface Design
        :link: generic-interface
        :link-type: ref

        All models share a common interface for consistent usage across different equations.

    .. grid-item-card:: 📊 Available Models
        :link: available-models
        :link-type: ref

        Browse the comprehensive list of implemented ground motion models.

    .. grid-item-card:: 🔧 Model Selection
        :link: model-selection-guide
        :link-type: ref

        Guidelines for choosing appropriate models for your analysis.

    .. grid-item-card:: 📝 Usage Examples
        :link: ../examples/index
        :link-type: doc

        Practical examples of using ground motion models.

.. _generic-interface:

Generic Interface
=================

All ground motion models in pyGMM inherit from :class:`~pygmm.model.GroundMotionModel`,
providing a consistent interface regardless of the underlying equation.

.. currentmodule:: pygmm.model

.. autoclass:: GroundMotionModel
   :members:
   :undoc-members:
   :show-inheritance:

   .. rubric:: Key Methods

   .. autosummary::

      ~GroundMotionModel.__init__
      ~GroundMotionModel.interp_ln_spec_accels
      ~GroundMotionModel.interp_spec_accels
      ~GroundMotionModel.interp_ln_stds

   .. rubric:: Properties

   .. autosummary::

      ~GroundMotionModel.periods
      ~GroundMotionModel.spec_accels
      ~GroundMotionModel.ln_stds
      ~GroundMotionModel.pga
      ~GroundMotionModel.pgv
      ~GroundMotionModel.pgd

Basic Usage
-----------

All models follow the same usage pattern:

.. code-block:: python

   import pygmm

   # Create scenario
   scenario = pygmm.Scenario(
       mag=6.5,
       dist_rup=20,
       v_s30=760,
       mechanism='strike_slip'
   )

   # Initialize model
   model = pygmm.CampbellBozorgnia2014(scenario)

   # Calculate ground motion
   ln_sa, ln_std = np.log(model.spec_accels), model.ln_stds

.. _available-models:

Available Models
================

pyGMM includes implementations of the following ground motion prediction equations:

Active Shallow Crust Models
----------------------------

.. grid:: 1 2 3 3
    :gutter: 2

    .. grid-item-card:: ASK14
        :class-title: text-center

        **Abrahamson, Silva & Kamai (2014)**

        NGA-West2 model for active tectonic regions

    .. grid-item-card:: BSSA14
        :class-title: text-center

        **Boore et al. (2014)**

        NGA-West2 BSSA model

    .. grid-item-card:: CB14
        :class-title: text-center

        **Campbell & Bozorgnia (2014)**

        NGA-West2 Campbell-Bozorgnia model

    .. grid-item-card:: CY14
        :class-title: text-center

        **Chiou & Youngs (2014)**

        NGA-West2 Chiou-Youngs model

    .. grid-item-card:: I14
        :class-title: text-center

        **Idriss (2014)**

        NGA-West2 Idriss model

    .. grid-item-card:: BA19
        :class-title: text-center

        **Bayless & Abrahamson (2019)**

        Updated ASK model implementation

Stable Continental Regions
--------------------------

.. grid:: 1 2 3 3
    :gutter: 2

    .. grid-item-card:: AB06
        :class-title: text-center

        **Atkinson & Boore (2006)**

        Eastern North America model

    .. grid-item-card:: PZT11
        :class-title: text-center

        **Pezeshk et al. (2011)**

        Central and Eastern US model

    .. grid-item-card:: TP05
        :class-title: text-center

        **Tavakoli & Pezeshk (2005)**

        Eastern North America model

Specialized Models
------------------

.. grid:: 1 2 2 2
    :gutter: 2

    .. grid-item-card:: Arias Intensity
        :class-title: text-center

        **Abrahamson, Shi & Yang (2016)**

        Arias intensity consistent with NGA-West2

    .. grid-item-card:: Vertical Components
        :class-title: text-center

        **Gulerce & Abrahamson (2011)**

        Vertical-to-horizontal ratios

    .. grid-item-card:: Conditional Spectra
        :class-title: text-center

        **Baker & Jayaram (2008)**

        Conditional mean spectrum

    .. grid-item-card:: Site Response
        :class-title: text-center

        **Kempton & Stewart (2006)**

        Site response modifications

    .. grid-item-card:: Duration
        :class-title: text-center

        **Coppersmith & Bommer (2014)**, **Afshari & Stewart (2016)**,
        **Kempton & Stewart (2006)**, **Pinilla-Ramos et al. (2023, 2024)**

        Significant duration models (crustal and subduction)

.. _detailed-model-list:

Detailed Model List
-------------------

Generated from the model registry. To reproduce::

    python -c "import pygmm; [print(i.key, sorted(i.provides)) for i in pygmm.find_models()]"

.. table:: Registered models and their capabilities
   :widths: auto

   ======  ===========================================  ==================  ==================
   Abbrev  Model                                        Provides            Tectonic setting
   ======  ===========================================  ==================  ==================
   AB20    Abrahamson and Bhasin (2020)                 pgv                 --
   AGA16   Abrahamson, Gregor, & Addo (2016)            pga, psa            subduction
   ASY16   Abrahamson, Shi, & Yang (2016)               arias               --
   AS96    Abrahamson Silva (1996)                      duration            --
   ASK14   Abrahamson, Silva, & Kamai (2014)            pga, pgv, psa       active_crustal
   AS16    Afshari and Stewart (2016)                   duration            --
   ASB14   Akkar, Sandikkaya, & Bommer (2014)           pga, pgv, psa       --
   --      AlemuEtAlSoilType                            soil_curves         --
   AB06    Atkinson and Boore (2006)                    pga, pgd, pgv, psa  stable_continental
   BA18    Bayless and Abrahamson (2018)                correlation         active_crustal
   BA19    Bayless and Abrahamson (2019)                fas                 active_crustal
   BSSA14  Boore, Stewart, Seyhan, and Atkinson (2014)  pga, pgv, psa       active_crustal
   C03     Campbell (2003)                              psa                 stable_continental
   CB14    Campbell & Bozorgnia (2014)                  pga, pgv, psa       active_crustal
   CY14    Chiou and Youngs (2014)                      pga, pgv, psa       active_crustal
   CB14    Coppersmith and Bommer (2014)                pga, psa            subduction
   --      DarendeliSoilType                            soil_curves         --
   DBC13   Derras, Bard & Cotton (2014)                 pga, pgv, psa       --
   GA11    Gülerce & Abrahamson (2011)                  vh_ratio            active_crustal
   HKR14   Hermkes, Kuehn, Riggelsen (2014)             pga, pgv, psa       --
   I14     Idriss (2014)                                pga, psa            active_crustal
   KS06    Kempton Stewart (2006)                       duration            --
   --      KishidaSoilType                              soil_curves         --
   MAL21   Macedo, Abrahamson, & Liu (2021)             cav                 --
   --      MenqSoilType                                 soil_curves         --
   Pea11   Pezeshk et al. (2011)                        pga, psa            stable_continental
   PR23    Pinilla-Ramos et al. (2023)                  duration            --
   PR24    Pinilla-Ramos et al. (2024)                  duration            subduction
   --      RollinsEtAlSoilType                          soil_curves         --
   ST      Single-corner Source Theory                  fas                 --
   PJS17   Stafford (2017)                              correlation         active_crustal
   Sea22   Stafford et al. (2022)                       fas                 --
   TP05    Tavakoli and Pezeshk (2005)                  pga, psa            stable_continental
   --      WangSoilType                                 soil_curves         --
   ======  ===========================================  ==================  ==================

Models are filed by the independent variable of what they predict -- shear
strain in :mod:`pygmm.soil_curves`, frequency in
:mod:`pygmm.fourier_spectrum`, period (or a scalar conditioned on a
:class:`~pygmm.model.Scenario`) in :mod:`pygmm.ground_motion`. A model's
*capabilities* are registry metadata rather than its location, because one
:class:`~pygmm.model.GroundMotionModel` commonly provides several.

.. currentmodule:: pygmm

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   ground_motion.abrahamson_bhasin_2020.AbrahamsonBhasin2020
   ground_motion.abrahamson_gregor_addo_2016.AbrahamsonGregorAddo2016
   ground_motion.abrahamson_shi_yang_2016.AbrahamsonShiYang2016
   ground_motion.abrahamson_silva_1996.AbrahamsonSilva1996
   ground_motion.abrahamson_silva_kamai_2014.AbrahamsonSilvaKamai2014
   ground_motion.afshari_stewart_2016.AfshariStewart2016
   ground_motion.akkar_sandikkaya_bommer_2014.AkkarSandikkayaBommer2014
   soil_curves.alemu_2025.AlemuEtAlSoilType
   ground_motion.atkinson_boore_2006.AtkinsonBoore2006
   fourier_spectrum.bayless_abrahamson_2018.BaylessAbrahamson2018
   fourier_spectrum.bayless_abrahamson_2019.BaylessAbrahamson2019
   ground_motion.boore_stewart_seyhan_atkinson_2014.BooreStewartSeyhanAtkinson2014
   ground_motion.campbell_2003.Campbell2003
   ground_motion.campbell_bozorgnia_2014.CampbellBozorgnia2014
   ground_motion.chiou_youngs_2014.ChiouYoungs2014
   ground_motion.coppersmith_bommer_2014.CoppersmithBommer2014
   soil_curves.darendeli_2001.DarendeliSoilType
   ground_motion.derras_bard_cotton_2014.DerrasBardCotton2014
   ground_motion.gulerce_abrahamson_2011.GulerceAbrahamson2011
   ground_motion.hermkes_kuehn_riggelsen_2014.HermkesKuehnRiggelsen2014
   ground_motion.idriss_2014.Idriss2014
   ground_motion.kempton_stewart_2006.KemptonStewart2006
   soil_curves.kishida.KishidaSoilType
   ground_motion.macedo_abrahamson_liu_2021.MacedoAbrahamsonLiu2021
   soil_curves.menq.MenqSoilType
   ground_motion.pezeshk_zandieh_tavakoli_2011.PezeshkZandiehTavakoli2011
   ground_motion.pinilla_ramos_et_al_2023.PinillaRamosEtAl2023
   ground_motion.pinilla_ramos_et_al_2024.PinillaRamosEtAl2024
   soil_curves.rollins_2020.RollinsEtAlSoilType
   fourier_spectrum.source_theory.SourceTheoryModel
   fourier_spectrum.stafford_2017.Stafford2017
   fourier_spectrum.stafford_2022.StaffordEtAl2022
   ground_motion.tavakoli_pezeshk_2005.TavakoliPezeshk05
   soil_curves.wang_stokoe_2022.WangSoilType

.. _model-selection-guide:

Model Selection Guide
=====================

Choosing the appropriate ground motion model depends on several factors:

.. tab-set::

    .. tab-item:: Tectonic Setting

        **Active Shallow Crust**

        - California, Japan, Turkey, New Zealand
        - Use NGA-West2 models (ASK14, BSSA14, CB14, CY14, I14)

        **Stable Continental Regions**

        - Eastern North America, Australia, Europe
        - Use AB06, PZT11, or TP05

    .. tab-item:: Magnitude Range

        **Small to Moderate (M < 6.5)**

        - Most models applicable
        - Consider regional calibration

        **Large Events (M > 7.0)**

        - NGA-West2 models well-constrained
        - Check applicable magnitude ranges

    .. tab-item:: Distance Range

        **Near-field (< 20 km)**

        - All NGA-West2 models
        - Consider directivity effects

        **Far-field (> 100 km)**

        - Check distance limits
        - Consider attenuation characteristics

Applicability Ranges
---------------------

Each model has specific ranges of applicability:

.. admonition:: Important
   :class: warning

   Using models outside their intended ranges may produce unreliable results.
   Always check the model documentation for applicable ranges.

.. dropdown:: Check Model Limits
   :class-title: sd-bg-info sd-text-white

   .. code-block:: python

      import pygmm

      model = pygmm.CampbellBozorgnia2014(scenario)
      print("Model limits:")
      for param, limits in model.LIMITS.items():
          print(f"  {param}: {limits}")

Multi-Model Analysis
====================

For robust analyses, consider using multiple models:

.. code-block:: python

   models = [
       pygmm.CampbellBozorgnia2014(scenario),
       pygmm.AbrahamsonSilvaKamai2014(scenario),
       pygmm.BooreStewartSeyhanAtkinson2014(scenario),
       pygmm.ChiouYoungs2014(scenario),
   ]

   # Calculate center, body, and range (CBR) statistics
   results = []
   for model in models:
       ln_sa, ln_std = np.log(model.spec_accels), model.ln_stds
       results.append(np.exp(ln_sa))

   # Central tendency and epistemic uncertainty
   median = np.median(results, axis=0)
   std_log = np.std(np.log(results), axis=0)

See Also
========

- :doc:`examples/model_comparison` - Detailed comparison examples
- :doc:`examples/basic_usage` - Getting started with models
- :doc:`modules` - Complete API reference

Mechanism Reference
==================

The following abbreviations are used for fault mechanism. Refer to each model
for the specific definition of the mechanism.

+--------------+--------------+
| Abbreviation | Name         |
+==============+==============+
| U            | Unspecified  |
+--------------+--------------+
| SS           | Strike-slip  |
+--------------+--------------+
| NS           | Normal slip  |
+--------------+--------------+
| RS           | Reverse slip |
+--------------+--------------+

Specific Models
---------------

Each supported ground motion model inherits from :class:`.Model`, which
provides the standard interface to access the calculated ground motion. The
following models have been implemented.

.. currentmodule:: pygmm
.. autosummary::
    :toctree: _autosummary
    :nosignatures:

    ~abrahamson_gregor_addo_2016.AbrahamsonGregorAddo2016
    ~abrahamson_shi_yang_2016.AbrahamsonShiYang2016
    ~abrahamson_silva_kamai_2014.AbrahamsonSilvaKamai2014
    ~akkar_sandikkaya_bommer_2014.AkkarSandikkayaBommer2014
    ~atkinson_boore_2006.AtkinsonBoore2006
    ~boore_stewart_seyhan_atkinson_2014.BooreStewartSeyhanAtkinson2014
    ~campbell_2003.Campbell2003
    ~campbell_bozorgnia_2014.CampbellBozorgnia2014
    ~chiou_youngs_2014.ChiouYoungs2014
    ~derras_bard_cotton_2014.DerrasBardCotton2014
    ~hermkes_kuehn_riggelsen_2014.HermkesKuehnRiggelsen2014
    ~idriss_2014.Idriss2014
    ~macedo_abrahamson_liu_2021.MacedoAbrahamsonLiu2021
    ~pezeshk_zandieh_tavakoli_2011.PezeshkZandiehTavakoli2011
    ~tavakoli_pezeshk_2005.TavakoliPezeshk05

If you are interested in contributing another model to the collection please see
:doc:`contributing`.

Conditional Spectrum Models
---------------------------

Conditional spectra models are used to create an acceleration response
spectrum conditioned on the response at one or multiple spectral periods. The
The :func:`~pygmm.baker_jayaram_2008.calc_cond_mean_spectrum`
provides functions for developing conditional spectra based on one conditioning
period, while the :func:`~pygmm.kishida_2017.calc_cond_mean_spectrum_vector`
uses the same correlation structure and permits conditioning on multiple
periods.

.. currentmodule:: pygmm
.. autosummary::
    :toctree: _autosummary
    :nosignatures:

    ~baker_jayaram_2008.calc_correls
    ~baker_jayaram_2008.calc_cond_mean_spectrum
    ~kishida_2017.calc_cond_mean_spectrum_vector

Vertical-to-Horizontal (V/H) Models
-----------------------------------

Vertical-to-horizontal models are used to compute the vertical acceleration
response spectrum from a horizontal response spectrum.

.. currentmodule:: pygmm
.. autosummary::
    :toctree: _autosummary
    :nosignatures:

    ~gulerce_abrahamson_2011.GulerceAbrahamson2011
