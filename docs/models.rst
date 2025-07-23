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
   model = pygmm.CampbellBozorgnia2014()

   # Calculate ground motion
   ln_sa, ln_std = model(scenario)

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

        **Coppersmith & Bommer (2014)**

        Significant duration models

Detailed Model List
-------------------

.. currentmodule:: pygmm

.. autosummary::
   :toctree: _autosummary
   :template: class.rst

   abrahamson_gregor_addo_2016.AbrahamsonGregorAddo2016
   abrahamson_silva_kamai_2014.AbrahamsonSilvaKamai2014
   akkar_sandikkaya_bommer_2014.AkkarSandikkayaBommer2014
   atkinson_boore_2006.AtkinsonBoore2006
   boore_stewart_seyhan_atkinson_2014.BooreStewartSeyhanAtkinson2014
   campbell_2003.Campbell2003
   campbell_bozorgnia_2014.CampbellBozorgnia2014
   chiou_youngs_2014.ChiouYoungs2014
   derras_bard_cotton_2014.DerrasBardCotton2014
   gulerce_abrahamson_2011.GulerceAbrahamson2011
   hermkes_kuehn_riggelsen_2014.HermkesKuehnRiggelsen2014
   idriss_2014.Idriss2014
   pezeshk_zandieh_tavakoli_2011.PezeshkZandiehTavakoli2011
   tavakoli_pezeshk_2005.TavakoliPezeshk05

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

      model = pygmm.CampbellBozorgnia2014()
      print("Model limits:")
      for param, limits in model.LIMITS.items():
          print(f"  {param}: {limits}")

Multi-Model Analysis
====================

For robust analyses, consider using multiple models:

.. code-block:: python

   models = [
       pygmm.CampbellBozorgnia2014(),
       pygmm.AbrahamsonSilvaKamai2014(),
       pygmm.BooreStewartSeyhanAtkinson2014(),
       pygmm.ChiouYoungs2014(),
   ]

   # Calculate center, body, and range (CBR) statistics
   results = []
   for model in models:
       ln_sa, ln_std = model(scenario)
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
