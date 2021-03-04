======
Models
======


Generic Interface
-----------------

The interfaces for models have been simplified to use same parameter names and
values where possible. Details of this interface are provided in
:class:`~pygmm.model.GroundMotionModel`.

.. currentmodule:: pygmm.model

.. autoclass:: Model

   .. automethod:: __init__

   .. rubric:: Summary of Methods

   .. autosummary::

      ~Model.__init__
      ~Model.interp_ln_stds
      ~Model.interp_spec_accels

   .. rubric:: Attributes

   .. autosummary::

      ~Model.ABBREV
      ~Model.INDEX_PGA
      ~Model.INDEX_PGD
      ~Model.INDEX_PGV
      ~Model.INDICES_PSA
      ~Model.LIMITS
      ~Model.NAME
      ~Model.PARAMS
      ~Model.PERIODS
      ~Model.PGD_SCALE
      ~Model.PGV_SCALE
      ~Model.ln_std_pga
      ~Model.ln_std_pgd
      ~Model.ln_std_pgv
      ~Model.ln_stds
      ~Model.periods
      ~Model.pga
      ~Model.pgd
      ~Model.pgv
      ~Model.spec_accels

Mechanism
.........

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
