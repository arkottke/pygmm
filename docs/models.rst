======
Models
======


Generic Interface
-----------------

The interfaces for models have been simplified to use same parameter names and
values where possible. Details of this interface are provided in
:class:`~pygmm.model.Model`.

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

.. automodule:: pygmm
.. autosummary::
    :toctree: 
    :nosignatures:

    AbrahamsonSilvaKamai2014
    AkkarSandikkayaBommer2014
    AtkinsonBoore2006
    BooreStewartSeyhanAtkinson2014
    Campbell2003
    CampbellBozorgnia2014
    ChiouYoungs2014
    DerrasBardCotton2014
    Idriss2014
    PezeshkZandiehTavakoli2011
    TavakoliPezeshk05

If you are interested in contributing another model to the collection please see
:doc:`contributing`.
