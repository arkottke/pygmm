======
Models
======


Generic Interface
-----------------

The interfaces for models have been simplified to use same parameter names and
values where possible.

.. autoclass:: pygmm.model.Model
    :members:

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

.. automodule:: pygmm
    :members:
    :special-members:
