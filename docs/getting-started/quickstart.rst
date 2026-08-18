Quickstart
==========

This guide walks through the basic workflow for evaluating ground motion predictions with pyGMM.

Defining a Scenario
-------------------

All ground motion models take a :class:`~pygmm.model.Scenario` that specifies the earthquake
and site parameters.  Only parameters relevant to a given model need to be provided.

.. code-block:: python

    import pygmm

    scenario = pygmm.Scenario(
        mag=6.5,          # moment magnitude
        dist_rup=20.0,    # rupture distance (km)
        v_s30=760.0,      # time-averaged Vs in top 30 m (m/s)
        mechanism="SS",   # fault mechanism: 'SS', 'NS', 'RS', or 'U'
    )

Running a Model
---------------

Construct a model with the scenario, then read the results off its attributes.
``spec_accels`` is the **median** spectral acceleration in g and ``ln_stds`` the
**log-standard deviation**, both over the model's spectral periods:

.. code-block:: python

    model = pygmm.CampbellBozorgnia2014(scenario)
    ln_sa, ln_std = np.log(model.spec_accels), model.ln_stds

    import numpy as np
    sa = np.exp(ln_sa)   # median spectral acceleration (g)

Access the periods with ``model.PERIODS``.

Comparing Models
----------------

Because all models share the same interface, comparing predictions is straightforward:

.. code-block:: python

    import matplotlib.pyplot as plt

    models = [
        pygmm.CampbellBozorgnia2014(scenario),
        pygmm.BooreStewartSeyhanAtkinson2014(scenario),
        pygmm.AbrahamsonSilvaKamai2014(scenario),
    ]

    fig, ax = plt.subplots()
    for m in models:
        ln_sa, _ = m(scenario)
        ax.semilogx(m.PERIODS, np.exp(ln_sa), label=m.NAME)

    ax.set_xlabel("Period (s)")
    ax.set_ylabel("Spectral Acceleration (g)")
    ax.legend()
    plt.show()

Duration Models
---------------

Duration models (e.g. significant duration D5-75) follow the same pattern:

.. code-block:: python

    dur_model = pygmm.AfshariStewart2016(scenario)
    ln_dur, ln_std = dur_model(scenario)
    d575_median = np.exp(ln_dur)   # seconds

Soil Curves
-----------

Strain-dependent modulus reduction and damping curves can be used directly with pyStrata:

.. code-block:: python

    curves = pygmm.DarendeliSoilType(
        p_atm=1.0,
        sigma_v0=50.0,   # vertical effective stress (kPa)
        PI=0,            # plasticity index
        OCR=1,           # over-consolidation ratio
        freq=1.0,
        n_cycles=10,
    )
    print(curves.strains, curves.mod_reduc, curves.damping)

Next Steps
----------

- :doc:`../user-guide/models` — complete list of available models
- :doc:`../user-guide/usage` — detailed usage patterns
- :doc:`../examples/index` — worked examples
