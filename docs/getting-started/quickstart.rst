Quickstart
==========

A scenario, a model and the model's attributes: the minimum to predict a response
spectrum. :doc:`../user-guide/index` describes how the other parts are added.

Define a scenario
-----------------

.. code-block:: python

    import pygmm

    scenario = pygmm.Scenario(
        mag=6.5,          # moment magnitude
        dist_rup=20.0,    # rupture distance (km)
        dist_jb=20.0,     # Joyner-Boore distance (km)
        dist_x=20.0,      # distance perpendicular to the strike (km)
        dip=90.0,         # fault dip (degrees)
        mechanism="SS",   # "SS", "NS", "RS" or "U"
        v_s30=760.0,      # time-averaged shear-wave velocity of the top 30 m (m/s)
    )

A model needs a subset of the parameters and raises ``ValueError`` for a required
parameter that is missing; see :doc:`../user-guide/scenario`.

Run a model
-----------

.. code-block:: python

    model = pygmm.ChiouYoungs2014(scenario)
    model.periods       # periods (s)
    model.spec_accels   # median 5%-damped spectral acceleration (g)
    model.ln_stds       # standard deviation of ln(spec_accels)
    model.pga, model.pgv

Compare models
--------------

All of the response-spectrum models have the same attributes.

.. code-block:: python

    import matplotlib.pyplot as plt

    classes = [
        pygmm.AbrahamsonSilvaKamai2014,
        pygmm.BooreStewartSeyhanAtkinson2014,
        pygmm.CampbellBozorgnia2014,
        pygmm.ChiouYoungs2014,
    ]

    fig, ax = plt.subplots()
    for cls in classes:
        m = cls(scenario)
        ax.plot(m.periods, m.spec_accels, label=cls.ABBREV)
    ax.set(xscale="log", yscale="log", xlabel="Period (s)", ylabel="Spectral acceleration (g)")
    ax.legend()

Other predictions
-----------------

A duration model is used the same way:

.. code-block:: python

    duration = pygmm.AfshariStewart2016(scenario).duration  # durations (s)

and a soil-type class returns modulus reduction and damping curves:

.. code-block:: python

    curves = pygmm.DarendeliSoilType(unit_wt=18.0, plas_index=20, stress_mean=100.0).curves()
    curves.strains, curves.mod_reduc, curves.damping

Next steps
----------

- :doc:`../user-guide/index`: how the parts fit together.
- :doc:`../examples/index`: worked examples.
- :doc:`../api/index`: the API reference.
