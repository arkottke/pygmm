.. _usage:

User Guide
==========

pyGMM evaluates empirical models. Each model takes the description of an
earthquake and a site and returns a predicted quantity (a response spectrum, a
Fourier amplitude spectrum, a duration, a soil curve, a velocity profile or a
slope displacement). This page describes the pieces and how they are combined;
the pages that follow describe each piece.

.. toctree::
   :maxdepth: 1

   scenario
   ground-motion
   fourier
   site
   slope-displacement

The pieces
----------

**Scenario** (:doc:`scenario`)
    :class:`~pygmm.model.Scenario` is a dictionary of the earthquake and site
    parameters, such as ``mag``, ``dist_rup`` and ``v_s30``. It is the input of
    the ground-motion, duration and Fourier models.

**Ground-motion models** (:doc:`ground-motion`)
    A class such as :class:`~pygmm.ChiouYoungs2014` is constructed with a
    scenario. The prediction is read from its attributes: ``periods`` and
    ``spec_accels`` (median, in g) and ``ln_stds`` (the log standard deviation)
    for a response-spectrum model, or ``duration`` for a duration model.
    :func:`pygmm.find_models` lists the models by what they provide.

**Fourier amplitude spectrum models** (:doc:`fourier`)
    Models that take a scenario (:class:`~pygmm.BaylessAbrahamson2019`) or
    seismological parameters (:class:`~pygmm.fourier_spectrum.SourceTheoryModel`)
    and expose ``freqs`` and ``fourier_amps``. The site attenuation :math:`\kappa`
    comes from the functions in :mod:`pygmm.kappa`.

**Site models** (:doc:`site`)
    Soil-type classes (:class:`~pygmm.DarendeliSoilType` and others) produce
    modulus reduction and damping curves with ``.curves()``. Profile functions
    (:func:`~pygmm.velocity_profile.kea16_profile` and others) produce a shear-wave
    velocity profile. Both return the dataclasses of :mod:`pygmm.contracts`.

**Slope displacement models** (:doc:`slope-displacement`)
    Functions of the yield coefficient and of ground-motion parameters, which
    are the quantities that the ground-motion models predict, and the transfer
    function of Hale and Abrahamson (2019).

How they fit together
---------------------

.. code-block:: text

   Scenario ──► ground-motion model ──► pga, pgv, spec_accels ──► slope displacement
       │                                                      └──► (pyRotD, hazard, ...)
       ├──────► duration model ───────► duration ───────────────► RvtMotion (pyRVT)
       ├──────► Fourier model ────────► freqs, fourier_amps ────► RvtMotion.from_fas
   site parameters ──► soil-type class ──► NonlinearSoilCurves ───► SoilType.from_curves (pyStrata)
   Vs30, depths ─────► profile function ─► VelocityProfile ───────► Profile.from_velocity_profile

A composed example: a ground-motion model predicts the spectral acceleration at
the degraded period of a slope, which is the input of a slope displacement model.

.. code-block:: python

   import pygmm
   from pygmm import slope_displacement

   scenario = pygmm.Scenario(
       mag=6.5, dist_rup=20.0, dist_jb=20.0, dist_x=20.0, dip=90.0,
       mechanism="SS", v_s30=760.0,
   )
   model = pygmm.ChiouYoungs2014(scenario)

   period_slide = 0.3  # s
   psa_dts = model.interp_spec_accels([1.5 * period_slide])[0]
   disp = slope_displacement.calc_disp_bt07(
       yield_coef=0.1, period_slide=period_slide, psa_dts=psa_dts, mag=6.5
   )

Other libraries
---------------

pyGMM does not depend on the other libraries. They use its output through the
dataclasses in :mod:`pygmm.contracts` (``NonlinearSoilCurves`` and
``VelocityProfile``), which pyStrata accepts by duck typing, and through any
object with ``freqs`` and ``fourier_amps`` for pyRVT.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - pyGMM
     - Library
     - Use
   * - Soil-type classes
     - pyStrata
     - ``SoilType.from_curves(model.curves())``
   * - Profile functions
     - pyStrata
     - ``Profile.from_velocity_profile(profile, soil_types)``
   * - Fourier models
     - pyRVT, pyStrata
     - ``RvtMotion.from_fas(model)``
   * - Response-spectrum models
     - pyStrata
     - target of ``CompatibleRvtMotion``

Examples
--------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Topic
     - Example
   * - Scenario spectra of the NGA-West2 models and their scaling
     - :doc:`../examples/example-01`
   * - Spread between models and within a model
     - :doc:`../examples/example-02`
   * - Fourier amplitude spectra and :math:`\kappa`
     - :doc:`../examples/example-03`
   * - Soil curves and velocity profiles
     - :doc:`../examples/example-04`
   * - Slope displacement models
     - :doc:`../examples/example-05`
   * - Duration (Abrahamson and Silva, 1996)
     - :doc:`../examples/example-06`
