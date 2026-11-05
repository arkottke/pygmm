Ground-motion models
====================

A ground-motion model is a class that takes a :class:`~pygmm.model.Scenario` and
exposes its prediction as attributes.

.. code-block:: python

   model = pygmm.ChiouYoungs2014(scenario)
   model.periods       # periods of the coefficients (s)
   model.spec_accels   # median 5%-damped spectral acceleration (g)
   model.ln_stds       # standard deviation of the natural log of spec_accels
   model.pga, model.pgv

``interp_spec_accels(periods)`` and ``interp_ln_stds(periods)`` interpolate the
spectrum, in log period and log amplitude, to other periods. The standard
deviation is of the natural logarithm; a spectrum at a number of standard
deviations from the median is ``spec_accels * exp(n * ln_stds)``.

Arrays of scenarios
-------------------

The ground-motion models accept arrays of scenario values, so many scenarios
are computed at once. For N scenarios, ``spec_accels`` and ``ln_stds`` have the
shape (N, periods) and ``pga`` has the shape (N,); scalar values are broadcast
against the arrays. The ``ims`` argument limits the computation to the requested
intensity measures, for example ``"pga"``, ``"psa_1p000"`` (1.0 s) or
``"psa_ngawest2_21"`` (the 21 NGA-West2 comparison periods).

.. code-block:: python

   scenario = pygmm.Scenario(
       mag=np.array([5.5, 6.5, 7.5]), dist_jb=20.0, v_s30=400.0, mechanism="SS"
   )
   model = pygmm.BooreStewartSeyhanAtkinson2014(scenario, ims=["pga", "psa_1p000"])
   model.spec_accels   # shape (3, 1)

USGS nshmp-lib variants
-----------------------

The NGA-West2, NGA-East and NGA-Subduction models, and Atkinson and Macias
(2009) and Zhao et al. (2006), follow the USGS nshmp-lib implementations used by
the National Seismic Hazard Model. Constructor options (for example ``basin``,
``m9`` or ``adjusted``) select the nshmp-lib variants, and the ``GMM_IDS``
attribute of each model maps the nshmp-lib identifiers to the options that
reproduce them.

Finding models
--------------

Every model is registered with what it provides, the tectonic setting it was
developed for, and the kind of input it takes. :func:`pygmm.find_models` filters
the registry and :func:`pygmm.get_model` returns a class from its name.

.. code-block:: python

   psa = pygmm.find_models(provides="psa", input="scenario")
   [info.key for info in psa]

   cls = pygmm.get_model("ChiouYoungs2014")

``provides`` is one of ``psa``, ``pga``, ``pgv``, ``pgd``, ``fas``, ``duration``,
``arias``, ``cav``, ``vh_ratio``, ``correlation`` and ``soil_curves``. ``input`` is
``scenario`` (a :class:`~pygmm.model.Scenario`), ``kwargs`` (keywords, as for the
source-theory model), ``conditional`` (the model is conditioned on another
quantity, such as PGA) or ``stateless`` (functions of the periods only).

Registered models
-----------------

.. list-table::
   :header-rows: 1
   :widths: 40 20 20 20

   * - Model
     - Provides
     - Tectonic setting
     - Input
   * - :class:`~pygmm.AbrahamsonSilvaKamai2014`
     - psa, pga, pgv
     - active crustal
     - scenario
   * - :class:`~pygmm.BooreStewartSeyhanAtkinson2014`
     - psa, pga, pgv
     - active crustal
     - scenario
   * - :class:`~pygmm.CampbellBozorgnia2014`
     - psa, pga, pgv
     - active crustal
     - scenario
   * - :class:`~pygmm.ChiouYoungs2014`
     - psa, pga, pgv
     - active crustal
     - scenario
   * - :class:`~pygmm.Idriss2014`
     - psa, pga
     - active crustal
     - scenario
   * - :class:`~pygmm.AkkarSandikkayaBommer2014`
     - psa, pga, pgv
     -
     - scenario
   * - :class:`~pygmm.DerrasBardCotton2014`
     - psa, pga, pgv
     -
     - scenario
   * - :class:`~pygmm.AtkinsonBoore2006`
     - psa, pga, pgv, pgd
     - stable continental
     - scenario
   * - :class:`~pygmm.Campbell2003`
     - psa
     - stable continental
     - scenario
   * - :class:`~pygmm.PezeshkZandiehTavakoli2011`
     - psa, pga
     - stable continental
     - scenario
   * - :class:`~pygmm.TavakoliPezeshk05`
     - psa, pga
     - stable continental
     - scenario
   * - :class:`~pygmm.NgaEast`
     - psa, pga, pgv
     - stable continental
     - scenario
   * - :class:`~pygmm.NgaEastSeeds`
     - psa, pga
     - stable continental
     - scenario
   * - :class:`~pygmm.AbrahamsonGregorAddo2016`
     - psa, pga
     - subduction
     - scenario
   * - :class:`~pygmm.CoppersmithBommer2014`
     - psa, pga
     - subduction
     - scenario
   * - :class:`~pygmm.AbrahamsonGulerce2020`
     - psa, pga, pgv
     - subduction
     - scenario
   * - :class:`~pygmm.KuehnEtAl2020`
     - psa, pga, pgv
     - subduction
     - scenario
   * - :class:`~pygmm.ParkerEtAl2020`
     - psa, pga, pgv
     - subduction
     - scenario
   * - :class:`~pygmm.AtkinsonMacias2009`
     - psa, pga, pgv
     - subduction
     - scenario
   * - :class:`~pygmm.ZhaoEtAl2006`
     - psa, pga, pgv
     - subduction
     - scenario

The registry also holds the duration, Fourier, velocity-to-horizontal, Arias
intensity, CAV and correlation models described below; ``pygmm.find_models()``
is the complete list.

Other quantities
----------------

**Duration**
    :class:`~pygmm.AfshariStewart2016`, :class:`~pygmm.KemptonStewart2006`,
    :class:`~pygmm.AbrahamsonSilva1996`, :class:`~pygmm.PinillaRamosEtAl2023` and
    :class:`~pygmm.PinillaRamosEtAl2024` take a scenario and expose the duration
    (``duration``, or ``interp`` for a normalized Arias intensity).

**Conditional models**
    :class:`~pygmm.AbrahamsonBhasin2020` (PGV), :class:`~pygmm.AbrahamsonShiYang2016`
    (Arias intensity) and :class:`~pygmm.MacedoAbrahamsonLiu2021` (CAV) are
    conditioned on other ground-motion quantities rather than on a scenario.

**Vertical-to-horizontal ratio**
    :class:`~pygmm.GulerceAbrahamson2011`.

**Spectral correlation**
    :class:`~pygmm.BaylessAbrahamson2018` and :class:`~pygmm.Stafford2017` are
    registered as ``correlation`` models that take no scenario.

Examples
--------

- :doc:`../examples/example-01`: the NGA-West2 models for one scenario and their
  scaling with magnitude, distance and :math:`V_{s30}`.
- :doc:`../examples/example-02`: the spread between the models and the standard
  deviation within them.
- :doc:`../examples/example-06`: a duration model.
