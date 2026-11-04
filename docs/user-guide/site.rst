Site models
===========

Two kinds of models describe the site: those that produce the nonlinear soil
behavior of a soil type, and those that produce a shear-wave velocity profile.

Soil curves
-----------

A soil-type class takes the properties of a soil, such as the unit weight, the
plasticity index and the mean effective stress, and ``.curves()`` returns a
:class:`~pygmm.contracts.NonlinearSoilCurves` with ``strains``, ``mod_reduc``
and ``damping``.

.. code-block:: python

   import pygmm

   soil = pygmm.DarendeliSoilType(
       unit_wt=18.0, plas_index=20, ocr=1, stress_mean=100.0
   )
   curves = soil.curves()
   curves.strains, curves.mod_reduc, curves.damping

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Class
     - Soil properties
   * - :class:`~pygmm.DarendeliSoilType`
     - plasticity index, OCR, mean effective stress, frequency, number of cycles
   * - :class:`~pygmm.MenqSoilType`
     - gravelly soils (Menq, 2003); coefficient of uniformity, mean grain size, mean effective stress
   * - :class:`~pygmm.RollinsEtAlSoilType`
     - gravels (Rollins et al., 2020); coefficient of uniformity, mean effective stress
   * - :class:`~pygmm.WangSoilType`
     - Wang and Stokoe (2022); a soil group (``clean_sand_and_gravel``, ``nonplastic_silty_sand`` or ``clayey_soil``)
   * - :class:`~pygmm.KishidaSoilType`
     - highly organic soils (Kishida, 2017); organic content, vertical stress
   * - :class:`~pygmm.AlemuEtAlSoilType`
     - transitional silts (Alemu et al., 2025); plasticity index, OCR, mean effective stress, fines content

All of the classes take ``strains`` to set the strains of the curves and
``damping_min`` (except where the model has no such parameter) to override the
small-strain damping.

Velocity profiles
-----------------

The profile functions take depths and :math:`V_{s30}` and return a
:class:`~pygmm.contracts.VelocityProfile` with ``depth``, ``vs_median`` (the median
shear-wave velocity), ``std_vs_ln`` (the standard deviation of its natural log)
and, where the model defines them, ``region`` and ``site_class``.

.. code-block:: python

   import numpy as np
   from pygmm import velocity_profile

   depth = np.linspace(0, 100, 21)
   profile = velocity_profile.sa18_profile(depth, vs30=400.0)

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Function
     - Profile
   * - :func:`~pygmm.velocity_profile.kea16_profile`
     - Median profile of Kamai et al. (2016), for ``"california"`` or ``"japan"``
   * - :func:`~pygmm.velocity_profile.sa18_profile`
     - Median profile of Shi and Asimaki (2018)
   * - :func:`~pygmm.velocity_profile.btc11_profile`
     - Generic profile of Boore et al. (2011)
   * - ``bj97gr_profile``, ``bj97gr760_profile``, ``bj97gvhr_profile``
     - Generic rock profiles (Boore and Joyner, 1997; Boore, 2016)

:func:`~pygmm.velocity_profile.calc_vel_shear_spt_wds12` estimates the shear-wave
velocity from SPT blow counts, and ``calc_density_bea16`` and ``calc_vel_comp_bea16``
the density and compression-wave velocity from the shear-wave velocity.

Use in pyStrata
---------------

.. code-block:: python

   import pystrata

   soil_type = pystrata.site.SoilType.from_curves(curves, unit_wt=18.0)
   profile = pystrata.site.Profile.from_velocity_profile(
       velocity_profile.sa18_profile(depth, vs30=400.0), soil_types=soil_type
   )

Examples
--------

- :doc:`../examples/example-04`: the soil curves and velocity profiles compared.
