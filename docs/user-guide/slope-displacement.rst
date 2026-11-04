Slope displacement
==================

The slope-displacement models estimate the permanent displacement of a sliding
mass. They are functions of the yield coefficient :math:`k_y`, the period of the
mass, and ground-motion parameters, which are what the ground-motion models of
:doc:`ground-motion` predict.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Function
     - Model
   * - :func:`~pygmm.slope_displacement.calc_disp_rs08`
     - Rathje and Saygili (2008), from PGA (and magnitude or PGV)
   * - :func:`~pygmm.slope_displacement.calc_disp_ra11`
     - Rathje and Antonakos (2011), from PGA and the periods of the mass and of the ground motion
   * - :func:`~pygmm.slope_displacement.calc_disp_bt07`
     - Bray and Travasarou (2007), from the spectral acceleration at the degraded period
   * - :func:`~pygmm.slope_displacement.calc_disp_bea17`
     - Bray et al. (2017), subduction events
   * - :func:`~pygmm.slope_displacement.calc_disp_wla06`
     - Watson-Lamprey and Abrahamson (2006), from PGA and :math:`S_a(1\,\mathrm{s})`
   * - ``calc_disp_cr20``, ``calc_disp_cr21``, ``calc_disp_cr22``
     - Cho and Rathje (2020, 2021, 2022)
   * - ``calc_prob_disp_bt07``, ``calc_prob_disp_bea17``
     - Probability of a non-zero displacement
   * - :func:`~pygmm.slope_displacement.calc_rigid_disp`
     - Displacement of a rigid block from an acceleration time series

The models differ in the ground-motion parameters they take, so the same
ground-motion model supplies each of them with a different quantity.

.. code-block:: python

   from pygmm import slope_displacement

   model = pygmm.ChiouYoungs2014(scenario)
   psa_dts = model.interp_spec_accels([1.5 * 0.3])[0]
   slope_displacement.calc_disp_bt07(
       yield_coef=0.1, period_slide=0.3, psa_dts=psa_dts, mag=6.5
   )

Earth dams
----------

:class:`~pygmm.slope_displacement.HaleAbrahamson19` is the transfer function of
Hale and Abrahamson (2019) from the input motion to the average acceleration of a
sliding mass in an earth dam. It is a function of the height and the shear-wave
velocity of the dam, and takes the 5%-damped spectral acceleration at its nonlinear
period. ``transfer_function(freqs, spec_accel)`` returns a
:class:`~pygmm.slope_displacement.HaleTransferFunction`. The period of a dam can be
estimated with :func:`~pygmm.slope_displacement.calc_dam_period_pk19`. In pyStrata
the same quantity is computed for a shear-wedge dam, and `example 21
<https://pystrata.readthedocs.io/en/latest/examples/example-21.html>`__ compares the two.

Examples
--------

- :doc:`../examples/example-05`: the models for one scenario and a range of yield
  coefficients.
