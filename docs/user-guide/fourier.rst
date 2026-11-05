Fourier amplitude spectra
=========================

A Fourier amplitude spectrum model predicts the amplitude of the Fourier transform
of the acceleration. The models expose ``freqs`` and ``fourier_amps`` (and
``ln_std`` where the model has a standard deviation).

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Model
     - Input
   * - :class:`~pygmm.BaylessAbrahamson2019`
     - :class:`~pygmm.model.Scenario`
   * - :class:`~pygmm.fourier_spectrum.SourceTheoryModel`
     - magnitude, distance and region, with optional stress drop and depth
   * - :class:`~pygmm.fourier_spectrum.StaffordEtAl2022`
     - magnitude, rupture or Joyner-Boore distance and mechanism (point-source model for :math:`V_{s30}` = 760 m/s)

.. code-block:: python

   from pygmm import fourier_spectrum

   fas = fourier_spectrum.SourceTheoryModel(magnitude=6.5, distance=20.0, region="wna")
   fas.freqs, fas.fourier_amps, fas.duration

Site amplification is included by default and ``disable_site_amp=True`` removes
it, so that the spectrum can be used as input at a different reference condition.

Site attenuation
----------------

:mod:`pygmm.kappa` has the relationships between :math:`V_{s30}` and the site
attenuation parameter :math:`\kappa_0`: :func:`~pygmm.kappa.calc_site_atten_vhea11`
and :func:`~pygmm.kappa.calc_site_atten_xr21`.

Correlation models
------------------

:class:`~pygmm.BaylessAbrahamson2018` and :class:`~pygmm.Stafford2017` are registered
as ``correlation`` models, which take no scenario (``input="stateless"``).

Use in pyRVT and pyStrata
-------------------------

A source-theory model has the Fourier amplitudes and the duration that pyRVT needs,
so it is passed to ``RvtMotion.from_fas`` (pyRVT or pyStrata) to compute the peak
response with random vibration theory.

Examples
--------

- :doc:`../examples/example-03`: Fourier spectra from the models, with the effect
  of :math:`\kappa`.
