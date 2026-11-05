API Reference
=============

The package is organized by the quantity that is predicted. Each model is a
class that takes a :class:`~pygmm.model.Scenario` (or the parameters of a
site) and exposes its results as attributes; see :doc:`../user-guide/index`
for how they fit together.

.. toctree::
   :maxdepth: 1

   core
   ground_motion
   fourier_spectrum
   kappa
   soil_curves
   velocity_profile
   slope_displacement
