Examples
========

The examples are notebooks in the ``examples`` directory of the repository, and
are executed by the test suite. They are grouped by the part of the library that
they use; :doc:`../user-guide/index` describes how the parts fit together.

Ground-motion models
--------------------

**NGA-West2 models for a scenario** (:doc:`example-01`)
    The response spectra and standard deviations of the five NGA-West2 models for
    one scenario, and their scaling with magnitude, distance and :math:`V_{s30}`.

**Spread between ground-motion models** (:doc:`example-02`)
    The spread of the medians of the models, as a function of period and distance,
    compared with the standard deviation of each model.

**Duration** (:doc:`example-06`)
    The duration model of Abrahamson and Silva (1996) as a function of the
    normalized Arias intensity.

Fourier amplitude spectra
-------------------------

**Fourier spectra and** :math:`\kappa` (:doc:`example-03`)
    Three Fourier amplitude spectrum models, the site attenuation relationships and
    the response spectrum of a source-theory motion computed with pyRVT.

Site models
-----------

**Soil curves and velocity profiles** (:doc:`example-04`)
    The modulus reduction and damping curves of the soil-type models, their
    dependence on confining stress, and two shear-wave velocity profile models.

Slope displacement
------------------

**Slope displacement models** (:doc:`example-05`)
    Five slope-displacement models, with the ground-motion parameters from a
    ground-motion model, for a range of yield coefficients, magnitudes and distances.

.. toctree::
   :maxdepth: 1
   :hidden:

   example-01
   example-02
   example-03
   example-04
   example-05
   example-06
