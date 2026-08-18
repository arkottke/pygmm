=============
API Reference
=============

This section provides detailed documentation for all modules, classes, and functions in pyGMM.

.. grid:: 1 2 2 2
    :gutter: 3

    .. grid-item-card:: 🎯 Core Module
        :link: pygmm
        :link-type: doc

        Main module containing all ground motion models and utilities.

    .. grid-item-card:: 🏗️ Model Framework
        :link: pygmm.model
        :link-type: ref

        Base classes and interfaces for ground motion models.

    .. grid-item-card:: 📊 Scenarios
        :link: pygmm.model.Scenario
        :link-type: ref

        Earthquake scenario definition and management.

    .. grid-item-card:: 🔎 Model Registry
        :link: pygmm.registry
        :link-type: ref

        Discover models by capability with ``find_models()``.

Quick Navigation
================

.. tab-set::

    .. tab-item:: By Category

        **Ground Motion Models**

        - :doc:`models` - Overview and selection guide
        - :ref:`pygmm.model.GroundMotionModel` - Base class

        **Scenario Definition**

        - :ref:`pygmm.model.Scenario` - Earthquake scenarios
        - :ref:`pygmm.model.Parameter` - Parameter validation

        **Discovery**

        - :func:`pygmm.find_models` - Filter models by capability
        - :func:`pygmm.get_model` - Look up a model by name

        **Utilities**

        - :mod:`pygmm.registry` - Capability-based model discovery
        - :mod:`pygmm.ground_motion.baker_jayaram_2008` - Conditional spectra

    .. tab-item:: Alphabetical

        .. currentmodule:: pygmm

        .. autosummary::

           contracts
           ground_motion
           fourier_spectrum
           model
           registry
           soil_curves

    .. tab-item:: By package

        .. currentmodule:: pygmm

        - :mod:`pygmm.ground_motion` -- response spectra, peak parameters,
          duration, Arias intensity, CAV, V/H ratios and inter-period
          correlation. Independent variable is period, or the output is a
          scalar conditioned on a :class:`~pygmm.model.Scenario`.
        - :mod:`pygmm.fourier_spectrum` -- Fourier amplitude spectra and
          inter-frequency correlation. Independent variable is frequency.
        - :mod:`pygmm.soil_curves` -- modulus-reduction and damping curves.
          Independent variable is shear strain.
        - :mod:`pygmm.registry` -- capability-based model discovery
          (:func:`~pygmm.find_models`, :func:`~pygmm.get_model`).
        - :mod:`pygmm.contracts` -- output dataclasses and the producer
          protocols consumers duck-type against.

Module Documentation
=====================

.. toctree::
   :maxdepth: 2

   pygmm
