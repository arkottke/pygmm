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

    .. grid-item-card:: 🔧 Utilities
        :link: pygmm.tools
        :link-type: ref

        Helper functions and tools for ground motion analysis.

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

        **Utilities**

        - :ref:`pygmm.tools` - Analysis tools
        - :ref:`pygmm.baker_jayaram_2008` - Conditional spectra

    .. tab-item:: Alphabetical

        .. currentmodule:: pygmm

        .. autosummary::

           abrahamson_bhasin_2020
           abrahamson_gregor_addo_2016
           abrahamson_silva_1996
           abrahamson_silva_kamai_2014
           afshari_stewart_2016
           akkar_sandikkaya_bommer_2014
           atkinson_boore_2006
           baker_jayaram_2008
           bayless_abrahamson_2018
           bayless_abrahamson_2019
           boore_stewart_seyhan_atkinson_2014
           campbell_2003
           campbell_bozorgnia_2014
           chiou_youngs_2014
           coppersmith_bommer_2014
           derras_bard_cotton_2014
           gulerce_abrahamson_2011
           hermkes_kuehn_riggelsen_2014
           idriss_2014
           kempton_stewart_2006
           kishida_2017
           model
           pezeshk_zandieh_tavakoli_2011
           stafford_2017
           tavakoli_pezeshk_2005
           tools

Module Documentation
=====================

.. toctree::
   :maxdepth: 2

   pygmm
