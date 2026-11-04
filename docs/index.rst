=====================================
pyGMM: Ground Motion Models in Python
=====================================

.. grid:: 1 2 2 2
    :gutter: 2

    .. grid-item-card:: 🚀 Quick Start
        :link: getting-started/quickstart
        :link-type: doc

        Get up and running with pyGMM in minutes. Install the package and run your first ground motion calculation.

    .. grid-item-card:: 📚 User Guide
        :link: user-guide/index
        :link-type: doc

        How the scenario, the models and their outputs fit together.

    .. grid-item-card:: 🔬 Ground Motion Models
        :link: user-guide/ground-motion
        :link-type: doc

        The ground-motion models, and how to find them by what they predict.

    .. grid-item-card:: 📖 API Reference
        :link: api/index
        :link-type: doc

        Detailed documentation of all classes, functions, and modules in pyGMM.

.. admonition:: What is pyGMM?
   :class: important

   pyGMM is a Python library for evaluating published empirical models: ground-motion
   models, Fourier amplitude spectrum models, soil curves, velocity profiles and
   slope-displacement models. The models share an interface, so they can be compared
   and used as input to site response and random vibration theory calculations.

Features
========

.. grid:: 1 2 3 3
    :gutter: 2

    .. grid-item::
        :columns: 12 6 4 4

        **🎯 Unified Interface**

        Consistent API across all ground motion models for easy comparison and analysis.

    .. grid-item::
        :columns: 12 6 4 4

        **📊 Multiple Models**

        Support for dozens of published ground motion prediction equations.

    .. grid-item::
        :columns: 12 12 4 4

        **🔧 Easy Integration**

        Simple installation and integration with existing Python workflows.

Getting Started
===============

Install pyGMM using pip:

.. code-block:: bash

   pip install pygmm

Quick example:

.. code-block:: python

   import pygmm

   scenario = pygmm.Scenario(
       mag=6.5, dist_rup=20.0, dist_jb=20.0, dist_x=20.0, dip=90.0,
       mechanism="SS", v_s30=760.0,
   )
   model = pygmm.ChiouYoungs2014(scenario)
   model.periods, model.spec_accels, model.ln_stds

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Getting Started

   getting-started/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   user-guide/index
   examples/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference

   api/index
   references

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Development

   developer-guide/index

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: About

   about/authors
   about/readme
   about/license
