=====================================
pyGMM: Ground Motion Models in Python
=====================================

.. grid:: 1 2 2 2
    :gutter: 2

    .. grid-item-card:: 🚀 Quick Start
        :link: installation
        :link-type: doc

        Get up and running with pyGMM in minutes. Install the package and run your first ground motion calculation.

    .. grid-item-card:: 📚 User Guide
        :link: usage
        :link-type: doc

        Learn how to use pyGMM effectively with comprehensive examples and tutorials.

    .. grid-item-card:: 🔬 Ground Motion Models
        :link: models
        :link-type: doc

        Explore the available ground motion prediction equations and their implementations.

    .. grid-item-card:: 📖 API Reference
        :link: modules
        :link-type: doc

        Detailed documentation of all classes, functions, and modules in pyGMM.

.. admonition:: What is pyGMM?
   :class: important

   pyGMM is a Python library for evaluating ground motion prediction equations (GMPEs).
   It provides a unified interface to numerous published ground motion models, making it
   easy to compare predictions and perform probabilistic seismic hazard analysis.

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

   # Create a scenario
   scenario = pygmm.Scenario(mag=6.5, dist_rup=20, v_s30=760)

   # Initialize a ground motion model
   gmpe = pygmm.CampbellBozorgnia2014()

   # Calculate ground motion
   mean, std = gmpe(scenario)

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: User Guide

   installation
   usage
   examples/index

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Reference

   models
   modules
   pygmm

.. toctree::
   :maxdepth: 1
   :caption: Development

   contributing
   changelog
   authors

.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: About

   readme
   zreferences
   license
