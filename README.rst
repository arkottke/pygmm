=====
pyGMM
=====

.. image:: https://img.shields.io/pypi/v/pygmm.svg
    :target: https://pypi.python.org/pypi/pygmm
    :alt: PyPi Cheese Shop

.. image:: https://img.shields.io/travis/arkottke/pygmm.svg
    :target: https://travis-ci.org/arkottke/pygmm
    :alt: Build Status

.. image:: https://readthedocs.org/projects/pygmm/badge/?version=latest&style=flat
    :target: https://pygmm.readthedocs.org
    :alt: Documentation Status

.. image:: https://coveralls.io/repos/github/arkottke/pygmm/badge.svg?branch=master
    :target: https://coveralls.io/github/arkottke/pygmm?branch=master
    :alt: Test Coverage

.. image:: https://landscape.io/github/arkottke/pygmm/master/landscape.svg?style=flat
    :target: https://landscape.io/github/arkottke/pygmm/master
    :alt: Code Health

.. image:: https://img.shields.io/badge/license-MIT-blue.svg
    :target: https://github.com/arkottke/pygmm/blob/master/LICENSE
    :alt: License
    
.. image:: https://zenodo.org/badge/21452/arkottke/pygmm.svg
   :target: https://zenodo.org/badge/latestdoi/21452/arkottke/pygmm
   :alt: DOI Information

Ground motion models implemented in Python.

* Free software: MIT license
* Documentation: https://pygmm.readthedocs.org.

I have recently learned that additional ground motion models have been implemented through GEM's OpenQuake Hazardlib_, which I recommend checking out.

.. _Hazardlib: https://github.com/gem/oq-hazardlib

Features
--------

Models currently supported:

* Akkar, Sandikkaya, & Bommer (2014) with unit tests

* Atkinson & Boore (2006)

* Abrahamson, Silva, & Kamai (2014) with unit tests

* Abrahamson, Gregor, & Addo (2016) with unit tests

* Boore, Stewart, Seyhan, & Atkinson (2014) with unit tests

* Campbell (2003)

* Campbell & Bozorgnia (2014) with unit tests

* Chiou & Youngs (2014) with unit tests

* Derras, Bard & Cotton (2013) with unit tests

* Idriss (2014) with unit tests

* Pezeshk, Zandieh, & Tavakoli (2001)

* Tavakoli & Pezeshk (2005)

Conditional spectra models:

* Baker & Jayaram (2008) with unit tests

* Kishida (2017) with unit tests

Unit tests means that each test cases are used to test the implemention of
the model.

Citation
--------

Please cite this software using the following DOI:

.. image:: https://zenodo.org/badge/21452/arkottke/pygmm.svg
   :target: https://zenodo.org/badge/latestdoi/21452/arkottke/pygmm
   :alt: DOI Information
