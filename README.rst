pyGMM
=====

|PyPi Cheese Shop| |Build Status| |Code Quality| |Test Coverage| |License| |DOI|

Ground motion models implemented in Python.

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

Please cite this software using the DOI_.

.. _DOI: https://zenodo.org/badge/latestdoi/21452/arkottke/pygmm


.. |PyPi Cheese Shop| image:: https://img.shields.io/pypi/v/pygmm.svg
   :target: https://img.shields.io/pypi/v/pygmm.svg
.. |Build Status| image:: https://travis-ci.com/arkottke/pygmm.svg?branch=master
   :target: https://travis-ci.com/arkottke/pygmm
.. |Code Quality| image:: https://api.codacy.com/project/badge/Grade/abc9878c890143c8b590e6f3602056b7
   :target: https://www.codacy.com/manual/arkottke/pygmm
.. |Test Coverage| image:: https://api.codacy.com/project/badge/Coverage/abc9878c890143c8b590e6f3602056b7
   :target: https://www.codacy.com/manual/arkottke/pygmm
.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
.. |DOI| image:: https://zenodo.org/badge/183696586.svg
   :target: https://zenodo.org/badge/latestdoi/183696586
