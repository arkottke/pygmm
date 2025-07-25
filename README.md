# pyGMM

[![Docs](https://readthedocs.org/projects/pygmm/badge/?version=latest)](https://pygmm.readthedocs.io/en/latest/?badge=latest)
[![PyPi Cheese
Shop](https://img.shields.io/pypi/v/pygmm.svg)](https://pypi.org/project/pygmm/)
[![Build
Status](https://github.com/arkottke/pygmm/actions/workflows/python-app.yml/badge.svg)](https://github.com/arkottke/pygmm/actions/workflows/python-app.yml)
[![Code
Quality](https://api.codacy.com/project/badge/Grade/abc9878c890143c8b590e6f3602056b7)](https://app.codacy.com/gh/arkottke/pygmm/dashboard)
[![Test
Coverage](https://api.codacy.com/project/badge/Coverage/abc9878c890143c8b590e6f3602056b7)](https://app.codacy.com/gh/arkottke/pygmm/dashboard)
![License](https://img.shields.io/badge/license-MIT-blue.svg)
[![DOI](https://zenodo.org/badge/53176693.svg)](https://zenodo.org/badge/latestdoi/53176693)

Ground motion models implemented in Python.

I have recently learned that additional ground motion models have been
implemented through GEM\'s OpenQuake
[Hazardlib](https://github.com/gem/oq-hazardlib), which I recommend
checking out.

## Features

Models currently supported:

-   Akkar, Sandikkaya, & Bommer (2014) with unit tests
-   Atkinson & Boore (2006)
-   Abrahamson, Silva, & Kamai (2014) with unit tests
-   Abrahamson, Gregor, & Addo (2016) with unit tests
-   Boore, Stewart, Seyhan, & Atkinson (2014) with unit tests
-   Campbell (2003)
-   Campbell & Bozorgnia (2014) with unit tests
-   Chiou & Youngs (2014) with unit tests
-   Derras, Bard & Cotton (2013) with unit tests
-   Idriss (2014) with unit tests
-   Pezeshk, Zandieh, & Tavakoli (2001)
-   Tavakoli & Pezeshk (2005)

Conditional spectra models:

-   Baker & Jayaram (2008) with unit tests
-   Kishida (2017) with unit tests

Duration models:

-   Kempton and Stewart (2006)
-   Afshari and Stewart (2016)

Most models are tested with unit tests that test the implemention of the
model.

## Citation

Please cite this software using the
[DOI](https://zenodo.org/badge/latestdoi/53176693).

## Contributors

-   Albert Kottke
-   Greg Lavrentiadis
-   Artie Rodgers
