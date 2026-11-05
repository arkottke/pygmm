# pyGMM

[![PyPI](https://img.shields.io/pypi/v/pygmm.svg)](https://pypi.org/project/pygmm/)
[![Python](https://img.shields.io/python/required-version-toml?tomlFilePath=https://raw.githubusercontent.com/arkottke/pygmm/main/pyproject.toml)](https://pypi.org/project/pygmm/)
[![Documentation](https://readthedocs.org/projects/pygmm/badge/?version=latest)](https://pygmm.readthedocs.io/?badge=latest)
[![Build Status](https://github.com/arkottke/pygmm/actions/workflows/python-app.yml/badge.svg)](https://github.com/arkottke/pygmm/actions/workflows/python-app.yml)
[![Code Quality](https://api.codacy.com/project/badge/Grade/abc9878c890143c8b590e6f3602056b7)](https://app.codacy.com/gh/arkottke/pygmm/dashboard)
[![Test Coverage](https://api.codacy.com/project/badge/Coverage/abc9878c890143c8b590e6f3602056b7)](https://app.codacy.com/gh/arkottke/pygmm/dashboard)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/arkottke/pygmm/blob/main/LICENSE)
[![DOI](https://zenodo.org/badge/53176693.svg)](https://zenodo.org/badge/latestdoi/53176693)

Ground motion models implemented in Python.

pyGMM evaluates published empirical models through a common interface, so they
can be compared with each other and used as input to site response and random
vibration theory calculations.

## Features

- Ground motion models for response spectra, PGA and PGV, including the
  NGA-West2, NGA-East and subduction models, and the variants used in the USGS
  National Seismic Hazard Model.
- Duration, Arias intensity, CAV, vertical-to-horizontal ratio and spectral
  correlation models.
- Fourier amplitude spectrum models, including source theory and Stafford et
  al. (2022).
- Nonlinear soil curves and velocity profiles.
- Slope displacement models.
- Vectorized: a scenario can hold arrays, and `ims` selects the intensity
  measures to compute.

See the [list of models](https://pygmm.readthedocs.io/en/latest/user-guide/ground-motion.html)
for what is available. Most models are tested with unit tests that check the
implementation of the model.

## Installation

pyGMM requires Python 3.13 or later and can be installed with
[pip](https://pip.pypa.io):

```bash
python -m pip install pygmm
```

The slope displacement models need the optional extra:

```bash
python -m pip install "pygmm[slope]"
```

## Example

```python
import pygmm

scenario = pygmm.Scenario(
    mag=6.5, dist_rup=20.0, dist_jb=20.0, dist_x=20.0, dip=90.0,
    mechanism="SS", v_s30=760.0,
)
model = pygmm.ChiouYoungs2014(scenario)
model.periods, model.spec_accels, model.ln_stds
```

## Documentation

Information on the usage can be found in the
[documentation](https://pygmm.readthedocs.io/). The
[changelog](https://github.com/arkottke/pygmm/blob/main/CHANGELOG.md) lists the
changes in each release.

The ground motion models of GEM's
[OpenQuake Hazardlib](https://github.com/gem/oq-hazardlib) cover additional
models and are also worth checking out.

## Related packages

pyGMM is developed together with:

- [pyRVT](https://github.com/arkottke/pyrvt): random vibration theory
- [pyRotd](https://github.com/arkottke/pyrotd): response spectra and rotated
  response spectra from acceleration time series
- [pyStrata](https://github.com/arkottke/pystrata): site response analysis

## Citation

When citing the software, please reference the
[DOI](https://zenodo.org/badge/latestdoi/53176693).

## Contributors

- Albert Kottke
- Greg Lavrentiadis
- Artie Rodgers

## License

pyGMM is made available under the MIT License. See
[LICENSE](https://github.com/arkottke/pygmm/blob/main/LICENSE).
