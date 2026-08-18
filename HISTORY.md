---
title: History
---

# 0.9.0 (unreleased)
- Changed: adopted [SPEC 0](https://scientific-python.org/specs/spec-0000/) for
  version support -- Python minors from the last 36 months, core dependencies
  from the last 24 months. **Minimum Python is now 3.12** (was 3.10), and
  dependency lower bounds are declared explicitly: `matplotlib>=3.10`,
  `numpy>=2.1`, `pandas>=2.3`, `pint>=0.25`, `scipy>=1.15`. CI now runs 3.12,
  3.13 and 3.14. The dated drop schedule is in `pyproject.toml`.

  The old `>=3.10` floor was already unmet in practice: numpy, scipy, pandas,
  matplotlib and pint all require >=3.11, so a 3.10 install silently resolved
  to a much older stack (pandas 2.3 rather than 3.0, numpy 2.2 rather than
  2.5). That forced the lockfile to fork on `python_full_version < '3.11'` into 5 duplicated packages.

- Added: Macedo, Abrahamson, & Liu (2021) conditional and scenario-based CAV models for shallow crustal settings
- Added: Abrahamson, Shi, & Yang (2016) Ground-Motion Prediction Equations for {A}rias Intensity Consistent with the NGA-West2 Ground-Motion Models
- Added: `pygmm.contracts` dataclasses (`FourierSpectrum`, `ResponseSpectrum`,
  `NonlinearSoilCurves`, `VelocityProfile`) so downstream packages can interoperate
  without depending on pygmm at runtime
- Added: `pygmm.fourier_spectrum` with the source-theory and Stafford et al. (2022)
  FAS models, moved out of pyRVT
- Added: `pygmm.soil_curves` with the empirical nonlinear soil-curve models
  (Darendeli, Menq, Wang, Alemu et al., Rollins et al., Kishida), moved out of pyStrata
- Added: `kea16_profile()`, moved out of pyStrata
- Added: `pygmm.registry` — `find_models()` and `get_model()` for capability-based
  discovery, e.g. `find_models(provides="psa", tectonic="active_crustal")`
- Added: producer protocols in `pygmm.contracts` (`SupportsResponseSpectrum`,
  `SupportsFourierSpectrum`, `SupportsDuration`, `SupportsSoilCurves`), replacing
  the seven `_base.py` abstract base classes — six of which never acquired a
  subclass. Conformance is structural, matching how pyStrata and pyRVT already
  consume pygmm.
- Added: `duration_model()` on the five duration models and `fourier_spectrum()`
  on the three FAS models, so they emit the corresponding contract
- Added: `fourier_amps` on `BaylessAbrahamson2019`, which previously exposed only
  `eas` and so could not be used with pyRVT despite living in `fourier_spectrum`
- Added: `Duration.ln_std`, `Duration.plus_sigma` and `Duration.minus_sigma`. Two
  conventions are carried because the models disagree: AS96/AS16/KS06 report a
  lognormal standard error, while Pinilla-Ramos et al. apply sigma in a
  transformed power space that is not symmetric in log space.
- Changed: models are organized into three packages, chosen by the independent
  variable of what they predict — shear strain (`soil_curves`), frequency
  (`fourier_spectrum`), or period / a scalar conditioned on a `Scenario`
  (`ground_motion`). The top-level names are unchanged. This replaces the nine
  category subpackages introduced earlier in this release, which were never
  published.
- Changed: `pint` is now a runtime dependency, and `soil_curves` is imported
  lazily. `import pygmm` previously failed outright on any install without pint.
- Fixed: `AbrahamsonBhasin2020.__init__` ended with `return ln_mean, ln_std`, so
  the model could not be constructed at all. Results are now exposed as
  `pgv`/`ln_pgv`/`ln_std`/`tau`/`phi`, matching the other conditional models.
- Fixed: `ChiouYoungs2014` no longer silently loses the `tau`/`phi` return added
  in 0.6.6; a duplicate copy of the module had reverted it.
- Fixed: `Stafford2017` referenced coefficient keys absent from its data file in
  one of two duplicate copies. Nine modules existed in duplicate; the stale
  copies are removed.
- Fixed: parameter-limit warnings printed `{self.min}`/`{self.max}` literally.
- Deprecated: `pygmm.models`. It mixed response-spectrum, duration, FAS and
  correlation models with no way to tell them apart, and 7 of its 20 entries had
  no `spec_accels`. Use `find_models()`.
- Removed: `pygmm.contracts.FaultDisplacement`, `CptSounding`,
  `SoilBehaviorProfile` and `LiquefactionTriggering`, which had no producer and
  no consumer
- Removed: `pygmm.tools` (an empty stub) and `pygmm.types` (renamed `_types`; it
  shadowed the stdlib `types` module, breaking any interpreter started with the
  package directory as the working directory)
- Changed: documentation reorganized to the Scientific Python layout

# 0.8.0 (2025-07-24)
- Added: Pinilla-Ramos et al. (2023) model for duration of crustal earthquakes
- Added: Pinilla-Ramos et al. (2024) model for duration of subduction earthquakes
- Added: Stafford (2017) model for FAS correlation

# 0.7.3 (2025-03-12)
-   Fixed Bayless and Abrahamson (2018) correlation model.

# 0.7.1 (2025-03-05)

-   Added compatibility with numpy 2.0

# 0.7.0 (2024-04-24)

-   Added: Abrahamson and Bhasin (2020)
-   Changed to Hatch build system

# 0.6.6 (2023-12-11)

-   Added: Return tau and phi in the standard deviation calculations

# 0.6.5 (2022-09-16)

-   Added: Afshari and Stewart (2016) duration model
-   Added: Kempton and Stewart (2006) duration model

# 0.6.4 (2022-01-24)

-   Added: Bayless and Abrahamson (2019)

# 0.6.3 (2021-12-08)

-   Fixed: error in ASK14 on a7 term

# 0.6.2 (2021-10-19)

-   Changed: Move site amplification to static functions on some GMPEs

# 0.6.1 (2020-06-03)

-   Added Coppersmith & Bommer (2014) model for Hanford
-   Factored tests

# 0.6.0 (2019-08-12)

-   Added Abrahamson, Gregor, Addo (2014)
-   Added Abrahamson & Gulerce (2011)
-   Added conditional mean spectra models.
-   Added Scenario objects.
-   Added typing for all classes.

# 0.4.0 (2016-04-08)

-   Added Hermkes et al. (2014).
-   Improved documentation.
-   Added Baker & Jayaram (2008), Kishida (2017)

# 0.3.2 (2016-03-30) {#section-1}

-   Nothing changed yet.

# 0.3.1 (2016-03-30) {#section-2}

-   First release on PyPI.
