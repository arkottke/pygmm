# Changelog

## 0.9.0 (unreleased)

- Added: `HaleAbrahamson19.transfer_function(freqs, spec_accel)` returns the
  transfer function to the average acceleration of the sliding mass, with the
  damping and first-mode frequency that define it, as a `HaleTransferFunction`.
  It takes the 5% damped spectral acceleration [g] at `period_nl`, so it works
  with any motion. It replaces the private `_calc_trans_func`.
- Changed: the documentation is reorganized: a user guide that describes how the
  scenario, the models and their outputs fit together, the API reference in
  `api/`, and six executed example notebooks (NGA-West2 models, the spread
  between models, Fourier spectra and kappa, soil curves and velocity profiles,
  slope displacement, duration) that replace the example scripts and the example
  pages, whose code did not run. The notebooks are executed by the tests.
- Added: `HaleAbrahamson19.from_metric` creates the model from a height in
  meters and a velocity in m/s; the model itself works in feet.

- Changed: adopted [SPEC 0](https://scientific-python.org/specs/spec-0000/) for
  version support -- Python minors from the last 36 months, core dependencies
  from the last 24 months. **Minimum Python is now 3.13** (was 3.10), and
  dependency lower bounds are declared explicitly: `matplotlib>=3.10`,
  `numpy>=2.2`, `pandas>=2.3`, `pint>=0.25`, `scipy>=1.15`. CI now runs 3.13
  and 3.14. The dated drop schedule is in `pyproject.toml`.

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
- Added: the USGS nshmp-lib versions of the Atkinson and Macias (2009) interface model
  (`AtkinsonMacias2009`; `GMM_IDS`: `AM_09_INTERFACE`, `_BASIN`, `_BASIN_M9`, `_BASIN_SITE_FIX`,
  and `_BASIN_M9_SITE_FIX`) and the Zhao et al. (2006) interface and intraslab model
  (`ZhaoEtAl2006`; `ZHAO_06_INTERFACE`, `_INTERFACE_BASIN`, `_INTERFACE_BASIN_M9`, `_SLAB`, and
  `_SLAB_BASIN`, with `event_type` as a scenario value), which are used by the Cascadia
  subduction logic trees of the current CONUS NSHM. The `basin` (USGS CB14 deep basin term), `m9`
  (Seattle M9 adjustment), and, for AM09, `site_fix` (the reference PGA instead of its natural log
  in the BA08 site term) options select the nshmp-lib variants. As in nshmp-lib, Zhao et al.
  (2006) selects the site class from `v_s30`, uses a depth of 20 km for interface events and
  `depth_tor` for intraslab events, and extrapolates 7.5 and 10 s with BC Hydro (2012) and AM09,
  and both models interpolate the periods without coefficients and compute PGV with the
  Abrahamson and Bhasin (2020) conditional model. Both reproduce the nshmp-lib reference results.
- Added: options of the NGA-West2 models (`AbrahamsonSilvaKamai2014`,
  `BooreStewartSeyhanAtkinson2014`, `CampbellBozorgnia2014`, and `ChiouYoungs2014`) for the USGS
  nshmp-lib variants used by the NSHM (see `pygmm.ground_motion.nshmp`): `nshmp` (the nshmp-lib
  implementation, e.g., no basin term for a `depth_1_0` of NaN, the hanging wall from `dist_x`,
  `width` from `depth_bor`, and `depth_hyp` (CB14) from the rupture geometry), `epistemic` (USGS
  NGA-West2 epistemic branches, collapsed as in nshmp-lib, with the branches in `ln_branches`),
  `basin` (USGS deep basin scaling), `cybershake` (2023 NSHM CyberShake basin terms), and `prvi`
  (2025 PRVI coefficients). `GMM_IDS` maps the nshmp-lib `Gmm` ids (`ASK_14`, `_BASE`, `_BASIN`,
  `_CYBERSHAKE`, `_VS30_MEASURED`, `_PRVI`, etc.) to the options and, for ASK14 and CY14, the
  `vs_source` scenario value. The active crust logic tree of the current CONUS NSHM uses the four
  `_BASIN` variants with weights of 0.25. The options are off by default, which gives the published
  models (unchanged). ASK14, CB14, and CY14 reproduce the nshmp-lib reference results. BSSA14 uses
  the published coefficients, which differ from the rounded nshmp-lib coefficients (and its `f_4`
  at 0.05, 0.075, 0.5, and 1.5 s), so it differs from nshmp-lib by up to about 0.04 in the natural
  log of the median at these periods for soft sites and less than 0.001 otherwise.
  `mechanism_from_rake` gives the nshmp-lib style of faulting of a rake.
- Added: `NgaEastSeeds`, the USGS NGA-East seed model logic tree of nshmp-lib (14 weighted seed
  models; `NGA_EAST_SEEDS_2018`, `_2023`, and `_2026` with the `adjusted` and `cpa` variants), which
  the current CONUS NSHM combines with `NgaEast` for stable crust. The model reproduces the
  nshmp-lib reference results. `dist_jb` is required (for the Shahjouei and Pezeshk (2016) seed).
  PGV is not provided (nshmp-lib uses the Abrahamson and Bhasin (2020) conditional PGV model for
  the seeds without PGV tables, which is not implemented).
- Added: `NgaEast`, the USGS nshmp-lib NGA-East versions used by the NSHM: "2018", "2023", and
  "2026" (default), with the `adjusted` (2023 short-period adjustment) and `cpa` (Chapman and Guo
  (2021) Gulf and Atlantic coastal plain amplification) variants. The 2023 and 2026 versions update
  the site amplification. The standard deviation combines the panel and EPRI models as in nshmp-lib
  (the square root of the weighted variances), and the model reproduces the nshmp-lib reference
  results.
- Added: `depth_sed` scenario value, the thickness of the coastal plain sediments (km), used by
  `NgaEast` with `adjusted` or `cpa`.
- Added: the final NGA-Subduction models of Abrahamson and Gülerce (2020, `AbrahamsonGulerce2020`),
  Kuehn et al. (2020, `KuehnEtAl2020`), and Parker et al. (2020, `ParkerEtAl2020`), following the
  USGS nshmp-lib implementations used by the current NSHM. `event_type` ("interface"/"intraslab")
  and `region` ("global", "alaska", "cascadia", "prvi") are scenario values, and constructor options
  (e.g., `basin`, `adjusted`, `ak_adjusted`, `seattle_basin`, `m9`, `epistemic`) select the
  nshmp-lib variants. The epistemic branches are collapsed as in nshmp-lib. Each model reproduces the
  nshmp-lib reference results. As in nshmp-lib, the Cascadia basin term is applied in all regions
  when `depth_2_5` is given, and the Kuehn et al. (2020) epistemic uncertainty does not depend on
  distance.
- Added: `ln_pga` property, the natural logarithm of PGA, which avoids computing `np.log(pga)`.
- Added: `ims` argument for all ground motion models (`GroundMotionModel` subclasses) to compute only
  the requested intensity measures: "pga", "pgv", individual spectral periods such as "psa_1p000" (1.0 s),
  "psa_ngawest2_21" (the 21 NGA-West2 comparison periods), and "psa_all" (all periods).
  `periods` and `spec_accels` contain only the computed periods, and `psa_ims` gives their names.
- Changed: models are organized into three packages, chosen by the independent
  variable of what they predict — shear strain (`soil_curves`), frequency
  (`fourier_spectrum`), or period / a scalar conditioned on a `Scenario`
  (`ground_motion`). The top-level names are unchanged. This replaces the nine
  category subpackages introduced earlier in this release, which were never
  published.
- Changed: `pint` is now a runtime dependency, and `soil_curves` is imported
  lazily. `import pygmm` previously failed outright on any install without pint.
- Changed: documentation reorganized to the Scientific Python layout
- Changed: `response_spectrum()` and `duration_model()` accept models of
  vectorized scenarios. For N scenarios, `ResponseSpectrum.spec_accels` has
  shape (N, periods) and the `Duration` values have shape (N,), and `ims`
  limits the response spectrum to the requested periods. The `Duration` values
  are typed as a float or an array (`contracts.FloatOrArray`).
- Changed: Idriss (2014) accepts arrays of scenario values (`mag`, `dist_rup`, `v_s30`, `mechanism`).
  Arrays of N scenarios give results with shape (N, periods), e.g., `pga` with shape (N,).
  Scalar scenarios give the same results as before.
- Changed: `NumericParameter` and `CategoricalParameter` check arrays of values, with one warning
  for the values outside of the limits or options.
- Changed: Boore, Stewart, Seyhan, and Atkinson (2014) accepts arrays of scenario values (`mag`,
  `dist_jb`, `v_s30`, `depth_1_0`, `mechanism`, `region`) and the `ims` argument. Scalar scenarios give
  the same results as before.
- Changed: Abrahamson, Silva, and Kamai (2014), Campbell and Bozorgnia (2014), and Chiou and Youngs
  (2014) accept arrays of scenario values and the `ims` argument. Values that are estimated when not
  provided (e.g., `depth_tor`, `width`, `depth_1_0`, `depth_2_5`, `depth_hyp`) are estimated for each
  scenario. Scalar scenarios give the same results as before.
- Changed: Abrahamson, Gregor, and Addo (2016) and Coppersmith and Bommer (2014) accept arrays of
  scenario values (`mag`, `dist_rup`, `dist_hyp`, `depth_hyp`, `v_s30`, `event_type`,
  `tectonic_region`) and the `ims` argument, which follows the other keyword arguments. `dist_hyp`
  and `depth_hyp` may be _None_ when no events are intraslab, and a missing distance or depth that is
  needed raises a `ValueError`.
- Changed: Akkar, Sandikkaya, and Bommer (2014), Derras, Bard, and Cotton (2014), and Hermkes, Kuehn,
  and Riggelsen (2014) accept arrays of scenario values and the `ims` argument. Derras et al. (2014)
  raises a `ValueError`, instead of a `KeyError`, for an unsupported mechanism.
- Changed: Atkinson and Boore (2006), Campbell (2003), Pezeshk et al. (2011), and Tavakoli and
  Pezeshk (2005) accept arrays of `mag` and `dist_rup` (and `v_s30` for Atkinson and Boore (2006))
  and the `ims` argument.
- Changed: `calc_width`, `calc_depth_tor`, `calc_depth_1_0`, and `calc_site_term` of Abrahamson et al.
  (2014); `calc_depth_2_5`, `calc_width`, `calc_depth_hyp`, `calc_depth_bor`, and `calc_site_term` of
  Campbell and Bozorgnia (2014); and `calc_depth_tor` and `calc_site_term` of Chiou and Youngs (2014)
  accept arrays.
- Changed: the duration models of Abrahamson and Silva (1996), Kempton and Stewart (2006), Afshari
  and Stewart (2016), and Pinilla-Ramos et al. (2023, 2024) accept arrays of scenario values. Arrays
  of N scenarios give durations and standard errors with shape (N,) (fields of the record arrays for
  Kempton and Stewart (2006) and Afshari and Stewart (2016)), and `AbrahamsonSilva1996.interp`
  returns shape (N, len(nias)), or (N, len(stds), len(nias)) with `stds`. Scalar scenarios give the
  same results as before. `AfshariStewart2016.calc_depth_1_0` accepts arrays of `v_s30` and
  `region`.
- Changed: `pygmm.ground_motion.pinilla_ramos_et_al_2024.duration_model` uses all of the values of array inputs. It
  previously used only the first value.
- Changed: `ChiouYoungs2014.calc_depth_1_0` accepts an array of regions.
- Fixed: `AbrahamsonBhasin2020.__init__` ended with `return ln_mean, ln_std`, so
  the model could not be constructed at all. Results are now exposed as
  `pgv`/`ln_pgv`/`ln_std`/`tau`/`phi`, matching the other conditional models.
- Fixed: `ChiouYoungs2014` no longer silently loses the `tau`/`phi` return added
  in 0.6.6; a duplicate copy of the module had reverted it.
- Fixed: `Stafford2017` referenced coefficient keys absent from its data file in
  one of two duplicate copies. Nine modules existed in duplicate; the stale
  copies are removed.
- Fixed: parameter-limit warnings printed `{self.min}`/`{self.max}` literally.
- Fixed: `AbrahamsonSilva1996.interp` and `calc_ln_dur_incr` no longer modify the `nias` array that
  is passed in, and accept integer values.
- Fixed: Boore et al. (2014) used the global relation for the reference depth to 1.0 km/s for the
  Japan region. It now uses the Japan relation, which changes the basin term at periods of 0.65 s
  and longer for `region="japan"` with `depth_1_0` specified.
- Fixed: Chiou and Youngs (2014) used the Boore et al. (2014) reference velocities (570.94 and
  412.39 m/s) in the average depth to 1.0 km/s (`calc_depth_1_0`), instead of the 571 and 412 m/s
  of Equations 1 and 2 in the paper. This changes the basin term by up to 0.002 in ln when
  `depth_1_0` is given. Boore et al. (2014) keeps its own 570.94 and 412.39 m/s.
- Fixed: Abrahamson et al. (2014) used a V1 of 801.2 m/s for periods of 3 s and longer, from
  extending the log-linear transition, instead of the 800 m/s of Equation 9. This changes the
  site term by up to 0.0014 in ln for `v_s30` above 800 m/s.
- Fixed: Pezeshk et al. (2011) used the standard deviation in log10 units as a natural log standard
  deviation. It is now converted to natural log units (multiplied by ln(10)), as in OpenQuake,
  e.g., 0.58 instead of 0.25 for PGA at M6.
- Fixed: Atkinson and Boore (2006) used the hard-rock coefficients only for `v_s30` of zero, so
  hard-rock sites used the B/C coefficients with the site amplification extrapolated to their
  velocity. Sites with `v_s30` of 2000 m/s or greater now use the hard-rock coefficients without site
  amplification, as in the paper and OpenQuake. At 2000 m/s, PGA increases by a factor of about 1.4
  to 2.1 and PSA at 1 s by about 1.5.
- Fixed: Atkinson and Boore (2006) divided PGV and PGD by g along with PGA and PSA, so `pgv` and
  `pgd` were 980.665 times too small. They are now in cm/sec and cm, as for the other models (e.g.,
  PGV of 44 cm/sec for M7 at 10 km on B/C).
- Fixed: Atkinson and Boore (2006) used the standard deviation of 0.30, which is in log10 units, as
  a natural log standard deviation. It is now 0.30 ln(10) = 0.691, as in nshmp-haz and OpenQuake.
- Fixed: Atkinson and Boore (2006) left out the "+ b_2" in the nonlinear site coefficient for
  180 < Vs30 <= 300 m/s (Eq. 8), so the site term was discontinuous at 180 and 300 m/s. This changes
  the results for Vs30 in this range, e.g., PGA decreases by about 28% for M7 at 10 km with
  Vs30 = 250 m/s.
- Fixed: `PinillaRamosEtAl2024.d5x_median` scaled the D5-75 median by approximate factors (e.g.,
  0.25 for D5-10 and 1.80 for D5-95), and `d5x_sigma` returned the D5-75 standard deviation for all
  thresholds. Both now use the conditional model of `duration_for_energy`. For example, D5-10 is
  about 4% of D5-75 instead of 25%, and D5-95 is about 2.0 to 2.6 times D5-75 instead of 1.8.
  Thresholds from "D5-10" to "D5-95" in steps of 5% are supported.
- Fixed: Pinilla-Ramos et al. (2023) had a0 = 1 instead of 0 for the 0.75 row of the conditional
  model, which doubled the durations from `duration_for_energy` for energies between about 0.725
  and 0.775 other than 0.75. These energies now give the D5-75 durations.
- Fixed: `AbrahamsonSilva1996.interp` added the increment for the 5 to 75% duration twice, which
  made the interpolated durations about 0.9% too large. `interp(0.75)` now equals `duration`, and
  the results match the reference spreadsheet to within its rounding.
- Fixed: the `adjust_c4` argument of Abrahamson, Gregor, and Addo (2016) was stored but not used. It
  is now added to C_4 (10 km) in the finite-fault distance term. The default of 0 gives the same
  results as before.
- Fixed: Campbell and Bozorgnia (2014) did not set PSA at periods shorter than 0.25 s to PGA when
  it is less than PGA, as specified by the model. This increases the short-period PSA for sites with
  strong nonlinear site response, e.g., by up to about 40% at 0.05 to 0.075 s for Vs30 = 150 m/s.
  The standard deviations are unchanged.
- Fixed: the Hermkes, Kuehn, and Riggelsen (2014) model data in the repository was a Dropbox web page
  instead of the model data, because the download URL used `dl=0`. The data (164 MB, Git LFS) is
  restored, the download uses `dl=1` and checks that the downloaded file is the model data, and a
  Git LFS pointer file (a clone without Git LFS) is replaced by the download. The data is not
  included in the package and is downloaded when the model is first imported.
- Fixed: Campbell and Bozorgnia (2014) used the global anelastic attenuation coefficient
  (Δc20 = 0) for `region="china"` because the region was compared with a list. It now uses the
  China coefficient (Δc20,CH), which increases the response at distances greater than 80 km
  (e.g., PGA by about 30% at 150 km for M6.5 and Vs30 = 760 m/s).
- Fixed: Afshari and Stewart (2016) applied the maximum basin term (a depth differential of 200 m)
  for `depth_1_0=np.nan`. An unknown depth (NaN) now has no basin term, as for `depth_1_0=None`, which
  allows arrays of scenarios with some unknown depths.
- Fixed: Afshari and Stewart (2016) passed the mechanism as the basin region, so the California
  relation for the mean depth to 1.0 km/s was always used. It now has an optional `region`
  parameter ("california", "global", or "japan"; default "california") and uses the Japan
  relation for `region="japan"`.
- Fixed: Afshari and Stewart (2016) computed the basin depth differential in km instead of m,
  which made the basin term about 1,000 times too small. Figure 12 of the paper shows that c_5
  applies to the differential in m (Table 2 lists c_5 in 1/km), consistent with the 200 m limit.
- Fixed: Afshari and Stewart (2016) used b_1 = 6.188 instead of 6.579 for the D20-80 duration of
  reverse earthquakes (Table 1 of the paper).
- Deprecated: `pygmm.models`. It mixed response-spectrum, duration, FAS and
  correlation models with no way to tell them apart, and 7 of its 20 entries had
  no `spec_accels`. Use `find_models()`.
- Removed: `pygmm.contracts.FaultDisplacement`, `CptSounding`,
  `SoilBehaviorProfile` and `LiquefactionTriggering`, which had no producer and
  no consumer
- Removed: `pygmm.tools` (an empty stub) and `pygmm.types` (renamed `_types`; it
  shadowed the stdlib `types` module, breaking any interpreter started with the
  package directory as the working directory)

## 0.8.0 (2025-07-24)

- Added: Pinilla-Ramos et al. (2023) model for duration of crustal earthquakes
- Added: Pinilla-Ramos et al. (2024) model for duration of subduction earthquakes
- Added: Stafford (2017) model for FAS correlation

## 0.7.3 (2025-03-12)

- Fixed Bayless and Abrahamson (2018) correlation model.

## 0.7.1 (2025-03-05)

- Added compatibility with numpy 2.0

## 0.7.0 (2024-04-24)

- Added: Abrahamson and Bhasin (2020)
- Changed to Hatch build system

## 0.6.6 (2023-12-11)

- Added: Return tau and phi in the standard deviation calculations

## 0.6.5 (2022-09-16)

- Added: Afshari and Stewart (2016) duration model
- Added: Kempton and Stewart (2006) duration model

## 0.6.4 (2022-01-24)

- Added: Bayless and Abrahamson (2019)

## 0.6.3 (2021-12-08)

- Fixed: error in ASK14 on a7 term

## 0.6.2 (2021-10-19)

- Changed: Move site amplification to static functions on some GMPEs

## 0.6.1 (2020-06-03)

- Added Coppersmith & Bommer (2014) model for Hanford
- Factored tests

## 0.6.0 (2019-08-12)

- Added Abrahamson, Gregor, Addo (2014)
- Added Abrahamson & Gulerce (2011)
- Added conditional mean spectra models.
- Added Scenario objects.
- Added typing for all classes.

## 0.4.0 (2016-04-08)

- Added Hermkes et al. (2014).
- Improved documentation.
- Added Baker & Jayaram (2008), Kishida (2017)

## 0.3.2 (2016-03-30) {#section-1}

- Nothing changed yet.

## 0.3.1 (2016-03-30) {#section-2}

- First release on PyPI.
