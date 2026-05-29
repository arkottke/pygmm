# Plan: Add Macedo, Abrahamson, & Liu (2021) CAV Model to pyGMM

## Paper

- **Citation:** Macedo, J., Abrahamson, N., and Liu, C. (2021). *New Scenario-Based
  Cumulative Absolute Velocity Models for Shallow Crustal Tectonic Settings.*
  Bulletin of the Seismological Society of America, **111**(1), 157–172.
- **DOI:** `10.1785/0120190321`
- **Local copy:** `bssa-2019321.1.pdf`

## Summary of what the paper provides

Two related models for cumulative absolute velocity (CAV, units of m/s) for shallow
crustal earthquakes, calibrated on the NGA-West2 database (~14,000 recordings,
287 events, M 3.0–7.9):

1. **Conditional Ground-Motion Model (CGMM)** — Eq. (12). Predicts `ln(CAV)`
   conditioned on PGA, Mw, Vs30, Rrup plus a hanging-wall term:

   ```
   ln(CAV) = c1 + c2·ln(PGA) + c3·Mw + c4·ln(Vs30) + c5·ln(Rrup)
            + c6·F_HW · T1(dip) · T2(Mw) · T5(R_JB)
   ```

   with `c1=1.79, c2=0.67, c3=0.57, c4=-0.47, c5=-0.0026, c6=0.17`,
   between-event τ = 0.17, within-event φ = 0.26, total σ = 0.35.

   `T1(dip)`, `T2(Mw)`, `T5(R_JB)` are the ASK14 hanging-wall tapers (Eqs. 13a–c).

2. **Scenario-based models** — Eq. (14). Same functional form but with the
   conditioning PGA replaced by the median PGA predicted by an NGA-West2 GMM
   (ASK14, BSSA14, CB14, CY14, or I14). The total σ is obtained via propagation
   of errors (Eqs. 16–17):

   ```
   σ_lnCAV = sqrt(0.097 + 0.4541 · σ²_lnPGA)
   ```

## Where this fits in pyGMM

pyGMM groups models by author/year using a snake_case module name and a
CamelCase class. The model is a scalar IM predictor (like AB20 PGV) rather than
a spectral GMM, so it follows the simpler single-output pattern.

Files to add or touch (mirroring existing structure):

| Path | Purpose |
|---|---|
| `src/pygmm/macedo_abrahamson_liu_2021.py` | Module containing the new class |
| `src/pygmm/__init__.py` | Import + register in `models` / `__all__` |
| `tests/test_macedo_abrahamson_liu_2021.py` | Pytest suite + reference values |
| `tests/data/macedo_abrahamson_liu_2021.csv` *(optional)* | Reference CAV values used by tests |
| `docs/references.bib` | BibTeX entry with key `macedo21` |
| `docs/models.rst` | One-line entry referencing the new model |
| `HISTORY.md` | Changelog note |

## Class design

Follow the conventions seen in `abrahamson_bhasin_2020.py`,
`stafford_2017.py`, and `pinilla_ramos_et_al_2024.py`:

```python
class MacedoAbrahamsonLiu2021(model.Model):
    NAME = "Macedo, Abrahamson, & Liu (2021)"
    ABBREV = "MAL21"

    PARAMS = [
        model.NumericParameter("mag",        True,  3.0, 8.0),
        model.NumericParameter("dist_rup",   True,  0,   300),
        model.NumericParameter("dist_jb",    True,  0,   300),
        model.NumericParameter("v_s30",      True,  150, 1500),
        model.NumericParameter("dip",        True,  15,  90),
        model.NumericParameter("on_hanging_wall", False),
        # PGA in g — required only for the conditional flavor
        model.NumericParameter("pga",        False, 1e-4, 3.0),
    ]

    # Coefficients from Table 1
    COEFFS = dict(c1=1.79, c2=0.67, c3=0.57, c4=-0.47,
                  c5=-0.0026, c6=0.17)
    TAU, PHI = 0.17, 0.26
    SIGMA = (TAU**2 + PHI**2) ** 0.5  # ≈ 0.31; paper reports 0.35

    SUPPORTED_PGA_GMMS = {
        "ASK14": AbrahamsonSilvaKamai2014,
        "BSSA14": BooreStewartSeyhanAtkinson2014,
        "CB14":  CampbellBozorgnia2014,
        "CY14":  ChiouYoungs2014,
        "I14":   Idriss2014,
    }
```

### Two evaluation modes

Expose both flavors through a single class with an optional `pga_model` argument:

- **Conditional mode** — user supplies `scenario["pga"]` (and optionally `pga_ln_std`).
  Returns `ln(CAV)` directly from Eq. (12) with sigma from Table 1.
- **Scenario mode** — user passes `pga_model="ASK14" | "BSSA14" | "CB14" | "CY14" | "I14"`
  (or a `model.Model` subclass). The class internally builds the PGA model from the
  same scenario, plugs `PGA_med` and `σ_lnPGA` into Eqs. (14) and (17).

### Properties to expose

Mirror other scalar-IM models:

- `ln_cav` — natural log of CAV (m/s)
- `cav` — median CAV (m/s)
- `ln_std` / `tau` / `phi` — aleatory components
- `cav_plus_sigma`, `cav_minus_sigma` — convenience accessors

### Helper functions (private)

- `_t1_dip(dip)` — Eq. (13a)
- `_t2_mag(mag)` — Eq. (13b)
- `_t5_rjb(dist_jb)` — Eq. (13c)
- `_hanging_wall_term(scenario)` — combines flag, T1, T2, T5
- `_check_pga_model(...)` — validates the chosen PGA backbone

## Validation strategy

The paper contains no published lookup table of test values, so reference
predictions must be generated and locked in:

1. Pick a small grid of scenarios that exercises every branch:
   - Magnitudes: 5.0, 6.0, 7.0, 7.5
   - Distances: 5, 30, 100 km
   - Vs30: 270, 425, 760 m/s
   - One HW and one FW site (dip = 45°, 90°)
2. Compute CAV by hand from Eq. (12) using the published coefficients —
   commit these to `tests/data/macedo_abrahamson_liu_2021.csv`.
3. For scenario mode, generate the PGA from each of the five NGA-West2 GMMs
   already in pyGMM and assert `ln_cav` matches a value computed with the
   same coefficients, plus σ computed from Eq. (17).
4. Compare overall trends qualitatively against Figures 8–12 in the paper
   (distance scaling, magnitude scaling, σ vs distance) in a notebook-style
   sanity check; not a unit test.

## Tests

Mirror `tests/test_pinilla_ramos_et_al_2024.py`. Cover:

- Initialization with required parameters.
- Conditional mode: pass `pga`, verify `cav`, `ln_std`, `tau`, `phi`.
- Scenario mode for each of the 5 PGA backbones, verifying `ln_cav` and the
  Eq. (17) σ.
- Hanging-wall on/off symmetry: F_HW = 0 must zero out the c6 term.
- Magnitude taper `T2`: assert continuity at M = 5.5 and 6.5.
- Dip taper `T1`: assert plateau for dip ≤ 30°.
- R_JB taper `T5`: assert zero for R_JB ≥ 15 km.
- Out-of-range parameter warnings via `model.Model.PARAMS` validation.
- Vector inputs should broadcast.

## Documentation

- Add a BibTeX entry to `docs/references.bib`:
  ```bibtex
  @article{macedo21,
    author = {Macedo, J. and Abrahamson, N. and Liu, C.},
    title  = {New Scenario-Based Cumulative Absolute Velocity Models
              for Shallow Crustal Tectonic Settings},
    journal= {Bulletin of the Seismological Society of America},
    volume = {111}, number = {1}, pages = {157--172}, year = {2021},
    doi    = {10.1785/0120190321},
  }
  ```
- Module docstring uses `:cite:`macedo21`` so Sphinx picks it up
  (same pattern as `abrahamson20`, `pinilla-ramos24`).
- Add bullet to `docs/models.rst` and one-line entry to `HISTORY.md`.

## Decisions

1. One class with `pga_model=None` meaning "conditional".
2. **PGA units in scenario.** Existing models pass PGA in g — keep that convention.
3. **HW handling.** Some pyGMM models infer HW from `dist_x`; others use an
   explicit flag (`on_hanging_wall`). The paper uses an explicit flag, so use
   that, but auto-set it from `dist_x > 0` when missing.
4. **Sigma reconciliation.** `sqrt(0.17² + 0.26²) ≈ 0.31`, but the paper
   states σ_total = 0.35 and Eq. (17) uses 0.097 = 0.31². Use the published
   values directly and add a comment noting the apparent rounding.
5. **VS30 range.** Paper does not state explicit bounds; use NGA-West2 typical
   range (150–1500 m/s) and emit a warning outside.

## Step-by-step implementation checklist

1. Read the paper sections "Conditional CAV Models" and "Conversion of the
   Conditional Models to Scenario-Based Models" once more for any nuances
   missed (footnotes on Table 1, supplemental material).
2. Create `src/pygmm/macedo_abrahamson_liu_2021.py` with helpers, class,
   conditional and scenario evaluation paths.
3. Register the class in `src/pygmm/__init__.py` (`from ... import ...`,
   add to `__all__`, append to `models`).
4. Hand-compute a small lookup table of expected `ln(CAV)` values and store
   as `tests/data/macedo_abrahamson_liu_2021.csv`.
5. Write `tests/test_macedo_abrahamson_liu_2021.py` covering the cases above.
6. Run `pytest tests/test_macedo_abrahamson_liu_2021.py` and full suite.
7. Update `docs/references.bib`, `docs/models.rst`, `HISTORY.md`.
8. Run `ruff`/`black` (per `pyproject.toml`) and `pytest` once more.
9. Open PR referencing the DOI.
