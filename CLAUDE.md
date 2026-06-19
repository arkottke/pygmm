# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Run all tests
/home/albert/envs/py312/bin/python -m pytest tests/

# Run a single test file
/home/albert/envs/py312/bin/python -m pytest tests/test_models.py

# Lint
uv run ruff check src/

# Format
uv run ruff format src/

# Build docs
uv run --group docs make -C docs html
```

## Architecture

pyGMM provides a unified interface to ground motion prediction equations (GMPEs).

### Core pattern

Every model is a callable class:
1. Instantiate with a :class:`~pygmm.model.Scenario` (earthquake + site parameters).
2. Call it; returns `(ln_mean, ln_std)` arrays over the model's spectral periods.

### Module layout

- `src/pygmm/model.py` — Base classes: `Model`, `GroundMotionModel`, `Scenario`, `Parameter`, `Coefficients`. All models inherit from `Model` (or its subclasses in `model.py`).
- `src/pygmm/contracts.py` — Frozen dataclasses `NonlinearSoilCurves` and `VelocityProfile` for pyStrata interop (no runtime dep on pystrata).
- `src/pygmm/__init__.py` — Public API. All user-facing names live here; the `models` list enumerates all GMPE classes.

### Subpackages (categorized by output type)

- `response_spectrum/` — NGA-West2 and international Sa models; `_base.py` holds shared base class
- `duration/` — Significant duration models (Afshari & Stewart 2016, Kempton & Stewart 2006, Pinilla-Ramos 2023/2024)
- `correlation/` — Inter-period and cross-model correlations (Baker-Jayaram, Bayless-Abrahamson, Stafford)
- `fourier_spectrum/` — FAS models: Bayless-Abrahamson 2019, Stafford 2022, SourceTheoryModel
- `soil_curves/` — Strain-dependent G/Gmax and damping curves (Darendeli, Menq, Kishida, Wang, Alemu, Rollins); `_base.py` has shared interpolation logic
- `velocity_profile/` — Kamai et al. (2016) reference Vs profiles (`kea16_profile`)

### Adding a new model

1. Create `src/pygmm/<author_year>.py` with a class inheriting from the appropriate base.
2. Add `NAME`, `ABBREV`, and `PARAMS` class attributes.
3. Implement `__init__(self, scenario)` and `__call__` (or rely on base class).
4. Import the class in `src/pygmm/__init__.py` and add to `__all__` and `models`.
5. Add an autosummary entry in `docs/user-guide/models.rst`.
