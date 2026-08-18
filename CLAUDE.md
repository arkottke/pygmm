# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environments

Two, and the difference matters:

- **`/home/albert/envs/py312`** (Python 3.12.9) — a shared environment with
  pygmm, **pyrvt and pystrata all installed editable**. Use this for anything
  touching cross-package interop; it is the only place the
  `tests/test_fourier_spectrum_interop.py` cases actually run rather than
  being skipped by `importorskip`.
- **`.venv`** (uv-managed) — synced from `[dependency-groups] dev`. Use via
  `uv run`. pyrvt/pystrata are absent here, so the interop module is skipped
  (645 passed instead of 654).

Coverage is **opt-in** (`--cov=pygmm`). The old `pytest.ini` hard-coded it in
`addopts`, which silently overrode `[tool.pytest.ini_options]` in
`pyproject.toml` and made pytest fail outright wherever pytest-cov was absent.
All pytest config now lives in `pyproject.toml`.

## Commands

```bash
# Run all tests (use the py312 env for interop coverage)
/home/albert/envs/py312/bin/python -m pytest tests/
uv run pytest                      # equivalent, minus the interop tests

# A single test file
/home/albert/envs/py312/bin/python -m pytest tests/test_registry.py

# Coverage
uv run pytest --cov=pygmm --cov-report=html

# Lint / format
uv run ruff check src/ tests/
uv run ruff format src/ tests/

# Type check (config and staged allowlist in pyproject.toml)
uv run mypy

# Build docs
uv run --group docs make -C docs html
```

## Architecture

pyGMM provides a unified interface to ground motion prediction equations (GMPEs).

### Core pattern

**No model defines `__call__`.** Construct with a `Scenario` and read results
off attributes:

```python
model = pygmm.ChiouYoungs2014(scenario)
model.periods, model.spec_accels, model.ln_stds   # median PSA in g
model.pga, model.pgv                              # where the model provides them
model.response_spectrum()                         # -> contracts.ResponseSpectrum
```

`Model` is a scenario filter/validator; `GroundMotionModel` adds the "one flat
`_ln_resp` vector sliced by `INDICES_PSA` / `INDEX_PGA` / `INDEX_PGV` /
`INDEX_PGD`" convention. That is why one model instance commonly provides
several intensity measures at once.

### Module layout

- `src/pygmm/model.py` — `Model`, `GroundMotionModel`, `Scenario`, `Parameter`, `Coefficients`.
- `src/pygmm/registry.py` — `@register`, `ModelInfo`, `find_models()`, `get_model()`. Capability-based discovery; replaces the deprecated `pygmm.models` list.
- `src/pygmm/contracts.py` — Output dataclasses (`ResponseSpectrum`, `FourierSpectrum`, `Duration`, `NonlinearSoilCurves`, `VelocityProfile`) plus the `Supports*` **Protocols** producers satisfy structurally. pystrata and pyrvt duck-type against these and take no runtime dependency on pygmm.
- `src/pygmm/kamai_2016.py` — `kea16_profile`; sole producer of `VelocityProfile`, so it has no package of its own.
- `src/pygmm/__init__.py` — Public API.

### The three model packages

**The independent variable decides where a model lives.** This rule is total
and requires no judgment, which is the point: the previous nine category
subpackages demanded a four-axis classification and consequently misfiled
`AbrahamsonBhasin2020` (a PGV model, filed under `correlation/`),
`GulerceAbrahamson2011` (a V/H ratio model, filed under `response_spectrum/`)
and `Stafford2017`.

| Varies with | Package | Contents |
|---|---|---|
| shear strain | `soil_curves/` | G/Gmax and damping curves (Darendeli, Menq, Kishida, Wang, Alemu, Rollins) |
| frequency | `fourier_spectrum/` | FAS models (BA19, Stafford 2022, SourceTheoryModel) and inter-frequency correlation (BA18, Stafford 2017) |
| period, or a scalar from a `Scenario` | `ground_motion/` | response spectra, peak parameters, duration, Arias, CAV, V/H ratio, inter-period correlation |

Anything else stays at the top level until a second producer of the same
quantity earns it a package.

A model's **capabilities are registry metadata, not its location** — a
directory can express only one classification, and `AtkinsonBoore2006` emits
PSA + PGA + PGV + PGD from a single response vector.

`ground_motion/__init__.py` deliberately does not import
`hermkes_kuehn_riggelsen_2014`: that module attempts a 430 kB download at
import time when its data file is absent.

### Adding a new model

1. Create `src/pygmm/<package>/<author_year>.py`, choosing the package by the rule above.
2. Add `NAME`, `ABBREV`, and `PARAMS` class attributes. `ABBREV` is display-only and **not unique** (`CampbellBozorgnia2014` and `CoppersmithBommer2014` both use `CB14`).
3. Implement `__init__(self, scenario)`. Store results on `self`; do not `return` from `__init__`.
4. Apply `@register(...)`. PSA/PGA/PGV/PGD are derived from the `INDEX_*` attributes for `GroundMotionModel` subclasses; anything else must pass `provides=(...)`.
5. Export from the package `__init__.py` and, if user-facing, from `src/pygmm/__init__.py` and `__all__`.
6. Regenerate the table in `docs/user-guide/models.rst` (it is generated from `find_models()`).

`tests/test_registry.py` instantiates every registered model and checks each
declared capability resolves to a real, finite attribute — so a model that is
registered but broken, or exported but unregistered, fails the suite.

## Version support policy

This repo follows [SPEC 0](https://scientific-python.org/specs/spec-0000/):

- **Python** — support the minors released in the last **36 months**.
- **Core dependencies** — support the minors released in the last **24 months**.

Floor as of 2026-08-03: **Python 3.12**; CI runs 3.12, 3.13, 3.14.

The dated drop schedule lives in `pyproject.toml` next to `[project.urls]`.
When advancing a floor, update all four in step: `requires-python`, the
dependency lower bounds, `[tool.ruff] target-version`, and the CI matrix.
Then `rm uv.lock && uv lock` — a stale lock is what left this repo pinned to a
scipy without wheels for the newest Python.
