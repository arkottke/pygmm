"""Smoke-test every response-spectrum model by plotting it.

Previously this file parametrized over ``pygmm.models`` in a function named
``plot_model_with_param``. Since ``pyproject.toml`` sets
``python_functions = ["test_*"]``, pytest never collected it, so it has never
run -- and it could not have passed: 7 of the 20 entries in ``pygmm.models``
have no ``spec_accels``.

It now parametrizes over the registry, which returns exactly the models that
declare a PSA capability and take a ``Scenario``. Figures are written into the
pytest ``tmp_path`` (kept for the last few runs under ``/tmp/pytest-of-*``)
rather than a ``figures/`` directory in the working directory.
"""

import matplotlib

matplotlib.use("agg")  # NOQA
import matplotlib.pyplot as plt
import numpy as np
import pytest

import pygmm

# `dist`, `flag_hw` and `flag_meas` used to be listed here; none is a
# Scenario.KNOWN_KEYS member, so they raise. They are leftovers from an older
# API, invisible while this file went uncollected.
DEFAULT_PROPS = dict(
    depth_2_5=5,
    depth_bor=15,
    depth_hyp=9,
    depth_tor=5,
    dip=90.0,
    dist_jb=30.0,
    dist_rup=30.0,
    dist_x=30.0,
    dpp_centered=0,
    # Required by the subduction models (AGA16, CB14). Each model filters the
    # scenario down to its own PARAMS, so this is inert for the others.
    event_type="interface",
    mag=6,
    mechanism="SS",
    on_hanging_wall=False,
    v_s30=500.0,
    width=10,
)

PSA_MODELS = pygmm.find_models(provides="psa", input="scenario")


@pytest.mark.parametrize("info", PSA_MODELS, ids=lambda i: i.key)
@pytest.mark.parametrize(
    "key,values,label",
    [
        ("mag", [5, 6, 7], "Magnitude"),
        (["dist_rup", "dist_jb", "dist_x"], [10, 50, 100], "Distance (km)"),
        ("v_s30", [300, 650, 1000], "$V_{s30}$ (m/s)"),
    ],
    ids=lambda a: a[0],
)
def test_plot_model_with_param(info, key, values, label, tmp_path):
    props = dict(DEFAULT_PROPS)

    fig, ax = plt.subplots()
    for v in values:
        if isinstance(key, str):
            props[key] = v
        else:
            for k in key:
                props[k] = v
        m = info.cls(pygmm.Scenario(**props))

        # The registry claims this model provides PSA; hold it to that.
        assert len(m.periods) == len(m.spec_accels)
        assert np.all(np.isfinite(m.spec_accels))
        assert np.all(m.spec_accels > 0)

        ax.plot(m.periods, m.spec_accels, label=f"{v:g}")

    ax.set_xlabel("Period (s)")
    ax.set_xscale("log")
    ax.set_ylabel("5% Damped, Spectral. Accel. (g)")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1e1)
    ax.legend(loc="upper right", title=label, fontsize="x-small")
    ax.grid()
    fig.tight_layout()

    prefix = key if isinstance(key, str) else key[0]
    fig.savefig(tmp_path / f"{prefix}-{info.key}.png")
    plt.close(fig)
