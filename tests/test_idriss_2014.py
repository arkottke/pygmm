"""Test the vectorized Idriss (2014) model against scalar scenarios."""

import itertools
import warnings

import numpy as np
import pytest

import pygmm

MAGS = [5.0, 6.0, 6.75, 6.76, 7.5, 8.2]
DISTS = [0.0, 10.0, 75.0, 150.0]
V_S30S = [450.0, 760.0, 1200.0]
MECHANISMS = ["SS", "RS"]
PERIODS = [0.05, 0.3, 1.0, 2.5]


def scalar_results(mag, dist_rup, v_s30, mechanism):
    m = pygmm.Idriss2014(
        pygmm.Scenario(mag=mag, dist_rup=dist_rup, v_s30=v_s30, mechanism=mechanism)
    )
    return {
        "pga": m.pga,
        "ln_std_pga": m.ln_std_pga,
        "spec_accels": m.spec_accels,
        "ln_stds": m.ln_stds,
        "interp_spec_accels": m.interp_spec_accels(PERIODS),
        "interp_ln_stds": m.interp_ln_stds(PERIODS),
    }


@pytest.fixture(scope="module")
def grid():
    rows = list(itertools.product(MAGS, DISTS, V_S30S, MECHANISMS))
    mag, dist_rup, v_s30, mechanism = (np.array(c) for c in zip(*rows))
    return rows, pygmm.Scenario(
        mag=mag, dist_rup=dist_rup, v_s30=v_s30, mechanism=mechanism
    )


@pytest.mark.parametrize(
    "key",
    [
        "pga",
        "ln_std_pga",
        "spec_accels",
        "ln_stds",
        "interp_spec_accels",
        "interp_ln_stds",
    ],
)
def test_vectorized_matches_scalar(grid, key):
    rows, scenario = grid
    m = pygmm.Idriss2014(scenario)
    if key.startswith("interp"):
        actual = getattr(m, key)(PERIODS)
    else:
        actual = getattr(m, key)
    expected = np.array([scalar_results(*row)[key] for row in rows])
    assert actual.shape == expected.shape
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, expected, rtol=1e-12)
    else:
        np.testing.assert_array_equal(actual, expected)


def test_vectorized_shapes(grid):
    rows, scenario = grid
    m = pygmm.Idriss2014(scenario)
    n = len(rows)
    n_periods = len(m.periods)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, n_periods)
    assert m.ln_stds.shape == (n, n_periods)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged():
    m = pygmm.Idriss2014(
        pygmm.Scenario(mag=6.5, dist_rup=20.0, v_s30=760.0, mechanism="SS")
    )
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (len(m.periods),)
    assert m.ln_stds.shape == (len(m.periods),)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_scalars_broadcast_with_arrays():
    # A scalar v_s30 and mechanism with arrays of magnitude and distance
    mag = np.array([5.5, 6.5, 7.5])
    dist_rup = np.array([5.0, 20.0, 80.0])
    m = pygmm.Idriss2014(
        pygmm.Scenario(mag=mag, dist_rup=dist_rup, v_s30=760.0, mechanism="RS")
    )
    expected = [
        scalar_results(mg, d, 760.0, "RS")["pga"] for mg, d in zip(mag, dist_rup)
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_multidimensional_inputs():
    # A grid of magnitudes and distances keeps its shape
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [5.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    m = pygmm.Idriss2014(
        pygmm.Scenario(mag=mag, dist_rup=dist_rup, v_s30=760.0, mechanism="SS")
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, len(m.periods))
    assert m.pga[2, 1] == scalar_results(7.5, 20.0, 760.0, "SS")["pga"]


def test_invalid_mechanism_entries_use_default():
    # As for a scalar scenario, an unsupported mechanism is replaced by the default (SS)
    mechanism = np.array(["SS", "NS", "RS", "U"])
    with pytest.warns(UserWarning, match="2 of 4 values"):
        m = pygmm.Idriss2014(
            pygmm.Scenario(
                mag=np.full(4, 6.5),
                dist_rup=np.full(4, 20.0),
                v_s30=760.0,
                mechanism=mechanism,
            )
        )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = [scalar_results(6.5, 20.0, 760.0, mech)["pga"] for mech in mechanism]
    np.testing.assert_array_equal(m.pga, expected)


def test_out_of_range_values_warn_once():
    v_s30 = np.array([300.0, 400.0, 760.0])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        pygmm.Idriss2014(
            pygmm.Scenario(mag=6.5, dist_rup=20.0, v_s30=v_s30, mechanism="SS")
        )
    messages = [str(w.message) for w in caught if "v_s30" in str(w.message)]
    assert len(messages) == 1
    assert "2 of 3 values" in messages[0] and "minimum of 300.0" in messages[0]


def test_lists_are_accepted():
    m = pygmm.Idriss2014(
        pygmm.Scenario(
            mag=[6.0, 7.0], dist_rup=[10.0, 50.0], v_s30=760.0, mechanism=["SS", "RS"]
        )
    )
    expected = [
        scalar_results(6.0, 10.0, 760.0, "SS")["pga"],
        scalar_results(7.0, 50.0, 760.0, "RS")["pga"],
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        pygmm.Idriss2014(
            pygmm.Scenario(
                mag=np.array([6.0, 7.0]),
                dist_rup=np.array([10.0, 20.0, 30.0]),
                v_s30=760.0,
                mechanism="SS",
            )
        )


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return pygmm.Scenario(
        mag=rng.uniform(5.0, 8.2, n),
        dist_rup=rng.uniform(0.0, 150.0, n),
        v_s30=rng.uniform(450.0, 1200.0, n),
        mechanism=rng.choice(["SS", "RS"], n),
    )


@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(ims):
    s = vector_scenario()
    full = pygmm.Idriss2014(s)
    pga_only = pygmm.Idriss2014(s, ims=ims)
    np.testing.assert_array_equal(pga_only.pga, full.pga)
    np.testing.assert_array_equal(pga_only.ln_std_pga, full.ln_std_pga)
    # Only one period is computed
    assert pga_only._ln_resp.shape == (500, 1)


def test_pga_only_scalar_matches_default():
    s = pygmm.Scenario(mag=6.5, dist_rup=20.0, v_s30=760.0, mechanism="RS")
    full = pygmm.Idriss2014(s)
    pga_only = pygmm.Idriss2014(s, ims=["pga"])
    assert pga_only.pga == full.pga
    assert pga_only.ln_std_pga == full.ln_std_pga
    assert isinstance(pga_only.pga, float)


@pytest.mark.parametrize("ims", [["psa_all"], ["pga", "psa_all"]])
def test_psa_matches_default(ims):
    s = vector_scenario()
    full = pygmm.Idriss2014(s)
    m = pygmm.Idriss2014(s, ims=ims)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels)
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds)
    np.testing.assert_allclose(
        m.interp_spec_accels(PERIODS), full.interp_spec_accels(PERIODS), rtol=1e-12
    )
    # PGA is the first spectral period (0.01 s) in Idriss (2014), so it is
    # available whenever the 0.01 s spectral acceleration is computed
    np.testing.assert_array_equal(m.pga, full.pga)


@pytest.mark.parametrize(
    "attr", ["spec_accels", "ln_stds", "interp_spec_accels", "interp_ln_stds"]
)
def test_values_not_computed_raise(attr):
    m = pygmm.Idriss2014(vector_scenario(), ims=["pga"])
    with pytest.raises(ValueError, match="include 'psa_all'"):
        value = getattr(m, attr)
        if callable(value):
            value(PERIODS)


@pytest.mark.parametrize(
    "ims,match",
    [
        (["pgv"], "does not provide 'pgv'"),
        (["sa"], "not a valid intensity measure"),
        ([], "at least one intensity measure"),
    ],
)
def test_invalid_ims_raise(ims, match):
    with pytest.raises(ValueError, match=match):
        pygmm.Idriss2014(vector_scenario(), ims=ims)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    s = vector_scenario()
    m = pygmm.Idriss2014(s, ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_array_equal(m.ln_pga, pygmm.Idriss2014(s).ln_pga)
    np.testing.assert_allclose(m.ln_pga, np.log(m.pga), rtol=1e-14)


def test_ln_pga_scalar():
    m = pygmm.Idriss2014(
        pygmm.Scenario(mag=6.5, dist_rup=20.0, v_s30=760.0, mechanism="SS")
    )
    assert isinstance(m.ln_pga, float)
    assert np.exp(m.ln_pga) == m.pga


def test_ln_pga_available_with_psa():
    m = pygmm.Idriss2014(vector_scenario(), ims=["psa_all"])
    # PGA is the first spectral period in Idriss (2014), so it is available with PSA
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
