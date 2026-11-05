#!/usr/bin/python
import gzip
import itertools
import json
import os
import warnings

import numpy as np
import pytest

from pygmm import DerrasBardCotton2014 as DBC14
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 1e-2

# Load the tests
fname = os.path.join(os.path.dirname(__file__), "data", "dbc14_tests.json.gz")
with gzip.open(fname, "rt") as fp:
    tests = json.load(fp)

testdata = [(t["params"], t["results"]) for t in tests]


@pytest.mark.parametrize("params,expected", testdata)
def test_spec_accels(params, expected):
    m = DBC14(Scenario(**params))
    np.testing.assert_allclose(
        m.interp_spec_accels(expected["periods"]),
        # Need to convert from m/sec to g
        np.array(expected["spec_accels"]) / m.GRAVITY,
        rtol=RTOL,
        err_msg="Spectral accelerations",
    )


@pytest.mark.parametrize("params,expected", testdata)
@pytest.mark.parametrize("key", ["pga", "pgv"])
def test_im_values(params, expected, key):
    m = DBC14(Scenario(**params))
    # PGA needs to be converted from m/sec² to g, and PGV need to be
    # converted from m/sec into cm/sec
    scale = m.GRAVITY if key == "pga" else 0.01

    np.testing.assert_allclose(
        getattr(m, key),
        expected[key] / scale,
        rtol=RTOL,
    )


@pytest.fixture
def model():
    """Instance of the DBC13 model."""
    return DBC14(Scenario(dist_jb=10, mag=6, v_s30=600, depth_hyp=10, mechanism="SS"))


def test_ln_std(model):
    # Log10 total standard deviations from the paper
    expected = [
        0.298,
        0.309,
        0.31,
        0.313,
        0.319,
        0.323,
        0.325,
        0.328,
        0.335,
        0.338,
        0.338,
        0.337,
        0.338,
        0.337,
        0.335,
        0.334,
        0.333,
        0.33,
        0.328,
        0.328,
        0.326,
        0.325,
        0.322,
        0.322,
        0.326,
        0.328,
        0.329,
        0.33,
        0.33,
        0.329,
        0.327,
        0.328,
        0.328,
        0.327,
        0.33,
        0.331,
        0.332,
        0.332,
        0.332,
        0.331,
        0.33,
        0.331,
        0.333,
        0.335,
        0.339,
        0.343,
        0.346,
        0.353,
        0.356,
        0.359,
        0.362,
        0.365,
        0.368,
        0.368,
        0.37,
        0.37,
        0.373,
        0.375,
        0.377,
        0.378,
        0.378,
        0.376,
        0.375,
        0.375,
    ]

    np.testing.assert_allclose(np.log10(np.exp(model._ln_std)), expected, rtol=RTOL)


# Vectorized scenarios
######################
MAGS = [4.0, 5.5, 7.0]
DISTS = [5.0, 30.0, 200.0]
V_S30S = [200.0, 400.0, 800.0]
DEPTHS = [0.0, 10.0, 25.0]
MECHANISMS = ["SS", "NS", "RS"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = ["pga", "pgv", "ln_pga", "ln_std_pga", "ln_std_pgv", "spec_accels", "ln_stds"]


def scalar_results(mag, dist_jb, v_s30, depth_hyp, mechanism):
    m = DBC14(
        Scenario(
            mag=mag,
            dist_jb=dist_jb,
            v_s30=v_s30,
            depth_hyp=depth_hyp,
            mechanism=mechanism,
        )
    )
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


@pytest.fixture(scope="module")
def grid():
    rows = list(itertools.product(MAGS, DISTS, V_S30S, DEPTHS, MECHANISMS))
    mag, dist_jb, v_s30, depth_hyp, mechanism = (np.array(c) for c in zip(*rows))
    scenario = Scenario(
        mag=mag, dist_jb=dist_jb, v_s30=v_s30, depth_hyp=depth_hyp, mechanism=mechanism
    )
    expected = [scalar_results(*row) for row in rows]
    return rows, scenario, expected


@pytest.mark.parametrize("key", KEYS + ["interp_spec_accels", "interp_ln_stds"])
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    m = DBC14(scenario)
    if key.startswith("interp"):
        actual = getattr(m, key)(PERIODS)
    else:
        actual = getattr(m, key)
    desired = np.array([e[key] for e in expected])
    assert actual.shape == desired.shape
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12)
    else:
        # Each scenario uses the same matrix-vector products as a scalar scenario
        np.testing.assert_array_equal(actual, desired)


def test_vectorized_shapes(grid):
    rows, scenario, _ = grid
    m = DBC14(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 62)
    assert m.ln_stds.shape == (n, 62)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged(model):
    assert isinstance(model.pga, float)
    assert isinstance(model.ln_pga, float)
    assert isinstance(model.ln_std_pga, float)
    assert model.spec_accels.shape == (62,)
    assert model.ln_stds.shape == (62,)
    assert model.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_scalars_broadcast_with_arrays():
    mag = np.array([4.5, 5.5, 6.5])
    dist_jb = np.array([10.0, 20.0, 80.0])
    m = DBC14(
        Scenario(mag=mag, dist_jb=dist_jb, v_s30=400.0, depth_hyp=10.0, mechanism="NS")
    )
    expected = [
        scalar_results(mg, d, 400.0, 10.0, "NS")["pga"] for mg, d in zip(mag, dist_jb)
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_multidimensional_inputs():
    mag, v_s30 = np.meshgrid(
        [4.5, 5.5, 6.5], [200.0, 400.0, 600.0, 800.0], indexing="ij"
    )
    m = DBC14(
        Scenario(mag=mag, dist_jb=30.0, v_s30=v_s30, depth_hyp=10.0, mechanism="RS")
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 62)
    assert m.pga[2, 1] == scalar_results(6.5, 30.0, 400.0, 10.0, "RS")["pga"]
    np.testing.assert_array_equal(
        m.spec_accels[1, 3],
        scalar_results(5.5, 30.0, 800.0, 10.0, "RS")["spec_accels"],
    )


@pytest.mark.parametrize("mechanism", ["U", np.array(["SS", "U", "RS"])])
def test_invalid_mechanism_raises(mechanism):
    # The mechanism has no default, so an unsupported mechanism raises an error
    # after a single warning
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        with pytest.raises(ValueError, match="mechanism must be one of"):
            DBC14(
                Scenario(
                    mag=6.0,
                    dist_jb=np.full(np.shape(mechanism), 20.0),
                    v_s30=400.0,
                    depth_hyp=10.0,
                    mechanism=mechanism,
                )
            )
    messages = [str(w.message) for w in caught if "is not one of" in str(w.message)]
    messages += [str(w.message) for w in caught if "mechanism has" in str(w.message)]
    assert len(messages) == 1


def vector_scenario(n=400, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(4.0, 7.0, n),
        dist_jb=rng.uniform(5.0, 200.0, n),
        v_s30=rng.uniform(200.0, 800.0, n),
        depth_hyp=rng.uniform(0.0, 25.0, n),
        mechanism=rng.choice(MECHANISMS, n),
    )


@pytest.mark.parametrize(
    "ims,keys",
    [
        (["pga"], ["pga", "ln_pga", "ln_std_pga"]),
        ("pga", ["pga", "ln_pga", "ln_std_pga"]),
        (["pgv"], ["pgv", "ln_std_pgv"]),
        (["psa_all"], ["spec_accels", "ln_stds"]),
        (["pga", "pgv", "psa_all"], KEYS),
    ],
)
def test_ims_match_default(ims, keys):
    s = vector_scenario()
    full = DBC14(s)
    m = DBC14(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)


def test_pga_only_computes_one_period():
    m = DBC14(vector_scenario(), ims=["pga"])
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


def test_pga_only_scalar_matches_default(model):
    m = DBC14(model.scenario, ims=["pga"])
    assert m.pga == model.pga
    assert m.ln_std_pga == model.ln_std_pga
    assert isinstance(m.pga, float)


def test_single_period():
    s = vector_scenario()
    full = DBC14(s)
    m = DBC14(s, ims=["psa_1p000"])
    np.testing.assert_array_equal(m.periods, [1.0])
    assert m.spec_accels.shape == (400, 1)
    np.testing.assert_array_equal(
        m.spec_accels[:, 0], full.spec_accels[:, full.periods == 1.0][:, 0]
    )
    with pytest.raises(ValueError, match="include 'pga'"):
        m.pga


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = DBC14(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (400,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
