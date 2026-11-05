#!/usr/bin/python
import gzip
import itertools
import json
import os
import warnings

import numpy as np
import pytest

from pygmm import AkkarSandikkayaBommer2014 as ASB14
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 1e-2

# Load the tests
fname = os.path.join(os.path.dirname(__file__), "data", "asb14_tests.json.gz")
with gzip.open(fname, "rt") as fp:
    tests = json.load(fp)

testdata = [
    (dist, t["params"], results) for t in tests for dist, results in t["results"]
]


def create_model(params, dist):
    params = dict(params)
    params[dist] = params.pop("dist")
    s = Scenario(**params)
    m = ASB14(s)
    return m


@pytest.mark.parametrize("dist,params,expected", testdata)
def test_spec_accels(dist, params, expected):
    m = create_model(params, dist)
    np.testing.assert_allclose(
        m.interp_spec_accels(expected["periods"]),
        expected["spec_accels"],
        rtol=RTOL,
        err_msg="Spectral accelerations",
    )


@pytest.mark.parametrize("dist,params,expected", testdata)
@pytest.mark.parametrize("key", ["pga", "pgv"])
def test_im_values(dist, params, expected, key):
    m = create_model(params, dist)
    np.testing.assert_allclose(
        getattr(m, key),
        expected[key],
        rtol=RTOL,
    )


# Vectorized scenarios
######################
# Magnitudes on both sides of c_1 (6.75 or 7.0), v_s30 below and above V_REF
# (750 m/s) and v_con (1000 m/s), and each mechanism
MAGS = [4.5, 6.0, 6.75, 6.76, 7.0, 7.5]
DISTS = [0.0, 10.0, 75.0, 200.0]
V_S30S = [150.0, 400.0, 750.0, 750.1, 900.0, 1000.0, 1150.0]
MECHANISMS = ["SS", "NS", "RS"]
DIST_KEYS = ["dist_jb", "dist_hyp", "dist_epi"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = ["pga", "pgv", "ln_pga", "ln_std_pga", "ln_std_pgv", "spec_accels", "ln_stds"]


def scalar_results(mag, dist, v_s30, mechanism, dist_key="dist_jb"):
    m = ASB14(Scenario(mag=mag, v_s30=v_s30, mechanism=mechanism, **{dist_key: dist}))
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


@pytest.fixture(scope="module", params=DIST_KEYS)
def grid(request):
    dist_key = request.param
    rows = list(itertools.product(MAGS, DISTS, V_S30S, MECHANISMS))
    mag, dist, v_s30, mechanism = (np.array(c) for c in zip(*rows))
    scenario = Scenario(mag=mag, v_s30=v_s30, mechanism=mechanism, **{dist_key: dist})
    expected = [scalar_results(*row, dist_key) for row in rows]
    return rows, scenario, expected


@pytest.mark.parametrize("key", KEYS + ["interp_spec_accels", "interp_ln_stds"])
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    m = ASB14(scenario)
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
        np.testing.assert_array_equal(actual, desired)


def test_vectorized_shapes(grid):
    rows, scenario, _ = grid
    m = ASB14(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 62)
    assert m.ln_stds.shape == (n, 62)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged():
    m = ASB14(Scenario(mag=6.5, dist_jb=20.0, v_s30=760.0, mechanism="SS"))
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (62,)
    assert m.ln_stds.shape == (62,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_scalars_broadcast_with_arrays():
    mag = np.array([5.5, 6.5, 7.5])
    dist_hyp = np.array([5.0, 20.0, 80.0])
    m = ASB14(Scenario(mag=mag, dist_hyp=dist_hyp, v_s30=400.0, mechanism="NS"))
    expected = [
        scalar_results(mg, d, 400.0, "NS", "dist_hyp")["pga"]
        for mg, d in zip(mag, dist_hyp)
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_multidimensional_inputs():
    mag, v_s30 = np.meshgrid(
        [5.5, 6.5, 7.5], [200.0, 500.0, 760.0, 1100.0], indexing="ij"
    )
    m = ASB14(Scenario(mag=mag, dist_epi=30.0, v_s30=v_s30, mechanism="RS"))
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 62)
    assert m.pga[2, 1] == scalar_results(7.5, 30.0, 500.0, "RS", "dist_epi")["pga"]
    np.testing.assert_array_equal(
        m.spec_accels[1, 3],
        scalar_results(6.5, 30.0, 1100.0, "RS", "dist_epi")["spec_accels"],
    )


def test_invalid_mechanism_entries_use_default():
    # As for a scalar scenario, an unsupported mechanism is replaced by the
    # default (None), which has no mechanism term, like strike-slip
    mechanism = np.array(["SS", "NS", "RS", "U"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = ASB14(
            Scenario(
                mag=6.5, dist_jb=np.full(4, 20.0), v_s30=400.0, mechanism=mechanism
            )
        )
    messages = [str(w.message) for w in caught if "mechanism has" in str(w.message)]
    assert len(messages) == 1 and "1 of 4 values" in messages[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = [scalar_results(6.5, 20.0, 400.0, mech)["pga"] for mech in mechanism]
    np.testing.assert_array_equal(m.pga, expected)
    assert m.pga[3] == m.pga[0]


def test_missing_distance_raises():
    with pytest.raises(NotImplementedError):
        ASB14(Scenario(mag=np.array([6.0, 7.0]), v_s30=400.0, mechanism="SS"))


def vector_scenario(n=400, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(4.0, 8.0, n),
        dist_jb=rng.uniform(0.0, 200.0, n),
        v_s30=rng.uniform(150.0, 1200.0, n),
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
    full = ASB14(s)
    m = ASB14(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)


def test_pga_only_computes_one_period():
    m = ASB14(vector_scenario(), ims=["pga"])
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(mag=6.5, dist_jb=20.0, v_s30=400.0, mechanism="RS")
    full = ASB14(s)
    m = ASB14(s, ims=["pga"])
    assert m.pga == full.pga
    assert m.ln_std_pga == full.ln_std_pga
    assert isinstance(m.pga, float)


def test_single_period():
    s = vector_scenario()
    full = ASB14(s)
    m = ASB14(s, ims=["psa_1p000"])
    np.testing.assert_array_equal(m.periods, [1.0])
    assert m.psa_ims == ["psa_1p000"]
    assert m.spec_accels.shape == (400, 1)
    np.testing.assert_array_equal(
        m.spec_accels[:, 0], full.spec_accels[:, full.periods == 1.0][:, 0]
    )
    with pytest.raises(ValueError, match="include 'pga'"):
        m.pga


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = ASB14(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (400,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)


def test_does_not_provide_ngawest2_21():
    with pytest.raises(ValueError, match="does not provide 'psa_ngawest2_21'"):
        ASB14(vector_scenario(), ims=["psa_ngawest2_21"])
