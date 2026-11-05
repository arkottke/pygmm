#!/usr/bin/python
import itertools
import warnings

import numpy as np
import pytest
from scipy.interpolate import NearestNDInterpolator

import pygmm.ground_motion.hermkes_kuehn_riggelsen_2014 as hkr
from pygmm.ground_motion.hermkes_kuehn_riggelsen_2014 import (
    HermkesKuehnRiggelsen2014 as HKR13,
)
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 1e-2

names = ("mag", "depth_hyp", "mechanism", "dist_jb", "v_s30")

# These test cases are extracted from the HKC13 data using the
# create_hkc14_data.py script
events = [
    [4.0, 5.0, "RS", 0.0, 100.0],  # Line 2 - first
    [8.0, 5.0, "SS", 29.0, 100.0],  # Line 9472
    [7.3, 15.0, "SS", 16.0, 900.0],  # Line 1245082
    [8.0, 30.0, "NS", 200.0, 760.0],  # Line 1780057 - last
]

# Predictions. The predictions are provided in pairs of predicted variable
# and variance for PGV, PGA, PSA(T=0.1s), PSA(T=0.5s), PSA(T=1.0s), and PSA(
# T=4.0s).
predictions = [
    [
        -4.07083934,
        0.67692656,
        -1.07231844,
        0.69332252,
        -0.33641352,
        0.7994422,
        -0.97409228,
        0.74814538,
        -1.82468142,
        0.92448446,
        -5.12087126,
        1.02068773,
    ],
    [
        -0.52099869,
        0.72063588,
        1.33454867,
        0.72938761,
        1.74372217,
        0.83268984,
        2.00780052,
        0.79246556,
        1.86559946,
        0.96531185,
        0.91670043,
        1.03366518,
    ],
    [
        -1.80883591,
        0.55683512,
        0.3117025,
        0.57986321,
        0.96330332,
        0.67941943,
        0.67846516,
        0.61702095,
        0.13766521,
        0.71312718,
        -0.60224743,
        0.90132807,
    ],
    [
        -3.50927882,
        0.96678294,
        -1.96656956,
        0.91876949,
        -1.75490823,
        0.97815117,
        -1.34961607,
        1.08108927,
        -0.96706897,
        1.30808437,
        -1.82250014,
        1.10876257,
    ],
]


@pytest.mark.slow
@pytest.mark.parametrize("event,prediction", zip(events, predictions))
@pytest.mark.parametrize(
    "indices,attr",
    [
        (0, "pgv"),
        (1, "ln_std_pgv"),
        (2, "pga"),
        (3, "ln_std_pga"),
        ([2, 4, 6, 8, 10], "spec_accels"),
        ([3, 5, 7, 9, 11], "ln_stds"),
    ],
)
def test_model(event, prediction, indices, attr):
    m = HKR13(Scenario(**dict(zip(names, event))))

    prediction = np.asarray(prediction)
    if "ln_" in attr:
        # Convert from variance to standard deviation
        prediction = np.sqrt(prediction[indices])
    else:
        # Convert from log to natural space.
        prediction = np.exp(prediction[indices])

    np.testing.assert_allclose(
        getattr(m, attr),
        prediction,
        rtol=RTOL,
    )


# Vectorized scenarios
######################
# The vectorized calculation is tested with a small synthetic interpolator, which
# has the same layout as the model data: events of (mag, depth_hyp, flag_rs,
# flag_ss, flag_ns, dist_jb, v_s30), and predictions of pairs of mean and
# variance for each period. This avoids loading the large model data.
MAGS = [4.0, 5.5, 7.0, 8.0]
DEPTHS = [5.0, 15.0, 30.0]
DISTS = [0.0, 20.0, 200.0]
V_S30S = [100.0, 400.0, 1200.0]
MECHANISMS = ["SS", "NS", "RS"]
PERIODS = [0.05, 0.3, 2.0]
KEYS = ["pga", "pgv", "ln_pga", "ln_std_pga", "ln_std_pgv", "spec_accels", "ln_stds"]


def synthetic_interpolator(n=5000, seed=0):
    rng = np.random.default_rng(seed)
    # Mechanism flags of reverse, strike-slip, normal, and none
    flags = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1], [0, 0, 0]])
    events = np.column_stack(
        [
            rng.uniform(4.0, 8.0, n),
            rng.uniform(0.0, 40.0, n),
            flags[rng.integers(0, 4, n)],
            rng.uniform(0.0, 200.0, n),
            rng.uniform(100.0, 1200.0, n),
        ]
    )
    predictions = np.empty((n, 12))
    predictions[:, 0::2] = rng.normal(-2.0, 1.0, (n, 6))
    predictions[:, 1::2] = rng.uniform(0.3, 1.2, (n, 6))
    return NearestNDInterpolator(events, predictions)


SYNTHETIC_INTERPOLATOR = synthetic_interpolator()


@pytest.fixture(autouse=True)
def synthetic(request, monkeypatch):
    # The slow tests use the model data
    if "slow" not in request.keywords:
        monkeypatch.setattr(hkr, "INTERPOLATOR", SYNTHETIC_INTERPOLATOR)


def scalar_results(mag, depth_hyp, dist_jb, v_s30, mechanism):
    m = HKR13(
        Scenario(
            mag=mag,
            depth_hyp=depth_hyp,
            dist_jb=dist_jb,
            v_s30=v_s30,
            mechanism=mechanism,
        )
    )
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


@pytest.fixture
def grid():
    rows = list(itertools.product(MAGS, DEPTHS, DISTS, V_S30S, MECHANISMS))
    mag, depth_hyp, dist_jb, v_s30, mechanism = (np.array(c) for c in zip(*rows))
    scenario = Scenario(
        mag=mag, depth_hyp=depth_hyp, dist_jb=dist_jb, v_s30=v_s30, mechanism=mechanism
    )
    expected = [scalar_results(*row) for row in rows]
    return rows, scenario, expected


@pytest.mark.parametrize("key", KEYS + ["interp_spec_accels", "interp_ln_stds"])
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    m = HKR13(scenario)
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
    m = HKR13(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 5)
    assert m.ln_stds.shape == (n, 5)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged():
    m = HKR13(Scenario(mag=6.5, dist_jb=20.0, v_s30=400.0, mechanism="SS"))
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (5,)
    assert m.ln_stds.shape == (5,)


def test_scalars_broadcast_with_arrays():
    # depth_hyp uses the default value
    mag = np.array([4.5, 5.5, 6.5])
    dist_jb = np.array([10.0, 20.0, 80.0])
    m = HKR13(Scenario(mag=mag, dist_jb=dist_jb, v_s30=400.0, mechanism="NS"))
    expected = [
        scalar_results(mg, 15.0, d, 400.0, "NS")["pga"] for mg, d in zip(mag, dist_jb)
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_multidimensional_inputs():
    mag, v_s30 = np.meshgrid(
        [4.5, 5.5, 6.5], [200.0, 400.0, 600.0, 800.0], indexing="ij"
    )
    m = HKR13(Scenario(mag=mag, dist_jb=30.0, v_s30=v_s30, mechanism="RS"))
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 5)
    assert m.pga[2, 1] == scalar_results(6.5, 15.0, 30.0, 400.0, "RS")["pga"]
    np.testing.assert_array_equal(
        m.ln_stds[1, 3], scalar_results(5.5, 15.0, 30.0, 800.0, "RS")["ln_stds"]
    )


def test_invalid_mechanism_entries_use_default():
    # As for a scalar scenario, an unsupported mechanism is replaced by the
    # default (None), which has all of the mechanism flags equal to zero
    mechanism = np.array(["SS", "NS", "RS", "U"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = HKR13(
            Scenario(
                mag=6.5, dist_jb=np.full(4, 20.0), v_s30=400.0, mechanism=mechanism
            )
        )
    messages = [str(w.message) for w in caught if "mechanism has" in str(w.message)]
    assert len(messages) == 1 and "1 of 4 values" in messages[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = [
            scalar_results(6.5, 15.0, 20.0, 400.0, mech)["pga"] for mech in mechanism
        ]
    np.testing.assert_array_equal(m.pga, expected)


def vector_scenario(n=400, seed=1):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(4.0, 8.0, n),
        depth_hyp=rng.uniform(0.0, 40.0, n),
        dist_jb=rng.uniform(0.0, 200.0, n),
        v_s30=rng.uniform(100.0, 1200.0, n),
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
    full = HKR13(s)
    m = HKR13(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)


def test_pga_only_computes_one_period():
    m = HKR13(vector_scenario(), ims=["pga"])
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(mag=6.5, dist_jb=20.0, v_s30=400.0, mechanism="RS")
    full = HKR13(s)
    m = HKR13(s, ims=["pga"])
    assert m.pga == full.pga
    assert m.ln_std_pga == full.ln_std_pga
    assert isinstance(m.pga, float)


def test_single_period():
    s = vector_scenario()
    full = HKR13(s)
    m = HKR13(s, ims=["psa_1p000"])
    np.testing.assert_array_equal(m.periods, [1.0])
    assert m.spec_accels.shape == (400, 1)
    np.testing.assert_array_equal(m.spec_accels[:, 0], full.spec_accels[:, 3])
    with pytest.raises(ValueError, match="include 'pga'"):
        m.pga


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = HKR13(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (400,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
