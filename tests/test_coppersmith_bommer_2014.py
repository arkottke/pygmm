import itertools
import warnings

import numpy as np
import pytest

from pygmm import CoppersmithBommer2014 as Model
from pygmm.model import Scenario

from . import load_tests

# Relative tolerance for all tests
RTOL = 2e-2

TESTS = load_tests("coppersmith_bommer_2014.json.gz")


def create_model(params):
    s = Scenario(**params)
    m = Model(s)
    return m


@pytest.mark.parametrize("test", TESTS)
@pytest.mark.parametrize("key", ["spec_accels", "ln_stds"])
def test_spec_accels(test, key):
    m = create_model(test["params"])
    np.testing.assert_allclose(
        getattr(m, key), test["results"][key], rtol=RTOL, err_msg=key
    )


# Vectorized scenarios
######################

# Magnitudes around the break points (C1 + adjust_c1): 7.5 for intraslab events,
# and 8.0 to 7.6 for interface events depending on the period
MAGS = [5.0, 7.5, 7.6, 7.8, 8.0, 8.4]
# Distances around the forearc/backarc minimum distance (40 km)
DISTS_RUP = [10.0, 40.0, 80.0, 500.0]
# Site conditions around v_lin (865 m/s at short periods) and 1000 m/s
V_S30S = [200.0, 800.0, 1000.0, 1200.0]
EVENT_TYPES = ["interface", "intraslab"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = ["pga", "ln_pga", "ln_std_pga", "spec_accels", "ln_stds", "adjust_c1"]


def scalar_results(mag, dist_rup, v_s30, event_type, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = Model(
            Scenario(mag=mag, dist_rup=dist_rup, v_s30=v_s30, event_type=event_type),
            **kwds,
        )
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def grid_scenario(rows):
    mag, dist_rup, v_s30, event_type = (np.array(c) for c in zip(*rows))
    return Scenario(mag=mag, dist_rup=dist_rup, v_s30=v_s30, event_type=event_type)


def assert_matches(actual, desired, key):
    assert actual.shape == desired.shape, key
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12, err_msg=key)
    else:
        np.testing.assert_array_equal(actual, desired, err_msg=key)


@pytest.mark.parametrize("scale_atten", [1.0, 0.5])
def test_vectorized_matches_scalar(scale_atten):
    rows = list(itertools.product(MAGS, DISTS_RUP, V_S30S, EVENT_TYPES))
    m = Model(grid_scenario(rows), scale_atten=scale_atten)
    expected = [scalar_results(*row, scale_atten=scale_atten) for row in rows]
    for key in KEYS + ["interp_spec_accels", "interp_ln_stds"]:
        actual = getattr(m, key)
        if key.startswith("interp"):
            actual = actual(PERIODS)
        assert_matches(actual, np.array([e[key] for e in expected]), key)


def test_scale_atten_positional():
    s = Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type="interface")
    m = Model(s, 0.5, ims=["pga"])
    assert m.pga == Model(s, 0.5).pga
    assert m.scale_atten == 0.5
    assert m.pga != Model(s).pga


def test_vectorized_shapes():
    rows = list(itertools.product(MAGS, DISTS_RUP, V_S30S, EVENT_TYPES))
    m = Model(grid_scenario(rows))
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 14)
    assert m.ln_stds.shape == (n, 14)


def test_scalar_shapes_are_unchanged():
    m = Model(Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type="intraslab"))
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (14,)
    assert m.ln_stds.shape == (14,)


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid([6.0, 7.0, 8.0], [20.0, 50.0, 200.0, 300.0])
    m = Model(Scenario(mag=mag, dist_rup=dist_rup, v_s30=400.0, event_type="interface"))
    assert m.pga.shape == (4, 3)
    assert m.spec_accels.shape == (4, 3, 14)
    expected = scalar_results(7.0, 200.0, 400.0, "interface")
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.spec_accels[2, 1], expected["spec_accels"])


def test_invalid_entries_use_default():
    # The event type has no default. As for a scalar scenario, the event is
    # treated as an interface event.
    event_type = np.array(["interface", "deep", "intraslab"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = Model(Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type=event_type))
    messages = [str(w.message) for w in caught if "not one of" in str(w.message)]
    assert len(messages) == 1
    assert "event_type has 1 of 3 values" in messages[0]
    expected = [scalar_results(7.0, 50.0, 400.0, e)["pga"] for e in event_type]
    np.testing.assert_array_equal(m.pga, expected)


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(5.0, 8.5, n),
        dist_rup=rng.uniform(10.0, 500.0, n),
        v_s30=rng.uniform(150.0, 1500.0, n),
        event_type=rng.choice(EVENT_TYPES, n),
    )


@pytest.mark.parametrize("ims", [["pga"], "pga"])
def test_pga_only_matches_default(ims):
    s = vector_scenario()
    full = Model(s)
    m = Model(s, ims=ims)
    for key in ["pga", "ln_pga", "ln_std_pga"]:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)
    # Only one period is computed
    assert m._ln_resp.shape == (500, 1)


@pytest.mark.parametrize("ims", [["psa_all"], ["psa_0p010", "psa_1p000"]])
def test_psa_matches_default(ims):
    s = vector_scenario()
    full = Model(s)
    m = Model(s, ims=ims)
    cols = np.isin(full.periods, m.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])
    # PGA is the response at 0.01 s, so it is available with that period
    np.testing.assert_array_equal(m.pga, full.pga)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = Model(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
