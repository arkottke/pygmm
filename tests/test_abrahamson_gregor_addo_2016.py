import itertools
import warnings

import numpy as np
import pytest

from pygmm import AbrahamsonGregorAddo2016 as AGA16
from pygmm.model import Scenario

from . import load_tests

# Relative tolerance for all tests
RTOL = 2e-2

TESTS = load_tests("abrahamson_gregor_addo_2016.json.gz")


def create_model(params):
    s = Scenario(**params)
    m = AGA16(s)
    return m


@pytest.mark.parametrize("test", TESTS)
@pytest.mark.parametrize("key", ["spec_accels", "ln_stds", "pga", "ln_std_pga"])
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
# Distances around the forearc/backarc minimum distances (85 and 100 km)
DISTS_RUP = [10.0, 150.0]
DISTS_HYP = [20.0, 90.0]
# Depths around the maximum depth (120 km)
DEPTHS_HYP = [30.0, 150.0]
# Site conditions around v_lin (865 m/s at short periods) and 1000 m/s
V_S30S = [200.0, 800.0, 1000.0, 1200.0]
EVENT_TYPES = ["interface", "intraslab"]
TECTONIC_REGIONS = ["forearc", "backarc", "unknown"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = ["pga", "ln_pga", "ln_std_pga", "spec_accels", "ln_stds", "adjust_c1"]
NAMES = [
    "mag",
    "dist_rup",
    "dist_hyp",
    "depth_hyp",
    "v_s30",
    "event_type",
    "tectonic_region",
]


def scalar_results(*args, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = AGA16(Scenario(**dict(zip(NAMES, args))), **kwds)
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def grid_scenario(rows):
    return Scenario(**{n: np.array(c) for n, c in zip(NAMES, zip(*rows))})


@pytest.fixture(scope="module")
def grid():
    rows = list(
        itertools.product(
            MAGS,
            DISTS_RUP,
            DISTS_HYP,
            DEPTHS_HYP,
            V_S30S,
            EVENT_TYPES,
            TECTONIC_REGIONS,
        )
    )
    expected = [scalar_results(*row) for row in rows]
    return rows, grid_scenario(rows), expected


def assert_matches(actual, desired, key):
    assert actual.shape == desired.shape, key
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12, err_msg=key)
    else:
        np.testing.assert_array_equal(actual, desired, err_msg=key)


@pytest.mark.parametrize("key", KEYS + ["interp_spec_accels", "interp_ln_stds"])
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    m = AGA16(scenario)
    actual = getattr(m, key)
    if key.startswith("interp"):
        actual = actual(PERIODS)
    assert_matches(actual, np.array([e[key] for e in expected]), key)


@pytest.mark.parametrize(
    "kwds",
    [
        dict(adjust_c1=0.2),
        dict(adjust_c1=np.linspace(-0.2, 0.2, 23)),
        dict(scale_atten=0.5),
        dict(adjust_c1=-0.1, adjust_c4=1.0, scale_atten=1.5),
    ],
)
def test_options_match_scalar(kwds):
    rows = list(
        itertools.product(
            MAGS, [80.0], [90.0], [60.0], [400.0], EVENT_TYPES, ["forearc", "backarc"]
        )
    )
    m = AGA16(grid_scenario(rows), **kwds)
    expected = [scalar_results(*row, **kwds) for row in rows]
    for key in ["pga", "spec_accels", "ln_stds", "interp_spec_accels"]:
        actual = getattr(m, key)
        if key.startswith("interp"):
            actual = actual(PERIODS)
        assert_matches(actual, np.array([e[key] for e in expected]), key)


def test_options_positional():
    s = Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type="interface")
    m = AGA16(s, 0.1, 0, 0.5, ims=["pga"])
    assert m.pga == AGA16(s, 0.1, 0, 0.5).pga
    np.testing.assert_array_equal(m.adjust_c1, np.full(23, 0.1))
    assert m.scale_atten == 0.5


def test_vectorized_shapes(grid):
    rows, scenario, _ = grid
    m = AGA16(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 22)
    assert m.ln_stds.shape == (n, 22)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged():
    m = AGA16(
        Scenario(
            mag=7.0,
            dist_rup=50.0,
            dist_hyp=60.0,
            depth_hyp=40.0,
            v_s30=400.0,
            event_type="intraslab",
            tectonic_region="backarc",
        )
    )
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (22,)
    assert m.ln_stds.shape == (22,)
    assert m.adjust_c1.shape == (23,)


def test_scalars_broadcast_with_arrays():
    mag = np.array([6.0, 7.0, 8.0])
    dist_rup = np.array([20.0, 80.0, 200.0])
    m = AGA16(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            v_s30=400.0,
            event_type="interface",
            tectonic_region="backarc",
        )
    )
    expected = [
        scalar_results(mg, d, None, None, 400.0, "interface", "backarc")["pga"]
        for mg, d in zip(mag, dist_rup)
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_multidimensional_inputs():
    mag, dist_hyp = np.meshgrid([6.0, 7.0, 8.0], [50.0, 100.0, 200.0, 300.0])
    m = AGA16(
        Scenario(
            mag=mag,
            dist_hyp=dist_hyp,
            depth_hyp=60.0,
            v_s30=400.0,
            event_type="intraslab",
        )
    )
    assert m.pga.shape == (4, 3)
    assert m.spec_accels.shape == (4, 3, 22)
    assert m.adjust_c1.shape == (23,)
    expected = scalar_results(7.0, None, 200.0, 60.0, 400.0, "intraslab")
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.spec_accels[2, 1], expected["spec_accels"])


def test_unused_values_can_be_none():
    # Interface events do not need dist_hyp or depth_hyp
    mag = np.array([6.0, 7.0, 8.0])
    tectonic_region = ["forearc", "backarc", "unknown"]
    m = AGA16(
        Scenario(
            mag=mag,
            dist_rup=50.0,
            v_s30=400.0,
            event_type="interface",
            tectonic_region=np.array(tectonic_region),
        )
    )
    expected = [
        scalar_results(mg, 50.0, None, None, 400.0, "interface", t)["spec_accels"]
        for mg, t in zip(mag, tectonic_region)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)

    # Intraslab events do not need dist_rup
    depth_hyp = [40.0, 80.0, 130.0]
    m = AGA16(
        Scenario(
            mag=mag,
            dist_hyp=80.0,
            depth_hyp=np.array(depth_hyp),
            v_s30=400.0,
            event_type="intraslab",
        )
    )
    expected = [
        scalar_results(mg, None, 80.0, d, 400.0, "intraslab")["spec_accels"]
        for mg, d in zip(mag, depth_hyp)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


@pytest.mark.parametrize(
    "missing,match",
    [
        ("dist_hyp", "dist_hyp is required"),
        ("depth_hyp", "depth_hyp is required"),
        ("dist_rup", "dist_rup is required"),
    ],
)
def test_missing_values_raise(missing, match):
    params = dict(
        mag=7.0,
        dist_rup=50.0,
        dist_hyp=60.0,
        depth_hyp=40.0,
        v_s30=400.0,
        event_type=np.array(["interface", "intraslab"]),
    )
    del params[missing]
    with pytest.raises(ValueError, match=match):
        AGA16(Scenario(**params))


def test_invalid_entries_use_default():
    # As for a scalar scenario, an unsupported tectonic region is replaced by the
    # default ("unknown")
    tectonic_region = np.array(["forearc", "backarc", "middle", "unknown"])
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = AGA16(
            Scenario(
                mag=7.0,
                dist_rup=50.0,
                v_s30=400.0,
                event_type="interface",
                tectonic_region=tectonic_region,
            )
        )
    messages = [str(w.message) for w in caught if "not one of" in str(w.message)]
    assert len(messages) == 1
    assert "tectonic_region has 1 of 4 values" in messages[0]
    expected = [
        scalar_results(7.0, 50.0, None, None, 400.0, "interface", t)["pga"]
        for t in tectonic_region
    ]
    np.testing.assert_array_equal(m.pga, expected)


def test_invalid_event_types():
    # The event type has no default. As for a scalar scenario, the event is not
    # treated as an intraslab event, and the forearc/backarc scaling is not
    # defined (NaN) for a backarc site.
    params = dict(mag=7.0, dist_rup=50.0, dist_hyp=60.0, depth_hyp=40.0, v_s30=400.0)
    with pytest.warns(UserWarning, match="2 of 3 values"):
        m = AGA16(
            Scenario(
                event_type=np.array(["interface", "deep", "deep"]),
                tectonic_region=np.array(["forearc", "forearc", "backarc"]),
                **params,
            )
        )
    expected = scalar_results(7.0, 50.0, 60.0, 40.0, 400.0, "deep", "forearc")
    assert m.pga[1] == expected["pga"]
    np.testing.assert_array_equal(m.spec_accels[1], expected["spec_accels"])
    assert np.isnan(m.pga[2])
    with pytest.raises(NotImplementedError):
        scalar_results(7.0, 50.0, 60.0, 40.0, 400.0, "deep", "backarc")


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    dist_rup = rng.uniform(10.0, 300.0, n)
    return Scenario(
        mag=rng.uniform(5.0, 8.5, n),
        dist_rup=dist_rup,
        dist_hyp=dist_rup + rng.uniform(0.0, 50.0, n),
        depth_hyp=rng.uniform(20.0, 150.0, n),
        v_s30=rng.uniform(150.0, 1500.0, n),
        event_type=rng.choice(EVENT_TYPES, n),
        tectonic_region=rng.choice(TECTONIC_REGIONS, n),
    )


@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(ims):
    s = vector_scenario()
    full = AGA16(s)
    m = AGA16(s, ims=ims)
    for key in ["pga", "ln_pga", "ln_std_pga"]:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)
    # Only one period is computed
    assert m._ln_resp.shape == (500, 1)
    assert m._ln_std.shape == (500, 1)
    assert m.adjust_c1.shape == (500, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(
        mag=7.0,
        dist_rup=50.0,
        v_s30=400.0,
        event_type="interface",
        tectonic_region="backarc",
    )
    m = AGA16(s, ims=["pga"])
    assert m.pga == AGA16(s).pga
    assert isinstance(m.pga, float)


@pytest.mark.parametrize("ims", [["psa_all"], ["pga", "psa_1p000", "psa_0p200"]])
def test_psa_matches_default(ims):
    s = vector_scenario()
    full = AGA16(s)
    m = AGA16(s, ims=ims)
    cols = np.isin(full.periods, m.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])


def test_ngawest2_periods_are_not_provided():
    with pytest.raises(ValueError, match="does not provide 'psa_ngawest2_21'"):
        AGA16(vector_scenario(), ims=["psa_ngawest2_21"])


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = AGA16(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)


def test_adjust_c4_is_added_to_c4(monkeypatch):
    s = Scenario(
        mag=np.array([6.0, 7.5, 8.5]),
        dist_rup=np.array([20.0, 60.0, 150.0]),
        dist_hyp=np.array([60.0, 100.0, 200.0]),
        depth_hyp=50.0,
        v_s30=np.array([200.0, 400.0, 1000.0]),
        event_type=np.array(["interface", "intraslab", "interface"]),
        tectonic_region="forearc",
    )
    default = AGA16(s)
    adjusted = AGA16(s, adjust_c4=5.0)
    assert np.all(adjusted.spec_accels != default.spec_accels)

    # Same as increasing C_4 in the coefficients
    coeff = AGA16.COEFF.copy()
    coeff.c_4 = coeff.c_4 + 5.0
    monkeypatch.setattr(AGA16, "COEFF", coeff)
    np.testing.assert_allclose(adjusted._ln_resp, AGA16(s)._ln_resp, rtol=1e-12)
