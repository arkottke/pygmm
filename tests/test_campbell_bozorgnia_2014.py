#!/usr/bin/env python
"""Test calculation of CB14 static methods and the vectorized model."""

import itertools
import logging
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pygmm
from pygmm import CampbellBozorgnia2014 as CB14


def test_depth_2_5():
    # Value calculated from NGAW2 spreadsheet
    assert_allclose(CB14.calc_depth_2_5(600, "japan", None), 0.1844427)
    assert_allclose(CB14.calc_depth_2_5(600, "california", None), 0.7952589)


@pytest.mark.parametrize("dist_rup", [50.0, 150.0])
@pytest.mark.parametrize("vectorized", [False, True])
def test_china_anelastic_attenuation(dist_rup, vectorized):
    # China uses its own anelastic attenuation coefficient (dc_20ch) beyond
    # 80 km, and the global coefficient (dc_20ca = 0) otherwise
    regions = ["global", "china"]
    kwds = dict(mag=6.5, dist_rup=dist_rup, dist_jb=dist_rup, dist_x=-dist_rup)
    # The site term is linear above k_1, so it does not depend on the region
    kwds.update(dip=90.0, v_s30=1300.0, mechanism="SS")
    if vectorized:
        m = CB14(pygmm.Scenario(region=np.array(regions), **kwds))
        ln_resp_global, ln_resp_china = m._ln_resp
    else:
        ln_resp_global, ln_resp_china = (
            CB14(pygmm.Scenario(region=region, **kwds))._ln_resp for region in regions
        )
    c = CB14.COEFF
    assert np.any(c.dc_20ch != c.dc_20ca)
    assert_allclose(
        ln_resp_china - ln_resp_global,
        (c.dc_20ch - c.dc_20ca) * max(dist_rup - 80, 0),
        atol=1e-12,
    )


# Magnitudes cover the breakpoints of the magnitude, style of faulting,
# hanging wall, hypocentral depth, dip, and standard deviation terms
MAGS = [3.5, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.8]
# Rupture, Joyner-Boore, and Rx distances: zero distance, footwall, hanging
# wall over the rupture (Rx <= R1), hanging wall beyond the rupture, and
# distances beyond 80 km with anelastic attenuation
DISTS = [
    (0.0, 0.0, 0.0),
    (5.0, 3.0, -5.0),
    (10.0, 2.0, 8.0),
    (30.0, 25.0, 30.0),
    (120.0, 118.0, -10.0),
]
DIPS = [30.0, 90.0]
# The estimated Z2.5 is greater than 3 km, between 1 and 3 km, and less than
# 1 km. Japan has a site correction below 200 m/s, and the nonlinear site
# term applies below k_1.
V_S30S = [180.0, 400.0, 900.0, 1300.0]
MECHANISMS = ["SS", "NS", "RS"]
REGIONS = ["global", "japan", "italy", "china"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = ["pga", "pgv", "ln_pga", "ln_std_pga", "ln_std_pgv", "spec_accels", "ln_stds"]


def scalar_results(**kwds):
    m = CB14(pygmm.Scenario(**kwds))
    results = {key: getattr(m, key) for key in KEYS}
    results["tau"] = m._tau
    results["phi"] = m._phi
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    return results


@pytest.fixture(scope="module")
def grid():
    rows = [
        dict(
            mag=mag,
            dist_rup=dist[0],
            dist_jb=dist[1],
            dist_x=dist[2],
            dip=dip,
            v_s30=v_s30,
            mechanism=mechanism,
            region=region,
        )
        for mag, dist, dip, v_s30, mechanism, region in itertools.product(
            MAGS, DISTS, DIPS, V_S30S, MECHANISMS, REGIONS
        )
    ]
    scenario = pygmm.Scenario(
        **{key: np.array([row[key] for row in rows]) for key in rows[0]}
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = [scalar_results(**row) for row in rows]
    return rows, scenario, expected


def assert_matches(actual, desired, key):
    desired = np.array(desired)
    assert actual.shape == desired.shape
    if key == "interp_spec_accels":
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12)
    else:
        np.testing.assert_array_equal(actual, desired)


def get_value(m, key):
    if key == "interp_spec_accels":
        return m.interp_spec_accels(PERIODS)
    elif key in ("tau", "phi"):
        return getattr(m, "_" + key)
    else:
        return getattr(m, key)


@pytest.mark.parametrize("key", KEYS + ["tau", "phi", "interp_spec_accels"])
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = CB14(scenario)
    assert_matches(get_value(m, key), [e[key] for e in expected], key)


def test_vectorized_shapes(grid):
    rows, scenario, expected = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = CB14(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 21)
    assert m.ln_stds.shape == (n, 21)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


# Provided values that are otherwise estimated: Z_TOR below and above 16.66 km,
# Z2.5 less than 1, between 1 and 3, and greater than 3 km, and the width,
# hypocentral depth, and depth to the bottom of the seismogenic crust
OPTIONAL = dict(
    depth_tor=[None, 0.0, 5.0, 17.0],
    depth_2_5=[None, 0.5, 2.0, 5.0],
    width=[None, 8.0, 25.0],
    depth_hyp=[None, 4.0, 12.0],
    depth_bot=[None, 20.0],
    depth_1_0=[None, 0.3],
)


@pytest.mark.parametrize("name", list(OPTIONAL))
@pytest.mark.parametrize("region", ["global", "japan"])
def test_optional_values_match_scalar(name, region):
    # Each provided value, with the other optional values estimated
    rows = [
        dict(
            mag=mag,
            dist_rup=20.0,
            dist_jb=10.0,
            dist_x=12.0,
            dip=dip,
            v_s30=v_s30,
            mechanism=mechanism,
            region=region,
            **{name: value},
        )
        for mag, dip, v_s30, mechanism, value in itertools.product(
            [4.0, 5.8, 7.2], [40.0, 90.0], [190.0, 760.0], MECHANISMS, OPTIONAL[name]
        )
        if value is not None
    ]
    scenario = pygmm.Scenario(
        **{key: np.array([row[key] for row in rows]) for key in rows[0]}
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = CB14(scenario)
        expected = [scalar_results(**row) for row in rows]
    for key in ["pga", "spec_accels", "ln_stds"]:
        assert_matches(get_value(m, key), [e[key] for e in expected], key)


def test_estimated_values_match_scalar(grid):
    rows, scenario, expected = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = CB14(scenario).scenario
        for key in ["depth_2_5", "depth_tor", "width", "depth_bor", "depth_hyp"]:
            desired = [
                getattr(CB14(pygmm.Scenario(**row)).scenario, key) for row in rows[::7]
            ]
            np.testing.assert_array_equal(getattr(s, key)[::7], desired, err_msg=key)


def test_scalar_shapes_are_unchanged():
    m = CB14(
        pygmm.Scenario(
            mag=6.5,
            dist_rup=20.0,
            dist_jb=15.0,
            dist_x=10.0,
            dip=45.0,
            v_s30=400.0,
            mechanism="RS",
        )
    )
    assert isinstance(m.pga, float)
    assert isinstance(m.pgv, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (21,)
    assert m.ln_stds.shape == (21,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


SCALAR_KWDS = dict(dist_jb=10.0, dist_x=5.0, dip=50.0, v_s30=350.0, mechanism="RS")


def test_scalars_broadcast_with_arrays():
    mag = np.array([5.0, 6.5, 7.5])
    dist_rup = np.array([12.0, 20.0, 90.0])
    m = CB14(pygmm.Scenario(mag=mag, dist_rup=dist_rup, **SCALAR_KWDS))
    expected = [
        scalar_results(mag=mg, dist_rup=d, **SCALAR_KWDS)["spec_accels"]
        for mg, d in zip(mag, dist_rup)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [12.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    m = CB14(pygmm.Scenario(mag=mag, dist_rup=dist_rup, **SCALAR_KWDS))
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    expected = scalar_results(mag=7.5, dist_rup=20.0, **SCALAR_KWDS)
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.ln_stds[2, 1], expected["ln_stds"])


def vector_scenario(n=400, seed=0):
    rng = np.random.default_rng(seed)
    dist_jb = rng.uniform(0.0, 200.0, n)
    return pygmm.Scenario(
        mag=rng.uniform(4.0, 7.5, n),
        dist_rup=dist_jb + rng.uniform(0.0, 10.0, n),
        dist_jb=dist_jb,
        dist_x=rng.uniform(-50.0, 50.0, n),
        dip=rng.uniform(30.0, 90.0, n),
        v_s30=rng.uniform(180.0, 1300.0, n),
        mechanism=rng.choice(["SS", "NS", "RS"], n),
        region=rng.choice(["global", "california", "japan", "italy", "china"], n),
    )


@pytest.mark.parametrize(
    "ims,keys",
    [
        (["pga"], ["pga", "ln_pga", "ln_std_pga"]),
        ("pga", ["pga", "ln_pga", "ln_std_pga"]),
        (["pgv"], ["pgv", "ln_std_pgv"]),
        (["psa_all"], ["spec_accels", "ln_stds"]),
        (["psa_ngawest2_21"], ["spec_accels", "ln_stds"]),
        (["pga", "pgv"], ["pga", "pgv", "ln_std_pga", "ln_std_pgv"]),
        (["pgv", "psa_all", "pga"], ["pga", "pgv", "spec_accels", "ln_stds"]),
    ],
)
def test_ims_match_default(ims, keys):
    s = vector_scenario()
    full = CB14(s)
    m = CB14(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)


def test_single_period_matches_default():
    s = vector_scenario()
    full = CB14(s)
    m = CB14(s, ims=["psa_1p000"])
    assert m.psa_ims == ["psa_1p000"]
    np.testing.assert_array_equal(m.spec_accels[:, 0], full.spec_accels[:, 13])
    np.testing.assert_array_equal(m.ln_stds[:, 0], full.ln_stds[:, 13])


def test_pga_only_computes_one_period():
    m = CB14(vector_scenario(), ims=["pga"])
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


def test_pga_only_scalar_matches_default():
    s = pygmm.Scenario(mag=6.5, dist_rup=20.0, **SCALAR_KWDS)
    full = CB14(s)
    m = CB14(s, ims=["pga"])
    assert m.pga == full.pga
    assert m.ln_std_pga == full.ln_std_pga
    assert isinstance(m.pga, float)


@pytest.mark.parametrize(
    "ims,attr,im",
    [
        (["pga"], "pgv", "pgv"),
        (["pgv"], "pga", "pga"),
        (["pga"], "spec_accels", "psa_all"),
    ],
)
def test_values_not_computed_raise(ims, attr, im):
    m = CB14(vector_scenario(), ims=ims)
    with pytest.raises(ValueError, match=f"include '{im}'"):
        getattr(m, attr)


def test_cb14_does_not_provide_pgd():
    with pytest.raises(ValueError, match="does not provide 'pgd'"):
        CB14(vector_scenario(), ims=["pgd"])


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    s = vector_scenario()
    m = CB14(s, ims=ims)
    assert m.ln_pga.shape == (400,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_array_equal(m.ln_pga, CB14(s).ln_pga)


def test_invalid_region_entries_use_default():
    region = np.array(["global", "mars", "japan", "atlantis"])
    kwds = dict(mag=6.5, dist_rup=20.0, **SCALAR_KWDS)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = CB14(pygmm.Scenario(region=region, **kwds))
    messages = [
        str(w.message) for w in caught if "not one of the options" in str(w.message)
    ]
    assert len(messages) == 1
    assert "2 of 4 values" in messages[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = [scalar_results(region=r, **kwds)["spec_accels"] for r in region]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_invalid_mechanism_entries_match_scalar():
    # The mechanism has no default, so an invalid mechanism gets no style of
    # faulting term, as for a scalar scenario
    mechanism = np.array(["SS", "U", "RS", "NS"])
    kwds = dict(mag=6.5, dist_rup=20.0, dist_jb=10.0, dist_x=5.0, dip=50.0, v_s30=350.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = CB14(pygmm.Scenario(mechanism=mechanism, **kwds))
        expected = [scalar_results(mechanism=mech, **kwds)["pga"] for mech in mechanism]
    np.testing.assert_array_equal(m.pga, expected)


def test_magnitude_warnings_are_logged_once(caplog):
    mag = np.array([7.6, 7.8, 8.2, 6.0])
    mechanism = np.array(["NS", "NS", "SS", "NS"])
    with caplog.at_level(logging.WARNING):
        CB14(
            pygmm.Scenario(
                mag=mag,
                mechanism=mechanism,
                dist_rup=20.0,
                dist_jb=10.0,
                dist_x=5.0,
                dip=50.0,
                v_s30=350.0,
            )
        )
    messages = [r.getMessage() for r in caplog.records]
    assert len(messages) == 1
    assert "2 of 3 earthquakes, maximum of 7.8" in messages[0]
    assert "NS style" in messages[0]


def test_helpers_accept_arrays():
    mag = np.array([4.0, 5.5, 6.8, 7.9])
    dip = np.array([20.0, 45.0, 70.0, 90.0])
    depth_tor = np.array([0.0, 2.0, 5.0, 10.0])
    v_s30 = np.array([180.0, 400.0, 760.0, 1300.0])
    region = np.array(["japan", "global", "japan", "california"])

    def loop(func, *args):
        return [func(*row) for row in zip(*args)]

    np.testing.assert_array_equal(
        CB14.calc_depth_2_5(v_s30, region), loop(CB14.calc_depth_2_5, v_s30, region)
    )
    depth_1_0 = np.array([0.05, 0.2, 0.5, 1.0])
    np.testing.assert_array_equal(
        CB14.calc_depth_2_5(None, region, depth_1_0),
        [CB14.calc_depth_2_5(None, r, z) for r, z in zip(region, depth_1_0)],
    )
    width = CB14.calc_width(mag, dip, depth_tor)
    np.testing.assert_array_equal(width, loop(CB14.calc_width, mag, dip, depth_tor))
    depth_bor = CB14.calc_depth_bor(depth_tor, dip, width)
    np.testing.assert_array_equal(
        depth_bor, loop(CB14.calc_depth_bor, depth_tor, dip, width)
    )
    np.testing.assert_array_equal(
        CB14.calc_depth_hyp(mag, dip, depth_tor, depth_bor),
        loop(CB14.calc_depth_hyp, mag, dip, depth_tor, depth_bor),
    )

    pga_ref = np.array([0.05, 0.3, 0.8, 0.1])
    depth_2_5 = np.array([0.5, 2.0, 5.0, 1.0])
    site = CB14.calc_site_term(
        pga_ref[:, None], v_s30[:, None], depth_2_5[:, None], region[:, None]
    )
    np.testing.assert_array_equal(
        site, loop(CB14.calc_site_term, pga_ref, v_s30, depth_2_5, region)
    )


def test_helpers_scalar_types():
    assert isinstance(CB14.calc_depth_2_5(600.0, "japan"), float)
    assert isinstance(CB14.calc_width(6.0, 45.0, 2.0), float)
    assert isinstance(CB14.calc_depth_hyp(6.0, 45.0, 2.0, 12.0), float)
    assert CB14.calc_site_term(0.2, 400.0, 1.5).shape == (23,)


def soft_site_scenario():
    # Strong nonlinear site response reduces short-period PSA below PGA
    mag = np.array([5.5, 6.5, 7.5, 8.0])
    return pygmm.Scenario(
        mag=mag,
        dist_rup=np.full(4, 5.0),
        dist_jb=np.full(4, 5.0),
        dist_x=np.full(4, -5.0),
        dip=90.0,
        v_s30=150.0,
        mechanism="SS",
        region="global",
    )


def test_short_period_psa_not_less_than_pga():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = CB14(soft_site_scenario())
    is_short = m.periods < 0.25
    # The rule applies for these scenarios
    assert np.any(m.spec_accels[:, is_short] == m.pga[:, np.newaxis])
    assert np.all(m.spec_accels[:, is_short] >= m.pga[:, np.newaxis])


@pytest.mark.parametrize("ims", [["psa_0p050"], ["psa_0p050", "psa_1p000"]])
def test_short_period_psa_without_pga_matches_default(ims):
    # PGA is computed for the rule even if it is not requested
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        full = CB14(soft_site_scenario())
        m = CB14(soft_site_scenario(), ims=ims)
    cols = np.searchsorted(full.periods, m.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
