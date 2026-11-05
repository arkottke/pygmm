"""Test the vectorized Boore et al. (2014) model against scalar scenarios."""

import itertools
import logging

import numpy as np
import pytest

import pygmm

BSSA14 = pygmm.BooreStewartSeyhanAtkinson2014
MAGS = [3.5, 5.0, 5.7, 6.2, 7.0, 8.0]
DISTS = [0.0, 5.0, 30.0, 200.0]
V_S30S = [180.0, 360.0, 760.0, 1300.0]
MECHANISMS = ["U", "SS", "NS", "RS"]
REGIONS = ["global", "china", "japan", "taiwan"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = ["pga", "pgv", "ln_pga", "ln_std_pga", "ln_std_pgv", "spec_accels", "ln_stds"]


def scalar_results(mag, dist_jb, v_s30, mechanism, region, depth_1_0=None):
    m = BSSA14(
        pygmm.Scenario(
            mag=mag,
            dist_jb=dist_jb,
            v_s30=v_s30,
            mechanism=mechanism,
            region=region,
            depth_1_0=depth_1_0,
        )
    )
    results = {key: getattr(m, key) for key in KEYS}
    results["tau"] = m._tau
    results["phi"] = m._phi
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    return results


@pytest.fixture(scope="module")
def grid():
    rows = list(itertools.product(MAGS, DISTS, V_S30S, MECHANISMS, REGIONS))
    mag, dist_jb, v_s30, mechanism, region = (np.array(c) for c in zip(*rows))
    scenario = pygmm.Scenario(
        mag=mag, dist_jb=dist_jb, v_s30=v_s30, mechanism=mechanism, region=region
    )
    expected = [scalar_results(*row) for row in rows]
    return rows, scenario, expected


@pytest.mark.parametrize("key", KEYS + ["tau", "phi", "interp_spec_accels"])
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    m = BSSA14(scenario)
    if key == "interp_spec_accels":
        actual = m.interp_spec_accels(PERIODS)
    elif key in ("tau", "phi"):
        actual = getattr(m, "_" + key)
    else:
        actual = getattr(m, key)
    desired = np.array([e[key] for e in expected])
    assert actual.shape == desired.shape
    if key == "interp_spec_accels":
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12)
    else:
        np.testing.assert_array_equal(actual, desired)


@pytest.mark.parametrize("depth_1_0", [0.05, 0.5, 2.0])
def test_basin_depth_matches_scalar(depth_1_0):
    # Basin depth affects periods of 0.65 s and longer
    mag = np.array([5.5, 6.5, 7.5])
    v_s30 = np.array([200.0, 400.0, 760.0])
    m = BSSA14(
        pygmm.Scenario(
            mag=mag, dist_jb=20.0, v_s30=v_s30, mechanism="SS", depth_1_0=depth_1_0
        )
    )
    expected = [
        scalar_results(mg, 20.0, v, "SS", None, depth_1_0)["spec_accels"]
        for mg, v in zip(mag, v_s30)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_basin_depth_array_matches_scalar():
    depth_1_0 = np.array([0.05, 0.5, 2.0])
    m = BSSA14(pygmm.Scenario(mag=7.0, dist_jb=20.0, v_s30=300.0, depth_1_0=depth_1_0))
    expected = [
        scalar_results(7.0, 20.0, 300.0, None, None, z)["spec_accels"]
        for z in depth_1_0
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_scalar_shapes_are_unchanged():
    m = BSSA14(pygmm.Scenario(mag=6.5, dist_jb=20.0, v_s30=760.0, mechanism="SS"))
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert m.spec_accels.shape == (105,)
    assert m.ln_stds.shape == (105,)


def vector_scenario(n=400, seed=0):
    rng = np.random.default_rng(seed)
    return pygmm.Scenario(
        mag=rng.uniform(4.0, 8.0, n),
        dist_jb=rng.uniform(0.0, 200.0, n),
        v_s30=rng.uniform(180.0, 1300.0, n),
        mechanism=rng.choice(["U", "SS", "NS", "RS"], n),
    )


@pytest.mark.parametrize(
    "ims,keys",
    [
        (["pga"], ["pga", "ln_pga", "ln_std_pga"]),
        (["pgv"], ["pgv", "ln_std_pgv"]),
        (["psa_all"], ["spec_accels", "ln_stds"]),
        (["pga", "pgv"], ["pga", "pgv", "ln_std_pga", "ln_std_pgv"]),
    ],
)
def test_ims_match_default(ims, keys):
    s = vector_scenario()
    full = BSSA14(s)
    m = BSSA14(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)


def test_pga_only_computes_one_period():
    m = BSSA14(vector_scenario(), ims="pga")
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


@pytest.mark.parametrize(
    "ims,attr,im",
    [
        (["pga"], "pgv", "pgv"),
        (["pgv"], "pga", "pga"),
        (["pga"], "spec_accels", "psa_all"),
    ],
)
def test_values_not_computed_raise(ims, attr, im):
    m = BSSA14(vector_scenario(), ims=ims)
    with pytest.raises(ValueError, match=f"include '{im}'"):
        getattr(m, attr)


def test_bssa14_does_not_provide_pgd():
    with pytest.raises(ValueError, match="does not provide 'pgd'"):
        BSSA14(vector_scenario(), ims=["pgd"])


def test_calc_site_term_accepts_arrays():
    pga_ref = np.array([0.05, 0.3, 0.8])
    v_s30 = np.array([200.0, 400.0, 900.0])
    site = BSSA14.calc_site_term(pga_ref[:, None], v_s30[:, None], None)
    expected = [BSSA14.calc_site_term(p, v, None) for p, v in zip(pga_ref, v_s30)]
    np.testing.assert_array_equal(site, expected)


def test_magnitude_warnings_are_logged_once(caplog):
    mag = np.array([7.5, 7.8, 8.0, 6.0])
    mechanism = np.array(["NS", "NS", "SS", "NS"])
    with caplog.at_level(logging.WARNING):
        BSSA14(pygmm.Scenario(mag=mag, dist_jb=10.0, v_s30=760.0, mechanism=mechanism))
    messages = [r.getMessage() for r in caplog.records]
    assert len(messages) == 1
    assert "2 of 3 normal-slip earthquakes, 7.5 to 7.8" in messages[0]


@pytest.mark.parametrize("depth_1_0", [0.05, 0.3, 1.0])
def test_japan_basin_uses_japan_z1(depth_1_0):
    # Japan and Italy share the attenuation model, so they differ only in the
    # basin term, which uses the Japan relation for Z1.0 for the Japan region
    kwds = dict(mag=6.5, dist_jb=20.0, v_s30=300.0, mechanism="SS", depth_1_0=depth_1_0)
    japan = BSSA14(pygmm.Scenario(region="japan", **kwds))
    italy = BSSA14(pygmm.Scenario(region="italy", **kwds))
    long = japan.periods >= 0.65
    np.testing.assert_array_equal(japan.spec_accels[~long], italy.spec_accels[~long])

    c = BSSA14.COEFF[BSSA14.INDICES_PSA]
    pga_ref = np.exp(
        BSSA14(pygmm.Scenario(region="italy", **kwds))._calc_ln_resp(
            np.nan, BSSA14.COEFF[[BSSA14.INDEX_PGA]]
        )[0]
    )
    basin_japan = BSSA14._calc_site_term(c, pga_ref, 300.0, depth_1_0, "japan")
    basin_global = BSSA14._calc_site_term(c, pga_ref, 300.0, depth_1_0, "california")
    np.testing.assert_allclose(
        np.log(japan.spec_accels) - np.log(italy.spec_accels),
        basin_japan - basin_global,
        atol=1e-12,
    )
    assert np.any(np.abs(basin_japan - basin_global)[long] > 1e-3)


def test_japan_without_basin_depth_matches_italy():
    kwds = dict(mag=6.5, dist_jb=20.0, v_s30=300.0, mechanism="SS")
    japan = BSSA14(pygmm.Scenario(region="japan", **kwds))
    italy = BSSA14(pygmm.Scenario(region="italy", **kwds))
    np.testing.assert_array_equal(japan.spec_accels, italy.spec_accels)


def test_region_array_with_basin_depth_matches_scalar():
    region = np.array(["japan", "global", "italy", "japan"])
    v_s30 = np.array([250.0, 300.0, 400.0, 600.0])
    m = BSSA14(
        pygmm.Scenario(mag=7.0, dist_jb=30.0, v_s30=v_s30, region=region, depth_1_0=0.5)
    )
    expected = [
        scalar_results(7.0, 30.0, v, None, r, 0.5)["spec_accels"]
        for v, r in zip(v_s30, region)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_cy14_depth_1_0_accepts_region_array():
    v_s30 = np.array([200.0, 400.0, 760.0])
    region = np.array(["japan", "california", "japan"])
    expected = [
        pygmm.ChiouYoungs2014.calc_depth_1_0(v, r) for v, r in zip(v_s30, region)
    ]
    np.testing.assert_array_equal(
        pygmm.ChiouYoungs2014.calc_depth_1_0(v_s30, region), expected
    )
