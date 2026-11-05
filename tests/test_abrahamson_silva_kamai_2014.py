#!/usr/bin/env python
"""Test calculation of ASK14 static methods and the vectorized model."""

import itertools
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pygmm
from pygmm import AbrahamsonSilvaKamai2014 as ASK14


def test_depth_1_0():
    # Value calculated from NGAW2 spreadsheet
    assert_allclose(ASK14.calc_depth_1_0(600), 0.1424470, rtol=1e-5)


def test_depth_tor():
    # Value calculated from NGAW2 spreadsheet
    assert_allclose(ASK14.calc_depth_tor(6), 4.2545455)


def test_width():
    # Value calculated from NGAW2 spreadsheet
    assert_allclose(ASK14.calc_width(6, 50), 8.9125094)


MAGS = [3.5, 5.0, 5.5, 6.0, 6.75, 7.0, 8.0]
# Sites on the hanging wall, close and far from the surface trace of the
# rupture: (dist_rup, dist_jb, dist_x)
DISTS = [(3.0, 0.0, 2.0), (12.0, 5.0, 15.0), (60.0, 55.0, 80.0)]
V_S30S = [180.0, 400.0, 760.0, 1180.0]
MECHANISMS = ["SS", "NS", "RS"]
REGIONS = ["global", "china", "japan", "taiwan"]
ON_HANGING_WALL = [True, False]
PERIODS = [0.05, 0.3, 1.0, 2.5]
KEYS = [
    "pga",
    "pgv",
    "ln_pga",
    "ln_std_pga",
    "ln_std_pgv",
    "spec_accels",
    "ln_stds",
]


def scalar_results(**kwds):
    kwds.setdefault("dip", 50.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(pygmm.Scenario(**kwds))
    results = {key: getattr(m, key) for key in KEYS}
    results["tau"] = m._tau
    results["phi"] = m._phi
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def vector_results(scenario, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(scenario, **kwds)
    results = {key: getattr(m, key) for key in KEYS}
    results["tau"] = m._tau
    results["phi"] = m._phi
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def assert_results_match(actual, expected):
    for key in actual:
        desired = np.array([e[key] for e in expected])
        assert actual[key].shape == desired.shape, key
        if key.startswith("interp"):
            # interp1d evaluates 2-D arrays with slightly different floating
            # point round-off
            np.testing.assert_allclose(actual[key], desired, rtol=1e-12, err_msg=key)
        else:
            np.testing.assert_array_equal(actual[key], desired, err_msg=key)


@pytest.fixture(scope="module")
def grid():
    rows = [
        dict(
            mag=mag,
            dist_rup=dist[0],
            dist_jb=dist[1],
            dist_x=dist[2],
            v_s30=v_s30,
            mechanism=mechanism,
            region=region,
            on_hanging_wall=on_hanging_wall,
        )
        for mag, dist, v_s30, mechanism, region, on_hanging_wall in itertools.product(
            MAGS, DISTS, V_S30S, MECHANISMS, REGIONS, ON_HANGING_WALL
        )
    ]
    scenario = pygmm.Scenario(
        dip=50.0, **{k: np.array([r[k] for r in rows]) for k in rows[0]}
    )
    expected = [scalar_results(**row) for row in rows]
    return rows, scenario, expected


@pytest.mark.parametrize(
    "key", KEYS + ["tau", "phi", "interp_spec_accels", "interp_ln_stds"]
)
def test_vectorized_matches_scalar(grid, key):
    rows, scenario, expected = grid
    actual = vector_results(scenario)
    assert_results_match({key: actual[key]}, expected)


def random_rows(n, seed):
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n):
        rows.append(
            dict(
                mag=rng.uniform(3.5, 8.0),
                dist_rup=rng.uniform(0.0, 200.0),
                dist_jb=rng.uniform(0.0, 200.0),
                dist_x=rng.uniform(-20.0, 60.0),
                dist_y0=rng.uniform(0.0, 20.0),
                dist_crjb=rng.uniform(0.0, 20.0),
                v_s30=rng.choice([180.0, 300.0, 600.0, 1180.0, rng.uniform(180, 1000)]),
                depth_1_0=rng.uniform(0.0, 1.5),
                depth_tor=rng.uniform(0.0, 15.0),
                width=rng.uniform(2.0, 30.0),
                dip=rng.uniform(20.0, 90.0),
                mechanism=rng.choice(MECHANISMS),
                region=rng.choice(
                    ["global", "california", "china", "italy", "japan", "taiwan"]
                ),
                vs_source=rng.choice(["measured", "inferred"]),
                is_aftershock=rng.choice([True, False]),
                on_hanging_wall=rng.choice([True, False]),
            )
        )
    return rows


def as_scenario(rows):
    return pygmm.Scenario(**{k: np.array([r[k] for r in rows]) for k in rows[0]})


def test_all_scenario_values_match_scalar():
    rows = random_rows(300, 1)
    expected = [scalar_results(**row) for row in rows]
    assert_results_match(vector_results(as_scenario(rows)), expected)


@pytest.mark.parametrize(
    "missing",
    [
        ["width"],
        ["depth_tor"],
        ["dist_y0"],
        ["depth_1_0"],
        ["width", "depth_tor", "dist_y0", "depth_1_0", "dist_crjb"],
    ],
)
def test_estimated_values_match_scalar(missing):
    # Values that are not provided are estimated, or their terms are skipped
    rows = [
        {k: v for k, v in row.items() if k not in missing}
        for row in random_rows(200, 2)
    ]
    expected = [scalar_results(**row) for row in rows]
    assert_results_match(vector_results(as_scenario(rows)), expected)


def test_estimated_width_and_depth_tor():
    mag = np.array([4.5, 6.0, 7.5])
    dip = np.array([30.0, 50.0, 90.0])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(
            pygmm.Scenario(
                mag=mag,
                dip=dip,
                dist_rup=20.0,
                dist_jb=20.0,
                v_s30=400.0,
                mechanism="SS",
            )
        )
    np.testing.assert_array_equal(
        m.scenario.width, [ASK14.calc_width(mg, d) for mg, d in zip(mag, dip)]
    )
    np.testing.assert_array_equal(
        m.scenario.depth_tor, [ASK14.calc_depth_tor(mg) for mg in mag]
    )


def test_footwall_without_dist_x():
    # dist_x is only needed for sites on the hanging wall
    rows = [
        {k: v for k, v in row.items() if k != "dist_x"} for row in random_rows(50, 3)
    ]
    for row in rows:
        row["on_hanging_wall"] = False
    expected = [scalar_results(**row) for row in rows]
    assert_results_match(vector_results(as_scenario(rows)), expected)


def test_hanging_wall_requires_dist_x():
    on_hanging_wall = np.array([True, False])
    with pytest.raises(ValueError, match="dist_x"):
        ASK14(
            pygmm.Scenario(
                mag=6.5,
                dist_rup=10.0,
                dist_jb=5.0,
                dip=45.0,
                v_s30=400.0,
                mechanism="RS",
                on_hanging_wall=on_hanging_wall,
            )
        )


def test_vectorized_shapes(grid):
    rows, scenario, expected = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.spec_accels.shape == (n, 22)
    assert m.ln_stds.shape == (n, 22)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def scalar_scenario(**kwds):
    return pygmm.Scenario(
        **{
            **dict(
                mag=6.5,
                dist_rup=20.0,
                dist_jb=15.0,
                dist_x=10.0,
                dip=60.0,
                v_s30=400.0,
                mechanism="RS",
                on_hanging_wall=True,
            ),
            **kwds,
        }
    )


def test_scalar_shapes_are_unchanged():
    m = ASK14(scalar_scenario())
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.pgv, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (22,)
    assert m.ln_stds.shape == (22,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_scalars_broadcast_with_arrays():
    mag = np.array([5.5, 6.5, 7.5])
    dist_rup = np.array([5.0, 20.0, 80.0])
    m = ASK14(scalar_scenario(mag=mag, dist_rup=dist_rup, region="japan"))
    expected = [
        scalar_results(**{**scalar_scenario(mag=mg, dist_rup=d, region="japan").data})[
            "spec_accels"
        ]
        for mg, d in zip(mag, dist_rup)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [5.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    m = ASK14(scalar_scenario(mag=mag, dist_rup=dist_rup))
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 22)
    expected = scalar_results(**scalar_scenario(mag=7.5, dist_rup=20.0).data)
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.spec_accels[2, 1], expected["spec_accels"])


@pytest.mark.parametrize(
    "key,values,default",
    [
        ("mechanism", ["SS", "U", "RS", "XX"], None),
        ("region", ["japan", "mars", "taiwan", "nz"], "global"),
        ("vs_source", ["measured", "guess", "inferred", "other"], "measured"),
    ],
)
def test_invalid_entries_use_default(key, values, default):
    # As for a scalar scenario, an unsupported value is replaced by the default
    values = np.array(values)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = ASK14(scalar_scenario(**{key: values}))
    messages = [str(w.message) for w in caught if key in str(w.message)]
    assert len(messages) == 2
    assert "2 of 4 values" in messages[0]
    expected = [scalar_results(**scalar_scenario(**{key: v}).data) for v in values]
    for name in ["spec_accels", "ln_stds"]:
        np.testing.assert_array_equal(
            getattr(m, name), [e[name] for e in expected], err_msg=name
        )


def vector_scenario(n=400, seed=0):
    return as_scenario(random_rows(n, seed))


@pytest.mark.parametrize(
    "ims,keys",
    [
        (["pga"], ["pga", "ln_pga", "ln_std_pga"]),
        ("pga", ["pga", "ln_pga", "ln_std_pga"]),
        (["pgv"], ["pgv", "ln_std_pgv"]),
        (["psa_all"], ["spec_accels", "ln_stds"]),
        (["pga", "pgv"], ["pga", "pgv", "ln_std_pga", "ln_std_pgv"]),
        (["psa_1p000", "pga"], ["pga", "ln_std_pga"]),
    ],
)
def test_ims_match_default(ims, keys):
    s = vector_scenario()
    full = vector_results(s)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), full[key], err_msg=key)


def test_pga_only_computes_one_period():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(vector_scenario(), ims=["pga"])
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


def test_pga_only_scalar_matches_default():
    full = ASK14(scalar_scenario())
    m = ASK14(scalar_scenario(), ims=["pga"])
    assert m.pga == full.pga
    assert m.ln_std_pga == full.ln_std_pga
    assert isinstance(m.pga, float)


def test_single_period_matches_default():
    s = vector_scenario()
    full = vector_results(s)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(s, ims=["psa_1p000"])
    assert m.psa_ims == ["psa_1p000"]
    index = list(ASK14.PERIODS[ASK14.INDICES_PSA]).index(1.0)
    np.testing.assert_array_equal(m.spec_accels[:, 0], full["spec_accels"][:, index])


def test_ngawest2_21_periods():
    s = vector_scenario()
    full = vector_results(s)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(s, ims=["psa_ngawest2_21"])
    np.testing.assert_allclose(m.periods, ASK14.PERIODS_NGAWEST2_21)
    assert m.spec_accels.shape == (400, 21)
    # ASK14 also provides 6 s, which is not one of the 21 periods
    keep = ASK14.PERIODS[ASK14.INDICES_PSA] != 6.0
    np.testing.assert_array_equal(m.spec_accels, full["spec_accels"][:, keep])
    np.testing.assert_array_equal(m.ln_stds, full["ln_stds"][:, keep])


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    s = vector_scenario()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(s, ims=ims)
    assert m.ln_pga.shape == (400,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_allclose(m.ln_pga, np.log(m.pga), rtol=1e-14)


def test_ln_pga_scalar():
    m = ASK14(scalar_scenario())
    assert isinstance(m.ln_pga, float)
    assert np.exp(m.ln_pga) == m.pga


@pytest.mark.parametrize(
    "ims,attr,im", [(["pga"], "pgv", "pgv"), (["pga"], "spec_accels", "psa_all")]
)
def test_values_not_computed_raise(ims, attr, im):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = ASK14(vector_scenario(), ims=ims)
    with pytest.raises(ValueError, match=f"include '{im}'"):
        getattr(m, attr)


def test_does_not_provide_pgd():
    with pytest.raises(ValueError, match="does not provide 'pgd'"):
        ASK14(scalar_scenario(), ims=["pgd"])


def test_calc_depth_1_0_accepts_arrays():
    v_s30 = np.array([200.0, 400.0, 760.0, 1100.0])
    region = np.array(["japan", "california", "japan", "global"])
    expected = [ASK14.calc_depth_1_0(v, r) for v, r in zip(v_s30, region)]
    np.testing.assert_array_equal(ASK14.calc_depth_1_0(v_s30, region), expected)
    np.testing.assert_array_equal(
        ASK14.calc_depth_1_0(v_s30), [ASK14.calc_depth_1_0(v) for v in v_s30]
    )


def test_calc_width_and_depth_tor_accept_arrays():
    mag = np.array([4.0, 5.5, 6.5, 7.8])
    dip = np.array([30.0, 45.0, 70.0, 90.0])
    np.testing.assert_array_equal(
        ASK14.calc_width(mag, dip), [ASK14.calc_width(m, d) for m, d in zip(mag, dip)]
    )
    np.testing.assert_array_equal(
        ASK14.calc_depth_tor(mag), [ASK14.calc_depth_tor(m) for m in mag]
    )


@pytest.mark.parametrize("depth_1_0", [None, 0.3])
def test_calc_site_term_accepts_arrays(depth_1_0):
    resp_ref = np.exp(np.linspace(-4, -1, 22 + 2))
    v_s30 = np.array([150.0, 200.0, 400.0, 900.0, 1180.0])
    site = ASK14.calc_site_term(resp_ref, v_s30[:, None], depth_1_0, "japan")
    expected = [ASK14.calc_site_term(resp_ref, v, depth_1_0, "japan") for v in v_s30]
    np.testing.assert_array_equal(site, expected)


def test_interp_coeffs_matches_interp1d():
    from scipy.interpolate import interp1d

    c = ASK14.COEFF
    xp = [150, 250, 350, 450, 600, 850, 1150]
    fps = [c.a36, c.a37, c.a38, c.a39, c.a40, c.a41, c.a42]
    x = np.array([100.0, 150.0, 200.0, 250.0, 449.9, 1000.0, 1150.0, 1500.0])
    expected = interp1d(
        xp, np.c_[tuple(fps)], bounds_error=False, fill_value=(fps[0], fps[-1])
    )(x).T
    # Equal to within rounding: the arithmetic of interp1d differs by a unit in
    # the last place between SciPy versions (e.g., 1.17 and 1.18)
    np.testing.assert_allclose(
        ASK14._interp_coeffs(xp, fps, x[:, None]), expected, rtol=1e-12, atol=1e-15
    )
