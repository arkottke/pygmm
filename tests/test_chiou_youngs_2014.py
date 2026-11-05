#!/usr/bin/env python
"""Test calculation of CY14 static methods and the vectorized model."""

import itertools
import logging
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pygmm
from pygmm import ChiouYoungs2014 as CY14


def test_depth_1_0():
    # Value calculated from NGAW2 spreadsheet
    assert_allclose(CY14.calc_depth_1_0(600, "california"), 0.1259, 4)
    assert_allclose(CY14.calc_depth_1_0(600, "japan"), 0.0331, 4)


def test_depth_tor():
    # Value calculated from NGAW2 spreadsheet
    assert_allclose(CY14.calc_depth_tor(6, "SS"), 2.2587685)
    assert_allclose(CY14.calc_depth_tor(6, "RS"), 6.3447262)


# Magnitudes cover the breakpoints of the depth to top of rupture model (4.97,
# 5.849), the directivity term (5.5, 6.3), the uncertainty model (5, 6.5), and
# the regional adjustment for Japan and Italy (6 < M < 6.9)
MAGS = [3.5, 4.9, 5.6, 6.2, 6.95, 8.0]
DISTS = [0.5, 45.0, 120.0]
V_S30S = [180.0, 400.0, 1300.0]
MECHANISMS = ["U", "SS", "NS", "RS"]
REGIONS = ["california", "china", "italy", "japan"]
HANGING_WALL = [True, False]
VS_SOURCES = ["measured", "inferred"]
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


def scalar_model(**kwds):
    return CY14(pygmm.Scenario(**kwds))


def scalar_results(**kwds):
    m = scalar_model(**kwds)
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def grid_rows(estimated):
    rows = []
    for i, (mag, dist, v_s30, mechanism, region, hw, vs_source) in enumerate(
        itertools.product(
            MAGS, DISTS, V_S30S, MECHANISMS, REGIONS, HANGING_WALL, VS_SOURCES
        )
    ):
        row = dict(
            mag=mag,
            dist_rup=dist,
            dist_jb=0.8 * dist,
            dist_x=(0.5 if hw else -0.5) * dist,
            dip=[30.0, 60.0, 90.0][i % 3],
            v_s30=v_s30,
            mechanism=mechanism,
            region=region,
            on_hanging_wall=hw,
            vs_source=vs_source,
            dpp_centered=[0.0, -0.5, 1.2][i % 3],
        )
        if not estimated:
            row["depth_tor"] = [0.0, 3.0, 12.0][(i // 3) % 3]
            row["depth_1_0"] = [0.05, 0.4, 1.5][(i // 9) % 3]
        rows.append(row)
    return rows


@pytest.fixture(scope="module", params=[True, False], ids=["estimated", "given"])
def grid(request):
    # Depth to top of rupture and Z1.0 are either estimated by the model (None)
    # or specified for each scenario
    rows = grid_rows(request.param)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        logging.disable(logging.WARNING)
        try:
            expected = [scalar_results(**row) for row in rows]
            scenario = pygmm.Scenario(
                **{k: np.array([row[k] for row in rows]) for k in rows[0]}
            )
            m = CY14(scenario)
        finally:
            logging.disable(logging.NOTSET)
    return rows, m, expected


@pytest.mark.parametrize("key", KEYS + ["interp_spec_accels", "interp_ln_stds"])
def test_vectorized_matches_scalar(grid, key):
    rows, m, expected = grid
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
    rows, m, _ = grid
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 24)
    assert m.ln_stds.shape == (n, 24)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


SCENARIO = dict(
    mag=6.5,
    dist_rup=20.0,
    dist_jb=15.0,
    dist_x=10.0,
    dip=45.0,
    v_s30=360.0,
    mechanism="RS",
    on_hanging_wall=True,
)


def test_scalar_shapes_are_unchanged():
    m = scalar_model(**SCENARIO)
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (24,)
    assert m.ln_stds.shape == (24,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_scalars_broadcast_with_arrays():
    # Arrays of magnitude and distance with scalar values for the others
    mag = np.array([5.5, 6.5, 7.5])
    dist_rup = np.array([5.0, 20.0, 80.0])
    kwds = dict(SCENARIO, mag=mag, dist_rup=dist_rup, dist_jb=dist_rup)
    m = CY14(pygmm.Scenario(**kwds))
    expected = [
        scalar_results(**dict(SCENARIO, mag=mg, dist_rup=d, dist_jb=d))["spec_accels"]
        for mg, d in zip(mag, dist_rup)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_multidimensional_inputs():
    # A grid of magnitudes and distances keeps its shape
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [5.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    kwds = dict(SCENARIO, mag=mag, dist_rup=dist_rup, dist_jb=dist_rup)
    m = CY14(pygmm.Scenario(**kwds))
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 24)
    expected = scalar_results(**dict(SCENARIO, mag=7.5, dist_rup=20.0, dist_jb=20.0))
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.ln_stds[2, 1], expected["ln_stds"])


def test_lists_are_accepted():
    m = CY14(pygmm.Scenario(**dict(SCENARIO, mag=[6.0, 7.0], mechanism=["SS", "RS"])))
    expected = [
        scalar_results(**dict(SCENARIO, mag=6.0, mechanism="SS"))["pga"],
        scalar_results(**dict(SCENARIO, mag=7.0, mechanism="RS"))["pga"],
    ]
    np.testing.assert_array_equal(m.pga, expected)


@pytest.mark.parametrize(
    "name,values",
    [
        ("mechanism", ["SS", "XX", "RS", "U"]),
        ("region", ["japan", "taiwan", "china", "california"]),
        ("vs_source", ["inferred", "guess", "measured", "inferred"]),
    ],
)
def test_invalid_entries_use_default(name, values):
    # As for a scalar scenario, an invalid entry is replaced by the default
    values = np.array(values)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        m = CY14(pygmm.Scenario(**dict(SCENARIO, **{name: values})))
    messages = [str(w.message) for w in caught if name in str(w.message)]
    assert sum("1 of 4 values" in msg for msg in messages) == 1
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        expected = [
            scalar_results(**dict(SCENARIO, **{name: v}))["pga"] for v in values
        ]
    np.testing.assert_array_equal(m.pga, expected)


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        CY14(
            pygmm.Scenario(
                **dict(SCENARIO, mag=np.array([6.0, 7.0]), dist_rup=np.ones(3))
            )
        )


def test_magnitude_warnings_are_logged_once(caplog):
    mag = np.array([8.2, 8.3, 8.2, 6.0])
    mechanism = np.array(["RS", "NS", "SS", "NS"])
    with caplog.at_level(logging.WARNING):
        CY14(pygmm.Scenario(**dict(SCENARIO, mag=mag, mechanism=mechanism)))
    messages = [r.getMessage() for r in caplog.records]
    assert len(messages) == 1
    assert "2 of 3 RS or NS earthquakes, 8.2 to 8.3" in messages[0]


def test_calc_depth_tor_accepts_arrays():
    mag = np.array([4.0, 5.5, 6.5, 7.5, 8.5])
    mechanism = np.array(["SS", "RS", "NS", "RS", "U"])
    expected = [CY14.calc_depth_tor(mg, mech) for mg, mech in zip(mag, mechanism)]
    np.testing.assert_array_equal(CY14.calc_depth_tor(mag, mechanism), expected)
    # Scalar mechanism with an array of magnitudes
    np.testing.assert_array_equal(
        CY14.calc_depth_tor(mag, "RS"), [CY14.calc_depth_tor(mg, "RS") for mg in mag]
    )


def test_calc_site_term_accepts_arrays():
    resp_ref = np.array([0.05, 0.3, 0.8])
    v_s30 = np.array([200.0, 400.0, 900.0])
    depth_1_0 = np.array([0.1, 0.5, 0.02])
    region = np.array(["japan", "california", "japan"])
    site = CY14.calc_site_term(
        resp_ref[:, None], v_s30[:, None], depth_1_0[:, None], region[:, None]
    )
    expected = [
        CY14.calc_site_term(*args) for args in zip(resp_ref, v_s30, depth_1_0, region)
    ]
    np.testing.assert_array_equal(site, expected)


def vector_scenario(n=400, seed=0):
    rng = np.random.default_rng(seed)
    dist_rup = rng.uniform(0.0, 200.0, n)
    return pygmm.Scenario(
        mag=rng.uniform(4.0, 8.0, n),
        dist_rup=dist_rup,
        dist_jb=0.9 * dist_rup,
        dist_x=rng.uniform(-50.0, 50.0, n),
        dip=rng.uniform(30.0, 90.0, n),
        v_s30=rng.uniform(180.0, 1300.0, n),
        mechanism=rng.choice(["U", "SS", "NS", "RS"], n),
        region=rng.choice(["california", "china", "italy", "japan"], n),
        on_hanging_wall=rng.choice([True, False], n),
    )


@pytest.mark.parametrize(
    "ims,keys",
    [
        (["pga"], ["pga", "ln_pga", "ln_std_pga"]),
        ("pga", ["pga", "ln_pga", "ln_std_pga"]),
        (["pgv"], ["pgv", "ln_std_pgv"]),
        (["psa_all"], ["spec_accels", "ln_stds"]),
        (["pga", "pgv"], ["pga", "pgv", "ln_std_pga", "ln_std_pgv"]),
    ],
)
def test_ims_match_default(ims, keys):
    s = vector_scenario()
    full = CY14(s)
    m = CY14(s, ims=ims)
    for key in keys:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)


def test_pga_only_computes_one_period():
    m = CY14(vector_scenario(), ims=["pga"])
    assert m._ln_resp.shape == (400, 1)
    assert m._ln_std.shape == (400, 1)


def test_pga_only_scalar_matches_default():
    s = pygmm.Scenario(**SCENARIO)
    full = CY14(s)
    m = CY14(s, ims=["pga"])
    assert m.pga == full.pga
    assert m.ln_std_pga == full.ln_std_pga
    assert isinstance(m.pga, float)


def test_psa_ngawest2_21():
    s = vector_scenario()
    full = CY14(s)
    m = CY14(s, ims=["psa_ngawest2_21"])
    assert m.spec_accels.shape == (400, 21)
    np.testing.assert_array_equal(m.periods, m.PERIODS_NGAWEST2_21)
    cols = np.isin(full.periods, m.PERIODS_NGAWEST2_21)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])


@pytest.mark.parametrize(
    "ims,attr,im",
    [
        (["pga"], "pgv", "pgv"),
        (["pgv"], "pga", "pga"),
        (["pga"], "spec_accels", "psa_all"),
    ],
)
def test_values_not_computed_raise(ims, attr, im):
    m = CY14(vector_scenario(), ims=ims)
    with pytest.raises(ValueError, match=f"include '{im}'"):
        getattr(m, attr)


def test_does_not_provide_pgd():
    with pytest.raises(ValueError, match="does not provide 'pgd'"):
        CY14(vector_scenario(), ims=["pgd"])


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    s = vector_scenario()
    m = CY14(s, ims=ims)
    assert m.ln_pga.shape == (400,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_array_equal(m.ln_pga, CY14(s).ln_pga)
    np.testing.assert_allclose(m.ln_pga, np.log(m.pga), rtol=1e-14)
