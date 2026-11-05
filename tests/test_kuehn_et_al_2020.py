"""Test the Kuehn et al. (2020) model."""

import csv
import gzip
import os
import warnings

import numpy as np
import pytest

from pygmm.ground_motion.kuehn_et_al_2020 import KuehnEtAl2020 as KBCG20
from pygmm.model import Scenario

# nshmp-lib reference results
#############################

DATA_FILE = os.path.join(
    os.path.dirname(__file__), "data", "kuehn_et_al_2020-nshmp.csv.gz"
)

# nshmp-lib Gmm id suffix -> options
VARIANTS = {
    "": {},
    "AK_ADJUSTED": dict(ak_adjusted=True),
    "BASIN": dict(basin=True),
    "SEATTLE_BASIN": dict(basin=True, seattle_basin=True),
    "SEATTLE_BASIN_M9": dict(basin=True, seattle_basin=True, m9=True),
}


def gmm_options(gmm):
    """Scenario region and event type, and model options of a nshmp-lib Gmm id."""
    _, _, region, kind, *rest = gmm.split("_")
    event_type = {"INTERFACE": "interface", "SLAB": "intraslab"}[kind]
    options = dict(VARIANTS["_".join(rest)])
    if event_type == "intraslab" and "seattle_basin" in options:
        # KBCG_20_CASCADIA_SLAB_SEATTLE_BASIN does not use the USGS basin scaling
        del options["basin"]
    return region.lower(), event_type, options


def load_reference():
    with gzip.open(DATA_FILE, "rt") as fp:
        rows = list(csv.DictReader(line for line in fp if not line.startswith("#")))
    by_gmm = {}
    for row in rows:
        by_gmm.setdefault(row["gmm"], []).append(row)
    return by_gmm


def imt_name(period):
    # nshmp-lib names, e.g., SA0P01 and SA1P0
    return (
        "SA"
        + f"{period:g}".replace(".", "P")
        + ("P0" if period >= 1 and period == int(period) else "")
    )


REFERENCE = load_reference()
# nshmp-lib IMTs in the order of the pygmm results
IMTS = ["PGA"] + [imt_name(p) for p in KBCG20.PERIODS[KBCG20.INDICES_PSA]]


def test_reference_ids():
    assert len(REFERENCE) == 14
    assert sum(len(rows) for rows in REFERENCE.values()) == 5720
    assert {r["imt"] for rows in REFERENCE.values() for r in rows} == set(IMTS)


@pytest.mark.parametrize("gmm", sorted(REFERENCE))
def test_nshmp_reference(gmm):
    rows = REFERENCE[gmm]
    region, event_type, options = gmm_options(gmm)
    inputs = sorted({int(r["index"]): r for r in rows}.items())
    scenario = Scenario(
        mag=np.array([float(r["mag"]) for _, r in inputs]),
        dist_rup=np.array([float(r["dist_rup"]) for _, r in inputs]),
        depth_tor=np.array([float(r["depth_tor"]) for _, r in inputs]),
        v_s30=np.array([float(r["v_s30"]) for _, r in inputs]),
        depth_2_5=np.array([float(r["depth_2_5"]) for _, r in inputs]),
        event_type=event_type,
        region=region,
    )
    with warnings.catch_warnings():
        # Some inputs are outside the recommended limits
        warnings.simplefilter("ignore")
        m = KBCG20(scenario, **options)
    median = np.column_stack([m.pga, m.spec_accels])
    ln_std = np.column_stack([m.ln_std_pga, m.ln_stds])

    expected = np.full(median.shape + (2,), np.nan)
    order = [i for i, _ in inputs]
    for r in rows:
        expected[order.index(int(r["index"])), IMTS.index(r["imt"])] = [
            float(r["median"]),
            float(r["sigma"]),
        ]
    assert not np.any(np.isnan(expected))
    # The reference values are printed with 10 decimals
    np.testing.assert_allclose(median, expected[..., 0], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(ln_std, expected[..., 1], rtol=0, atol=5.1e-11)


def test_epistemic_uses_10_km_values():
    # nshmp-lib compares log10(dist_rup) with the table distances, so the
    # epistemic uncertainty does not depend on the distance
    params = dict(mag=7.3, depth_tor=20.0, v_s30=760.0, event_type="interface")
    dist_rup = np.array([10.0, 50.0, 300.0, 1000.0])
    epi = KBCG20(Scenario(dist_rup=dist_rup, **params))
    center = KBCG20(Scenario(dist_rup=dist_rup, **params), epistemic=False)
    offset = np.log(epi.spec_accels) - np.log(center.spec_accels)
    np.testing.assert_allclose(offset, offset[[0]] * np.ones((4, 1)), rtol=1e-12)


# Model features
################

KEYS = ["pga", "ln_pga", "ln_std_pga", "pgv", "ln_std_pgv", "spec_accels", "ln_stds"]
PERIODS = [0.05, 0.3, 1.0, 2.5]
NAMES = [
    "mag",
    "dist_rup",
    "depth_tor",
    "v_s30",
    "depth_2_5",
    "event_type",
    "region",
]
EVENT_TYPES = ["interface", "intraslab"]
REGIONS = ["global", "alaska", "cascadia", "prvi"]
OPTIONS = [
    {},
    dict(epistemic=False),
    dict(basin=True),
    dict(basin=True, seattle_basin=True, m9=True),
    dict(seattle_basin=True),
    dict(ak_adjusted=True),
]


def scalar_results(*args, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = KBCG20(Scenario(**dict(zip(NAMES, args))), **kwds)
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def assert_matches(actual, desired, key):
    assert actual.shape == desired.shape, key
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12, err_msg=key)
    else:
        np.testing.assert_array_equal(actual, desired, err_msg=key)


def grid_rows():
    rng = np.random.default_rng(1)
    n = 160
    return list(
        zip(
            rng.choice([5.0, 6.5, 7.2, 7.6, 8.0, 8.6, 9.2, 9.5], n),
            rng.choice([10.0, 40.0, 150.0, 600.0], n),
            rng.choice([5.0, 20.0, 45.0, 80.0, 120.0], n),
            rng.choice([150.0, 260.0, 450.0, 760.0, 1100.0, 1500.0], n),
            rng.choice([np.nan, 0.5, 2.0, 4.0, 7.0], n),
            rng.choice(EVENT_TYPES, n),
            rng.choice(REGIONS, n),
        )
    )


def grid_scenario(rows):
    return Scenario(**{n: np.array(c) for n, c in zip(NAMES, zip(*rows))})


@pytest.fixture(scope="module")
def grid():
    rows = grid_rows()
    return rows, grid_scenario(rows)


@pytest.mark.parametrize("options", OPTIONS)
def test_vectorized_matches_scalar(grid, options):
    rows, scenario = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = KBCG20(scenario, **options)
    expected = [scalar_results(*row, **options) for row in rows]
    for key in KEYS + ["interp_spec_accels", "interp_ln_stds"]:
        if key.startswith("interp"):
            actual = getattr(m, key)(PERIODS)
        else:
            actual = getattr(m, key)
        desired = np.array([e[key] for e in expected])
        assert_matches(actual, desired, key)
    assert np.all(np.isfinite(m.spec_accels))


def test_vectorized_shapes(grid):
    rows, scenario = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = KBCG20(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.spec_accels.shape == (n, 21)
    assert m.ln_stds.shape == (n, 21)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged():
    for event_type in EVENT_TYPES:
        for region in REGIONS:
            m = KBCG20(
                Scenario(
                    mag=7.0,
                    dist_rup=80.0,
                    depth_tor=30.0,
                    v_s30=400.0,
                    event_type=event_type,
                    region=region,
                )
            )
            assert isinstance(m.pga, float)
            assert isinstance(m.ln_pga, float)
            assert isinstance(m.ln_std_pga, float)
            assert isinstance(m.pgv, float)
            assert m.spec_accels.shape == (21,)
            assert m.ln_stds.shape == (21,)
            assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_region_defaults_to_global():
    params = dict(mag=7.0, dist_rup=80.0, depth_tor=30.0, v_s30=400.0)
    m = KBCG20(Scenario(event_type="interface", **params))
    g = KBCG20(Scenario(event_type="interface", region="global", **params))
    np.testing.assert_array_equal(m.spec_accels, g.spec_accels)


def test_scalars_broadcast_with_arrays():
    mag = np.array([6.0, 7.0, 8.0])
    dist_rup = np.array([20.0, 80.0, 200.0])
    m = KBCG20(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            depth_tor=60.0,
            v_s30=300.0,
            depth_2_5=3.0,
            event_type="intraslab",
            region="cascadia",
        ),
        basin=True,
    )
    for i in range(3):
        expected = scalar_results(
            mag[i], dist_rup[i], 60.0, 300.0, 3.0, "intraslab", "cascadia", basin=True
        )
        assert m.pga[i] == expected["pga"]
        np.testing.assert_array_equal(m.spec_accels[i], expected["spec_accels"])
        np.testing.assert_array_equal(m.ln_stds[i], expected["ln_stds"])


def test_only_basin_depth_is_an_array():
    # The standard deviation does not depend on depth_2_5, but has its shape
    m = KBCG20(
        Scenario(
            mag=8.0,
            dist_rup=100.0,
            depth_tor=15.0,
            v_s30=400.0,
            depth_2_5=np.array([np.nan, 2.0, 5.0]),
            event_type="interface",
            region="cascadia",
        )
    )
    assert m.spec_accels.shape == m.ln_stds.shape == (3, 21)
    assert m.ln_std_pga.shape == (3,)


def test_mixed_event_types_and_regions():
    rows = [
        (9.0, 100.0, 15.0, 360.0, 6.5, "interface", "cascadia"),
        (7.0, 100.0, 50.0, 360.0, np.nan, "intraslab", "alaska"),
        (6.0, 50.0, 80.0, 760.0, 1.5, "intraslab", "prvi"),
        (8.0, 200.0, 20.0, 260.0, np.nan, "interface", "global"),
    ]
    for options in OPTIONS:
        m = KBCG20(grid_scenario(rows), **options)
        for i, row in enumerate(rows):
            expected = scalar_results(*row, **options)
            np.testing.assert_array_equal(m.spec_accels[i], expected["spec_accels"])
            np.testing.assert_array_equal(m.ln_stds[i], expected["ln_stds"])
            assert m.pga[i] == expected["pga"]


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [6.0, 7.0, 8.0], [20.0, 50.0, 100.0, 300.0], indexing="ij"
    )
    m = KBCG20(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            depth_tor=60.0,
            v_s30=400.0,
            event_type="intraslab",
        )
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    assert m.ln_stds.shape == (3, 4, 21)
    expected = scalar_results(8.0, 50.0, 60.0, 400.0, None, "intraslab", "global")
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.spec_accels[2, 1], expected["spec_accels"])


@pytest.mark.parametrize(
    "event_type", ["intraslab", np.array(["interface", "intraslab"])]
)
def test_missing_depth_raises(event_type):
    with pytest.raises(ValueError, match="depth_tor is required"):
        KBCG20(Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type=event_type))


def test_m9_requires_seattle_basin():
    s = Scenario(
        mag=9.0, dist_rup=50.0, depth_tor=10.0, v_s30=400.0, event_type="interface"
    )
    with pytest.raises(ValueError, match="requires seattle_basin"):
        KBCG20(s, m9=True)


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_event_types():
    # The event type has no default, so results for invalid entries are NaN
    params = dict(mag=7.0, dist_rup=50.0, depth_tor=40.0, v_s30=400.0)
    with pytest.warns(UserWarning, match="2 of 4 values"):
        m = KBCG20(
            Scenario(
                event_type=np.array(["interface", "deep", "intraslab", "crustal"]),
                **params,
            )
        )
    assert np.all(np.isnan(m.spec_accels[[1, 3]]))
    assert np.all(np.isnan(m.ln_stds[[1, 3]]))
    assert np.all(np.isnan(m.pga[[1, 3]]))
    assert np.all(np.isfinite(m.spec_accels[[0, 2]]))
    assert m.pga[0] == scalar_results(*params.values(), None, "interface", None)["pga"]
    assert m.pga[2] == scalar_results(*params.values(), None, "intraslab", None)["pga"]
    # Same for a scalar scenario
    assert np.isnan(scalar_results(*params.values(), None, "deep", None)["pga"])


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_regions_use_global():
    params = dict(mag=7.0, dist_rup=50.0, depth_tor=40.0, v_s30=400.0)
    with pytest.warns(UserWarning, match="region has 1 of 2 values"):
        m = KBCG20(
            Scenario(
                event_type="intraslab",
                region=np.array(["japan", "alaska"]),
                **params,
            )
        )
    g = scalar_results(*params.values(), None, "intraslab", "global")
    a = scalar_results(*params.values(), None, "intraslab", "alaska")
    np.testing.assert_array_equal(m.spec_accels[0], g["spec_accels"])
    np.testing.assert_array_equal(m.spec_accels[1], a["spec_accels"])


def test_options_apply_to_their_regions():
    # Each option only changes the regions and event types of the nshmp-lib
    # variants
    rows = [
        (9.0, 60.0, 10.0, 260.0, 6.5, event_type, region)
        for event_type in EVENT_TYPES
        for region in REGIONS
    ]
    s = grid_scenario(rows)
    default = KBCG20(s).spec_accels
    changed = {
        "basin": [("interface", "cascadia"), ("intraslab", "cascadia")],
        "seattle_basin": [("interface", "cascadia"), ("intraslab", "cascadia")],
        "ak_adjusted": [("interface", "global")],
    }
    for option, expected in changed.items():
        m = KBCG20(s, **{option: True}).spec_accels
        for i, (*_, event_type, region) in enumerate(rows):
            same = np.array_equal(m[i], default[i])
            assert same == ((event_type, region) not in expected), (option, i)

    # M9 only changes Cascadia interface events at long periods
    seattle = KBCG20(s, seattle_basin=True, basin=True)
    m9 = KBCG20(s, seattle_basin=True, basin=True, m9=True)
    diff = m9.spec_accels != seattle.spec_accels
    assert np.all(diff[2, m9.periods >= 2.0])
    assert not np.any(diff[2, m9.periods < 2.0])
    assert not np.any(np.delete(diff, 2, axis=0))
    # The basin term is ln(2) instead of the Seattle basin mean residual
    long = m9.periods >= 2.0
    c = KBCG20.COEFF[KBCG20.INDICES_PSA]
    np.testing.assert_allclose(
        np.log(m9.spec_accels[2, long] / seattle.spec_accels[2, long]),
        np.log(2.0) - c.mean_residual_Seattle_basin[long],
        rtol=1e-10,
    )


def test_basin_term():
    # Without the basin depth, there is no basin term
    params = dict(mag=8.0, dist_rup=100.0, depth_tor=15.0, v_s30=400.0)
    none = KBCG20(Scenario(event_type="interface", **params))
    nan = KBCG20(Scenario(event_type="interface", depth_2_5=np.nan, **params))
    np.testing.assert_array_equal(none.spec_accels, nan.spec_accels)

    # The USGS scaling removes the basin term for shallow basins and short
    # periods
    shallow = KBCG20(
        Scenario(event_type="interface", depth_2_5=1.0, region="cascadia", **params),
        basin=True,
    )
    no_basin = KBCG20(Scenario(event_type="interface", region="cascadia", **params))
    np.testing.assert_array_equal(shallow.spec_accels, no_basin.spec_accels)
    deep = KBCG20(
        Scenario(event_type="interface", depth_2_5=5.0, region="cascadia", **params),
        basin=True,
    )
    short = deep.periods <= 0.5
    np.testing.assert_array_equal(deep.spec_accels[short], no_basin.spec_accels[short])
    assert deep.pga == no_basin.pga


def test_short_periods_not_less_than_pga():
    m = KBCG20(
        Scenario(
            mag=np.array([5.0, 7.0, 9.0]),
            dist_rup=np.array([10.0, 100.0, 500.0]),
            depth_tor=20.0,
            v_s30=np.array([150.0, 760.0, 1500.0]),
            event_type="interface",
        ),
        epistemic=False,
    )
    short = m.periods <= 0.1
    assert np.all(m.spec_accels[:, short] >= m.pga[:, np.newaxis])


def test_epistemic_branches():
    s = Scenario(
        mag=np.array([7.0, 8.5, 6.5, 5.0, 9.5]),
        dist_rup=np.array([60.0, 30.0, 100.0, 200.0, 50.0]),
        depth_tor=np.array([20.0, 15.0, 50.0, 90.0, 10.0]),
        v_s30=np.array([300.0, 760.0, 200.0, 1000.0, 400.0]),
        event_type=np.array(
            ["interface", "interface", "intraslab", "intraslab", "interface"]
        ),
        region=np.array(["global", "cascadia", "alaska", "prvi", "global"]),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        epi = KBCG20(s)
        center = KBCG20(s, epistemic=False)
    assert epi.epistemic and not center.epistemic

    ln_epi = np.column_stack([epi.ln_pga, np.log(epi.spec_accels)])
    ln_center = np.column_stack([center.ln_pga, np.log(center.spec_accels)])
    # By Jensen's inequality, the collapsed median is larger than the central median
    assert np.all(ln_epi > ln_center)
    # The standard deviation is the same for the branches
    np.testing.assert_array_equal(epi.ln_stds, center.ln_stds)


def test_epistemic_magnitude_limits():
    params = dict(dist_rup=50.0, depth_tor=10.0, v_s30=400.0, event_type="interface")
    mag = np.array([3.5, 4.0, 9.5, 9.9])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        epi = KBCG20(Scenario(mag=mag, **params))
        center = KBCG20(Scenario(mag=mag, **params), epistemic=False)
    offset = np.log(epi.spec_accels) - np.log(center.spec_accels)
    np.testing.assert_allclose(offset[0], offset[1], rtol=1e-12)
    np.testing.assert_allclose(offset[2], offset[3], rtol=1e-12)


def test_options_positional():
    s = Scenario(
        mag=7.0, dist_rup=50.0, depth_tor=10.0, v_s30=400.0, event_type="interface"
    )
    assert (
        KBCG20(s, False, False, False, False, False, ["pga"]).pga
        == KBCG20(s, epistemic=False).pga
    )


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(5.0, 9.5, n),
        dist_rup=rng.uniform(10.0, 500.0, n),
        depth_tor=rng.uniform(0.0, 150.0, n),
        v_s30=rng.uniform(150.0, 1500.0, n),
        depth_2_5=np.where(rng.uniform(size=n) < 0.3, np.nan, rng.uniform(0.0, 8.0, n)),
        event_type=rng.choice(EVENT_TYPES, n),
        region=rng.choice(REGIONS, n),
    )


@pytest.mark.parametrize("options", OPTIONS)
@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(ims, options):
    s = vector_scenario()
    full = KBCG20(s, **options)
    m = KBCG20(s, ims=ims, **options)
    for key in ["pga", "ln_pga", "ln_std_pga"]:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)
    # Only one period is computed
    assert m._ln_resp.shape == (500, 1)
    assert m._ln_std.shape == (500, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(
        mag=7.0, dist_rup=50.0, depth_tor=70.0, v_s30=300.0, event_type="intraslab"
    )
    m = KBCG20(s, ims=["pga"])
    assert m.pga == KBCG20(s).pga
    assert isinstance(m.pga, float)


@pytest.mark.parametrize(
    "ims",
    [
        ["psa_all"],
        ["psa_ngawest2_21"],
        ["pga", "psa_1p000", "psa_0p200"],
        ["psa_0p050"],
        ["pgv", "psa_10p000"],
    ],
)
def test_psa_matches_default(ims):
    s = vector_scenario()
    full = KBCG20(s, basin=True, seattle_basin=True, m9=True)
    m = KBCG20(s, basin=True, seattle_basin=True, m9=True, ims=ims)
    cols = np.isin(full.periods, m.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])


def test_psa_ngawest2_21_periods():
    m = KBCG20(vector_scenario(), ims=["psa_ngawest2_21"])
    np.testing.assert_array_equal(m.periods, KBCG20.PERIODS_NGAWEST2_21)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = KBCG20(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)


def test_pgv_matches_default():
    s = vector_scenario()
    full = KBCG20(s)
    m = KBCG20(s, ims=["pgv"])
    np.testing.assert_array_equal(m.pgv, full.pgv)
    np.testing.assert_array_equal(m.ln_std_pgv, full.ln_std_pgv)
