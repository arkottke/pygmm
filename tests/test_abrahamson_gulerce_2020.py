"""Test the Abrahamson and Gülerce (2020) model.

The model is compared with the nshmp-lib (commit 44728a7d) reference results for
``AbrahamsonGulerce_2020`` (``nga-sub-interface-results.csv`` and
``nga-sub-slab-results.csv``), which nshmp-lib checks with a tolerance of 1e-10.
"""

import csv
import gzip
import itertools
import math
import os
import warnings

import numpy as np
import pytest

from pygmm.ground_motion.abrahamson_gulerce_2020 import AbrahamsonGulerce2020 as AG20
from pygmm.model import Scenario

# The reference inputs include V_s30 values above the recommended limit
pytestmark = pytest.mark.filterwarnings("ignore:v_s30 .* recommended limit")

FPATH_REF = os.path.join(
    os.path.dirname(__file__), "data", "abrahamson_gulerce_2020-nshmp.csv.gz"
)

# nshmp-lib Gmm id -> (region, options)
GMMS = {
    "AG_20_GLOBAL_INTERFACE": ("global", {}),
    "AG_20_GLOBAL_INTERFACE_AK_ADJUSTED": ("global", {"ak_adjusted": True}),
    "AG_20_GLOBAL_SLAB": ("global", {}),
    "AG_20_CASCADIA_INTERFACE": ("cascadia", {}),
    "AG_20_CASCADIA_INTERFACE_BASIN": ("cascadia", {"basin": True}),
    "AG_20_CASCADIA_INTERFACE_ADJUSTED": ("cascadia", {"adjusted": True}),
    "AG_20_CASCADIA_INTERFACE_ADJUSTED_BASIN": (
        "cascadia",
        {"adjusted": True, "basin": True},
    ),
    "AG_20_CASCADIA_SLAB": ("cascadia", {}),
    "AG_20_CASCADIA_SLAB_BASIN": ("cascadia", {"basin": True}),
    "AG_20_CASCADIA_SLAB_ADJUSTED": ("cascadia", {"adjusted": True}),
    "AG_20_CASCADIA_SLAB_ADJUSTED_BASIN": (
        "cascadia",
        {"adjusted": True, "basin": True},
    ),
    "AG_20_ALASKA_INTERFACE": ("alaska", {}),
    "AG_20_ALASKA_INTERFACE_ADJUSTED": ("alaska", {"adjusted": True}),
    "AG_20_ALASKA_SLAB": ("alaska", {}),
    "AG_20_ALASKA_SLAB_ADJUSTED": ("alaska", {"adjusted": True}),
    "AG_20_PRVI_INTERFACE": ("prvi", {}),
    "AG_20_PRVI_SLAB": ("prvi", {}),
}

IMTS = ["PGA"] + [
    f"SA{p}".replace(".", "P")
    for p in [
        "0.01", "0.02", "0.03", "0.05", "0.075", "0.1", "0.15", "0.2", "0.25",
        "0.3", "0.4", "0.5", "0.75", "1.0", "1.5", "2.0", "3.0", "4.0", "5.0",
        "7.5", "10.0",
    ]
]  # fmt: skip


def load_reference():
    """Reference inputs and results by Gmm id."""
    with gzip.open(FPATH_REF, "rt") as fp:
        lines = [line for line in fp if not line.startswith("#")]
    ref = {}
    for row in csv.DictReader(lines):
        gmm = ref.setdefault(row["gmm"], {})
        index = int(row["index"])
        case = gmm.setdefault(
            index,
            {
                k: float(row[k])
                for k in ["mag", "dist_rup", "depth_tor", "v_s30", "depth_2_5"]
            },
        )
        case[row["imt"]] = (float(row["median"]), float(row["sigma"]))
    return ref


REFERENCE = load_reference()


def test_reference_contents():
    assert set(REFERENCE) == set(GMMS)
    n = 0
    for gmm, cases in REFERENCE.items():
        assert len(cases) == (18 if "SLAB" in gmm else 19)
        for case in cases.values():
            assert set(IMTS) <= set(case)
            n += len(IMTS)
    assert n == 6930


@pytest.mark.parametrize("gmm", sorted(GMMS))
def test_nshmp_reference(gmm):
    region, options = GMMS[gmm]
    cases = REFERENCE[gmm]
    indices = sorted(cases)

    def values(key):
        return np.array([cases[i][key] for i in indices])

    s = Scenario(
        mag=values("mag"),
        dist_rup=values("dist_rup"),
        depth_tor=values("depth_tor"),
        depth_2_5=values("depth_2_5"),
        v_s30=values("v_s30"),
        event_type="intraslab" if "SLAB" in gmm else "interface",
        region=region,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = AG20(s, **options)
    median = np.column_stack([m.pga, m.spec_accels])
    sigma = np.column_stack([m.ln_std_pga, m.ln_stds])
    expected = np.array([[cases[i][imt] for imt in IMTS] for i in indices])
    # The reference values are printed with 10 decimals
    np.testing.assert_allclose(median, expected[..., 0], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(sigma, expected[..., 1], rtol=0, atol=5.1e-11)


# Literal transcription of nshmp-lib
####################################
# The reference inputs do not include intraslab events deeper than 50 km or
# V_s30 between 1000 m/sec and V_lin (up to 1085.7 m/sec), so the model is also
# compared with a deliberately literal, scalar transcription of
# ``AbrahamsonGulerce_2020.calc``.

DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "pygmm", "data")


def load_java_coeffs(name):
    with open(os.path.join(DATA_DIR, name)) as fp:
        lines = [line for line in fp if not line.startswith("#")]
    coeffs = {}
    for row in csv.DictReader(lines):
        period = float(row.pop("period"))
        coeffs["PGA" if period == 0 else period] = {k: float(v) for k, v in row.items()}
    return coeffs


def load_java_ak():
    with open(os.path.join(DATA_DIR, "nga_sub_ak_interface_adjustment.csv")) as fp:
        lines = [line for line in fp if not line.startswith("#")]
    return {
        (row["T"] if row["T"] in ("PGA", "PGV") else float(row["T"])): float(
            row["adj_ak"]
        )
        for row in csv.DictReader(lines)
    }


JAVA_COEFFS = load_java_coeffs("abrahamson_gulerce_2020.csv")
JAVA_COEFFS_PRVI = load_java_coeffs("abrahamson_gulerce_2020-prvi.csv")
JAVA_AK = load_java_ak()


def java_coefficients(imt, zone):
    c = dict(JAVA_COEFFS[imt])
    a1r = c["a1"]
    a1_adj = 0.0
    a6r = 0.0
    a12r = 0.0
    c1sr = 7.5
    eps = 0.26
    if imt != "PGA" and imt > 3.0:
        eps += 0.15 * math.log(imt / 3.0)
    if zone == "alaska":
        a1r, a1_adj, a6r, a12r, c1sr, eps = (
            c["a31"], c["cAk"], c["a24"], c["a17"], 7.9, 0.21
        )  # fmt: skip
    elif zone == "cascadia":
        a1r, a1_adj, a6r, a12r, c1sr, eps = (
            c["a32"], c["cCasc"], c["a25"], c["a18"], 7.1, 0.27
        )  # fmt: skip
    elif zone == "prvi":
        a1r = JAVA_COEFFS_PRVI[imt]["a1"]
    c.update(c1s=c1sr, eps=eps, a1=a1r, a1Adj=a1_adj, a6=c["a6"] + a6r)
    c["a12"] = JAVA_COEFFS_PRVI[imt]["a12"] if zone == "prvi" else c["a12"] + a12r
    c["ak_adjust"] = JAVA_AK[imt]
    c["imt"] = imt
    return c


def java_delta_z25_scale(imt, z2p5):
    if not math.isnan(z2p5) and z2p5 > 1.0 and imt != "PGA" and imt > 0.5:
        z_scale = (min(max(z2p5, 1.0), 3.0) - 1.0) / 2.0
        return z_scale * 0.585 if imt == 0.75 else z_scale
    return 0.0


def java_calc_mean(
    c, slab, usgs_basin, adjust, ak_adjusted, pga_rock, Mw, rRup, zTor, vs30, z2p5
):
    lnRrupHff = math.log(rRup + 10.0 * math.exp((Mw - 6.0) * c["a9"]))
    C1 = c["c1s"] if slab else c["c1i"]
    c13m = c["a13"] * (10 - Mw) * (10 - Mw)
    c4s = c["a4"] + (c["a45"] if slab else 0.0)
    fMag = (c4s if Mw <= C1 else c["a5"]) * (Mw - C1) + c13m
    fSlab = 0.0
    if slab:
        fSlab = (
            c["a10"] + (c["a4"] + c["a45"]) * (c["c1s"] - 7.5) + c["a14"] * lnRrupHff
        )
    fDepth = 0.0
    if slab:
        dzTor = (zTor - 50.0) if zTor <= 200.0 else 150.0
        fDepth += (c["a8"] if zTor <= 50.0 else c["a11"]) * dzTor
    vsS = min(vs30, 1000.0)
    fSite = c["a12"] * math.log(vsS / c["vlin"])
    if vsS < c["vlin"]:
        fSite += -c["b"] * math.log(pga_rock + 1.88) + c["b"] * math.log(
            pga_rock + 1.88 * math.pow(vsS / c["vlin"], 1.18)
        )
    else:
        fSite += c["b"] * 1.18 * math.log(vsS / c["vlin"])
    fBasin = 0.0
    if not math.isnan(z2p5):
        lnZ2p5ref = 8.52
        if vs30 >= 570:
            lnZ2p5ref = 7.6
        elif vs30 > 200:
            lnZ2p5ref = 8.52 - 0.88 * math.log(vs30 / 200.0)
        lnzPrime = math.log((z2p5 * 1000.0 + 50.0) / (math.exp(lnZ2p5ref) + 50.0))
        if lnzPrime > 0:
            if usgs_basin:
                lnzPrime *= java_delta_z25_scale(c["imt"], z2p5)
            fBasin = c["a39"] * lnzPrime
    mu = (
        c["a1"]
        + (c["a2"] + c["a3"] * (Mw - 7.0)) * lnRrupHff
        + c["a6"] * rRup
        + fMag
        + fDepth
        + fSite
        + fSlab
        + fBasin
    )
    if adjust:
        mu += c["a1Adj"]
    if not slab and ak_adjusted:
        mu += c["ak_adjust"]
    return mu


def java_phi_lin_sq(d1, d2, rRup):
    phi = d1
    if rRup > 450.0:
        phi += d2
    elif rRup >= 150.0:
        phi += d2 * (rRup - 150.0) / 300.0
    return phi


def java_calc_sigma(c, c_pga, pga_rock, rRup, vs30):
    tau_sq = 0.47 * 0.47
    phi_sq = java_phi_lin_sq(c["d1"], c["d2"], rRup)
    if min(vs30, 1000.0) < c["vlin"]:
        phi_b_sq = phi_sq - 0.09
        phi_b_sq_pga = java_phi_lin_sq(c_pga["d1"], c_pga["d2"], rRup) - 0.09
        d_site = (
            c["b"]
            * pga_rock
            * (
                1.0 / (pga_rock + 1.88 * math.pow(vs30 / c["vlin"], 1.18))
                - 1.0 / (pga_rock + 1.88)
            )
        )
        phi_sq += (
            d_site * d_site * phi_b_sq_pga
            + 2.0 * d_site * math.sqrt(phi_b_sq_pga) * math.sqrt(phi_b_sq) * c["rhoW"]
        )
        tau_sq += d_site * d_site * tau_sq + 2.0 * d_site * tau_sq * c["rhoB"]
    return math.sqrt(tau_sq + phi_sq)


def java_calc(imt, zone, slab, options, Mw, rRup, zTor, vs30, z2p5):
    """Return (ln median, sigma) as computed by nshmp-lib (combined branches)."""
    flags = (
        slab,
        options.get("basin", False),
        options.get("adjusted", False),
        options.get("ak_adjusted", False),
    )
    c = java_coefficients(imt, zone)
    c_pga = java_coefficients("PGA", zone)
    pga_rock = 0.0
    if vs30 < c["vlin"]:
        pga_rock = math.exp(
            java_calc_mean(c_pga, *flags, 0.0, Mw, rRup, zTor, 1000.0, math.nan)
        )
    mu = java_calc_mean(c, *flags, pga_rock, Mw, rRup, zTor, vs30, z2p5)
    sigma = java_calc_sigma(c, c_pga, pga_rock, rRup, vs30)
    if not options.get("epistemic", True):
        return mu, sigma
    mus = [mu - c["eps"] * 1.645, mu, mu + c["eps"] * 1.645]
    weights = [0.185, 0.63, 0.185]
    # GroundMotions.combine
    ln_median = math.log(sum(math.exp(m) * w for m, w in zip(mus, weights)))
    sigma = math.sqrt(sum(sigma * sigma * w for w in weights))
    return ln_median, sigma


EVENT_TYPES = ["interface", "intraslab"]
REGIONS = ["global", "alaska", "cascadia", "prvi"]
NAMES = ["mag", "dist_rup", "depth_tor", "depth_2_5", "v_s30", "event_type", "region"]
# Magnitudes around the break points, distances around the phi breaks (150
# and 450 km), depths around 50 and 200 km, sites around V_lin and 1000 m/sec,
# and basin depths around the reference and the USGS basin scaling (1 to 3 km)
GRID = dict(
    mag=[5.0, 7.1, 7.5, 7.9, 8.15, 9.0],
    dist_rup=[0.0, 100.0, 300.0, 600.0],
    depth_tor=[30.0, 50.0, 120.0, 250.0],
    depth_2_5=[np.nan, 0.5, 2.0, 6.5],
    v_s30=[150.0, 400.0, 760.0, 1050.0, 1500.0],
    event_type=EVENT_TYPES,
    region=REGIONS,
)
OPTIONS = [
    {},
    {"basin": True, "adjusted": True},
    {"ak_adjusted": True, "epistemic": False},
]


def grid_rows():
    return list(itertools.product(*GRID.values()))


def grid_scenario(rows):
    return Scenario(**{n: np.array(c) for n, c in zip(NAMES, zip(*rows))})


@pytest.fixture(scope="module")
def grid():
    rows = grid_rows()
    return rows, grid_scenario(rows)


@pytest.mark.parametrize("options", OPTIONS)
def test_matches_nshmp_transcription(grid, options):
    rows, scenario = grid
    m = AG20(scenario, **options)
    imts = ["PGA"] + list(AG20.PERIODS[AG20.INDICES_PSA])
    # A subset of the rows keeps the scalar transcription fast
    subset = np.arange(0, len(rows), 7)
    expected = np.array(
        [
            [
                java_calc(
                    imt,
                    region,
                    event_type == "intraslab",
                    options,
                    mag,
                    dist,
                    ztor,
                    vs,
                    z2p5,
                )
                for imt in imts
            ]
            for mag, dist, ztor, z2p5, vs, event_type, region in [
                rows[i] for i in subset
            ]
        ]
    )
    ln_resp = np.column_stack([m.ln_pga, np.log(m.spec_accels)])[subset]
    ln_std = np.column_stack([m.ln_std_pga, m.ln_stds])[subset]
    np.testing.assert_allclose(ln_resp, expected[..., 0], rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(ln_std, expected[..., 1], rtol=1e-12)


# Vectorized scenarios
######################

KEYS = ["pga", "ln_pga", "ln_std_pga", "spec_accels", "ln_stds"]
PERIODS = [0.05, 0.3, 1.0, 2.5]


def scalar_results(*args, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = AG20(Scenario(**dict(zip(NAMES, args))), **kwds)
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


@pytest.mark.slow
@pytest.mark.parametrize("options", OPTIONS)
def test_vectorized_matches_scalar(grid, options):
    rows, scenario = grid
    m = AG20(scenario, **options)
    expected = [scalar_results(*row, **options) for row in rows]
    for key in KEYS + ["interp_spec_accels", "interp_ln_stds"]:
        if key.startswith("interp"):
            actual = getattr(m, key)(PERIODS)
        else:
            actual = getattr(m, key)
        desired = np.array([e[key] for e in expected])
        assert_matches(actual, desired, key)


def test_vectorized_matches_scalar_subset(grid):
    rows, scenario = grid
    subset = np.arange(0, len(rows), 37)
    m = AG20(grid_scenario([rows[i] for i in subset]), basin=True)
    expected = [scalar_results(*rows[i], basin=True) for i in subset]
    for key in KEYS:
        desired = np.array([e[key] for e in expected])
        assert_matches(getattr(m, key), desired, key)


def test_vectorized_shapes(grid):
    rows, scenario = grid
    m = AG20(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 21)
    assert m.ln_stds.shape == (n, 21)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


@pytest.mark.parametrize("region", REGIONS)
@pytest.mark.parametrize("event_type", EVENT_TYPES)
def test_scalar_shapes_are_unchanged(event_type, region):
    m = AG20(
        Scenario(
            mag=7.0,
            dist_rup=80.0,
            depth_tor=50.0,
            depth_2_5=3.0,
            v_s30=400.0,
            event_type=event_type,
            region=region,
        )
    )
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (21,)
    assert m.ln_stds.shape == (21,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


def test_scalars_broadcast_with_arrays():
    mag = np.array([6.0, 7.0, 8.0])
    dist_rup = np.array([20.0, 80.0, 200.0])
    m = AG20(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            depth_tor=60.0,
            v_s30=300.0,
            event_type="intraslab",
            region="cascadia",
        )
    )
    for i in range(3):
        expected = scalar_results(
            mag[i], dist_rup[i], 60.0, None, 300.0, "intraslab", "cascadia"
        )
        assert m.pga[i] == expected["pga"]
        np.testing.assert_array_equal(m.spec_accels[i], expected["spec_accels"])
        np.testing.assert_array_equal(m.ln_stds[i], expected["ln_stds"])


def test_mixed_event_types_and_regions():
    event_type = np.array(["interface", "intraslab", "intraslab", "interface"])
    region = np.array(["alaska", "prvi", "global", "cascadia"])
    mag = np.array([9.0, 7.0, 6.0, 8.0])
    depth_tor = np.array([np.nan, 50.0, 120.0, np.nan])
    depth_2_5 = np.array([np.nan, 2.0, np.nan, 5.0])
    m = AG20(
        Scenario(
            mag=mag,
            dist_rup=100.0,
            depth_tor=depth_tor,
            depth_2_5=depth_2_5,
            v_s30=360.0,
            event_type=event_type,
            region=region,
        ),
        basin=True,
        adjusted=True,
    )
    for i in range(4):
        expected = scalar_results(
            mag[i],
            100.0,
            depth_tor[i],
            depth_2_5[i],
            360.0,
            event_type[i],
            region[i],
            basin=True,
            adjusted=True,
        )
        np.testing.assert_array_equal(m.spec_accels[i], expected["spec_accels"])
        np.testing.assert_array_equal(m.ln_stds[i], expected["ln_stds"])
        assert m.pga[i] == expected["pga"]
    assert np.all(np.isfinite(m.spec_accels))


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [6.0, 7.0, 8.0], [20.0, 50.0, 100.0, 300.0], indexing="ij"
    )
    m = AG20(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            depth_tor=60.0,
            v_s30=400.0,
            event_type="intraslab",
            region=np.array(["global", "alaska", "cascadia", "prvi"]),
        )
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    assert m.ln_stds.shape == (3, 4, 21)
    expected = scalar_results(8.0, 50.0, 60.0, None, 400.0, "intraslab", "alaska")
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.ln_stds[2, 1], expected["ln_stds"])


def test_interface_depth_can_be_none():
    params = dict(
        mag=np.array([7.0, 8.5]), dist_rup=np.array([30.0, 90.0]), v_s30=500.0
    )
    m = AG20(Scenario(event_type="interface", **params))
    expected = AG20(
        Scenario(event_type="interface", depth_tor=np.array([10.0, 500.0]), **params)
    )
    np.testing.assert_array_equal(m.spec_accels, expected.spec_accels)


@pytest.mark.parametrize(
    "event_type", ["intraslab", np.array(["interface", "intraslab"])]
)
def test_missing_depth_raises(event_type):
    with pytest.raises(ValueError, match="depth_tor is required"):
        AG20(Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type=event_type))


def test_nan_depth_tor_gives_nan_for_slab():
    m = AG20(
        Scenario(
            mag=7.0,
            dist_rup=50.0,
            v_s30=400.0,
            depth_tor=np.array([np.nan, np.nan]),
            event_type=np.array(["interface", "intraslab"]),
        )
    )
    assert np.isfinite(m.pga[0]) and np.isnan(m.pga[1])


def test_depth_tor_limited_to_200_km():
    params = dict(mag=7.0, dist_rup=250.0, v_s30=400.0, event_type="intraslab")
    deep = AG20(Scenario(depth_tor=np.array([200.0, 300.0, 600.0]), **params))
    assert deep.pga[0] == deep.pga[1] == deep.pga[2]


def test_basin_term():
    # The reference depth is exp(7.6) m (about 2 km) for V_s30 >= 570 m/sec
    params = dict(mag=8.0, dist_rup=100.0, v_s30=760.0, event_type="interface")
    depth_2_5 = np.array([np.nan, 1.5, 2.5, 6.0])
    none = AG20(Scenario(**params))
    plain = AG20(Scenario(depth_2_5=depth_2_5, **params))
    usgs = AG20(Scenario(depth_2_5=depth_2_5, **params), basin=True)
    for m in [plain, usgs]:
        # No basin term for NaN depths or depths less than the reference
        np.testing.assert_array_equal(m.spec_accels[:2], [none.spec_accels] * 2)
        # Only amplifies
        assert np.all(m.spec_accels[2:] >= none.spec_accels)
        # The standard deviation does not depend on the basin
        np.testing.assert_array_equal(m.ln_stds, [none.ln_stds] * 4)
    # The basin coefficient (a39) is 0 for periods shorter than 0.25 s
    has_basin = plain.periods >= 0.25
    assert np.all(plain.spec_accels[2:, has_basin] > none.spec_accels[has_basin])
    np.testing.assert_array_equal(
        plain.spec_accels[2:, ~has_basin], [none.spec_accels[~has_basin]] * 2
    )
    # The USGS basin scaling only applies to periods longer than 0.5 s and
    # reduces the basin term for depths less than 3 km
    long = usgs.periods > 0.5
    np.testing.assert_array_equal(
        usgs.spec_accels[:, ~long], [none.spec_accels[~long]] * 4
    )
    assert np.all(usgs.spec_accels[2:, long] > none.spec_accels[long])
    assert np.all(usgs.spec_accels[2, long] < plain.spec_accels[2, long])
    # At 6 km, the scale factor is 1 except at 0.75 s (0.585)
    np.testing.assert_allclose(
        usgs.spec_accels[3, usgs.periods >= 1.0],
        plain.spec_accels[3, usgs.periods >= 1.0],
        rtol=1e-14,
    )
    i = np.flatnonzero(usgs.periods == 0.75)[0]
    np.testing.assert_allclose(
        np.log(usgs.spec_accels[3, i] / none.spec_accels[i]),
        0.585 * np.log(plain.spec_accels[3, i] / none.spec_accels[i]),
        rtol=1e-12,
    )


def test_options_without_effect():
    s = Scenario(
        mag=np.array([7.0, 8.5, 6.5, 7.5]),
        dist_rup=np.array([60.0, 30.0, 100.0, 200.0]),
        depth_tor=np.array([np.nan, np.nan, 50.0, 90.0]),
        v_s30=np.array([300.0, 760.0, 200.0, 1000.0]),
        event_type=np.array(["interface", "interface", "intraslab", "intraslab"]),
    )
    for region in ["global", "prvi"]:
        # Adjusted only applies to Alaska and Cascadia
        r = s.copy_with(region=region)
        np.testing.assert_array_equal(
            AG20(r, adjusted=True).spec_accels, AG20(r).spec_accels
        )
    # The Alaska adjustment only applies to interface events
    ak = AG20(s, ak_adjusted=True)
    base = AG20(s)
    np.testing.assert_array_equal(ak.spec_accels[2:], base.spec_accels[2:])
    assert np.all(ak.spec_accels[:2] != base.spec_accels[:2])
    # For V_s30 > V_lin (linear site response), it shifts the mean. Otherwise,
    # it also changes the reference PGA of the nonlinear site term.
    linear = s.copy_with(v_s30=1100.0)
    np.testing.assert_allclose(
        np.log(AG20(linear, ak_adjusted=True).spec_accels[:2])
        - np.log(AG20(linear).spec_accels[:2]),
        np.broadcast_to(AG20.ADJ_AK[1:], (2, 21)),
        rtol=1e-12,
        atol=1e-14,
    )
    for region in ["alaska", "cascadia"]:
        r = s.copy_with(region=region)
        assert np.all(AG20(r, adjusted=True).pga != AG20(r).pga)


def test_epistemic_branches():
    s = Scenario(
        mag=np.array([7.0, 8.5, 6.5, 7.5]),
        dist_rup=np.array([60.0, 30.0, 100.0, 200.0]),
        depth_tor=np.array([np.nan, np.nan, 50.0, 90.0]),
        v_s30=np.array([300.0, 760.0, 200.0, 1000.0]),
        event_type=np.array(["interface", "interface", "intraslab", "intraslab"]),
        region=np.array(["global", "alaska", "cascadia", "prvi"]),
    )
    epi = AG20(s)
    center = AG20(s, epistemic=False)
    assert epi.epistemic and not center.epistemic
    ln_epi = np.column_stack([epi.ln_pga, np.log(epi.spec_accels)])
    ln_center = np.column_stack([center.ln_pga, np.log(center.spec_accels)])

    def offset(eps):
        eps = 1.645 * np.asarray(eps)
        return np.log(0.185 * np.exp(-eps) + 0.63 + 0.185 * np.exp(eps))

    eps_global = 0.26 + 0.15 * np.log(np.maximum(AG20.PERIODS, 3.0) / 3.0)
    expected = np.vstack(
        [
            offset(eps_global),
            offset(0.21) * np.ones(22),
            offset(0.27) * np.ones(22),
            offset(eps_global),
        ]
    )
    np.testing.assert_allclose(ln_epi - ln_center, expected, rtol=1e-12, atol=1e-14)
    np.testing.assert_array_equal(epi.ln_stds, center.ln_stds)


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_event_types():
    # The event type has no default, so results for invalid entries are NaN
    with pytest.warns(UserWarning, match="2 of 4 values"):
        m = AG20(
            Scenario(
                mag=7.0,
                dist_rup=50.0,
                depth_tor=60.0,
                v_s30=400.0,
                event_type=np.array(["interface", "deep", "intraslab", "crustal"]),
            )
        )
    assert np.all(np.isnan(m.spec_accels[[1, 3]]))
    assert np.all(np.isnan(m.pga[[1, 3]]))
    assert np.all(np.isfinite(m.spec_accels[[0, 2]]))
    assert m.pga[0] == scalar_results(7.0, 50.0, 60.0, None, 400.0, "interface")["pga"]
    assert m.pga[2] == scalar_results(7.0, 50.0, 60.0, None, 400.0, "intraslab")["pga"]
    assert np.isnan(scalar_results(7.0, 50.0, 60.0, None, 400.0, "deep")["pga"])


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_regions_use_global():
    with pytest.warns(UserWarning, match="1 of 3 values"):
        m = AG20(
            Scenario(
                mag=8.0,
                dist_rup=50.0,
                v_s30=400.0,
                event_type="interface",
                region=np.array(["cascadia", "japan", "global"]),
            )
        )
    assert m.pga[1] == m.pga[2] != m.pga[0]
    with pytest.warns(UserWarning, match="not one of the options"):
        scalar = AG20(
            Scenario(
                mag=8.0, dist_rup=50.0, v_s30=400.0, event_type="interface", region="x"
            )
        )
    assert scalar.pga == m.pga[2]


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(5.0, 9.5, n),
        dist_rup=rng.uniform(0.0, 500.0, n),
        depth_tor=rng.uniform(20.0, 250.0, n),
        depth_2_5=rng.uniform(0.0, 8.0, n),
        v_s30=rng.uniform(150.0, 1000.0, n),
        event_type=rng.choice(EVENT_TYPES, n),
        region=rng.choice(REGIONS, n),
    )


@pytest.mark.parametrize("options", OPTIONS)
@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(ims, options):
    s = vector_scenario()
    full = AG20(s, **options)
    m = AG20(s, ims=ims, **options)
    for key in ["pga", "ln_pga", "ln_std_pga"]:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)
    # Only one period is computed
    assert m._ln_resp.shape == (500, 1)
    assert m._ln_std.shape == (500, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(
        mag=7.0,
        dist_rup=50.0,
        depth_tor=70.0,
        v_s30=300.0,
        event_type="intraslab",
        region="alaska",
    )
    m = AG20(s, ims=["pga"])
    assert m.pga == AG20(s).pga
    assert isinstance(m.pga, float)


@pytest.mark.parametrize(
    "ims", [["psa_all"], ["psa_ngawest2_21"], ["pga", "psa_1p000", "psa_0p200"]]
)
def test_psa_matches_default(ims):
    s = vector_scenario()
    full = AG20(s, basin=True)
    m = AG20(s, basin=True, ims=ims)
    cols = np.isin(full.periods, m.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])


def test_psa_ngawest2_21_is_all_periods():
    m = AG20(vector_scenario(), ims=["psa_ngawest2_21"])
    np.testing.assert_array_equal(m.periods, AG20.PERIODS_NGAWEST2_21)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = AG20(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)


def test_pgv_is_not_provided():
    with pytest.raises(ValueError, match="does not provide 'pgv'"):
        AG20(vector_scenario(), ims=["pgv"])


def test_coefficients():
    c = AG20.COEFF
    assert c.period[AG20.INDEX_PGA] == 0
    np.testing.assert_array_equal(AG20.COEFF_PRVI.period, c.period)
    np.testing.assert_array_equal(
        AG20.PERIODS[AG20.INDICES_PSA], AG20.PERIODS_NGAWEST2_21
    )
    # The 0.01 s coefficients are used for PGA (the PRVI a12 values differ)
    pga = c[AG20.INDEX_PGA]
    psa_0p01 = c[1]
    for name in c.dtype.names:
        if name != "period":
            assert pga[name] == psa_0p01[name], name
    assert AG20.ADJ_AK[0] == -0.26167 and AG20.ADJ_AK[-1] == -0.22860
