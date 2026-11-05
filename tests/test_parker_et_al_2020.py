"""Test the Parker et al. (2020) NGA-Subduction model.

The model is compared with the nshmp-lib reference results (``NgaSubInterface``
and ``NgaSubSlab`` tests) and with a deliberately literal, scalar transcription
of ``ParkerEtAl_2020`` in nshmp-lib (Java) for scenarios and options that the
reference results do not cover.
"""

import csv
import itertools
import math
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from pygmm.ground_motion.parker_et_al_2020 import ParkerEtAl2020 as PSBAH20
from pygmm.model import Scenario

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
COEFF_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "pygmm", "data")

# nshmp-lib reference results
#############################

NSHMP = pd.read_csv(
    os.path.join(DATA_DIR, "parker_et_al_2020-nshmp.csv.gz"), comment="#"
)


def gmm_options(gmm):
    """Scenario values and model options of an nshmp-lib Gmm id."""
    parts = gmm.split("_")
    assert parts[:2] == ["PSBAH", "20"]
    values = dict(
        region=parts[2].lower(),
        event_type={"INTERFACE": "interface", "SLAB": "intraslab"}[parts[3]],
    )
    options = dict(
        basin="BASIN" in parts[4:],
        m9="M9" in parts[4:],
        ak_adjusted=parts[4:] == ["AK", "ADJUSTED"],
        epistemic=parts[4:] != ["NO", "EPI"],
    )
    assert set(parts[4:]) <= {"BASIN", "M9", "AK", "ADJUSTED", "NO", "EPI"}
    return values, options


def period_columns(periods):
    return np.array(
        [np.flatnonzero(np.isclose(PSBAH20.PERIODS, p, atol=0))[0] for p in periods]
    )


def test_reference_results_coverage():
    gmms = set(NSHMP["gmm"])
    assert len(gmms) == 12
    assert len(NSHMP) == 7 * 19 * 22 + 5 * 18 * 22
    # All 22 IMTs (PGA and 21 periods) for each Gmm and input
    assert set(NSHMP["period"]) == {0.0} | set(PSBAH20.PERIODS_NGAWEST2_21)


@pytest.mark.parametrize("gmm", sorted(set(NSHMP["gmm"])))
def test_nshmp_reference_results(gmm):
    df = NSHMP[NSHMP["gmm"] == gmm]
    values, options = gmm_options(gmm)
    scenario = Scenario(
        mag=df["mag"].values,
        dist_rup=df["dist_rup"].values,
        depth_tor=df["depth_tor"].values,
        depth_2_5=df["depth_2_5"].values,
        v_s30=df["v_s30"].values,
        **values,
    )
    m = PSBAH20(scenario, **options)
    rows = np.arange(len(df))
    cols = period_columns(df["period"].values)
    median = np.exp(m._ln_resp[rows, cols])
    sigma = m._ln_std[rows, cols]
    # nshmp-lib writes the values with 10 decimals
    np.testing.assert_allclose(median, df["median"], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(sigma, df["sigma"], rtol=0, atol=5.1e-11)


def test_reference_variants_differ():
    # The reference inputs exercise the basin, M9, and Alaska adjustments
    def results(gmm):
        df = NSHMP[NSHMP["gmm"] == gmm]
        return df["median"].values

    for a, b in [
        ("PSBAH_20_CASCADIA_INTERFACE", "PSBAH_20_CASCADIA_INTERFACE_BASIN"),
        ("PSBAH_20_CASCADIA_INTERFACE_BASIN", "PSBAH_20_CASCADIA_INTERFACE_BASIN_M9"),
        ("PSBAH_20_CASCADIA_SLAB", "PSBAH_20_CASCADIA_SLAB_BASIN"),
        ("PSBAH_20_GLOBAL_INTERFACE", "PSBAH_20_GLOBAL_INTERFACE_AK_ADJUSTED"),
    ]:
        assert np.any(np.abs(results(a) - results(b)) > 1e-4), (a, b)


# Literal transcription of nshmp-lib
####################################


def load_java_coeffs(name):
    """Coefficients by period, read with the Java column names."""
    with open(os.path.join(COEFF_DIR, name)) as fp:
        lines = [line for line in fp if not line.startswith("#")]
    coeffs = {}
    for row in csv.DictReader(lines):
        coeffs[float(row.pop("period"))] = {k: float(v) for k, v in row.items()}
    return coeffs


COEFFS_INTERFACE = load_java_coeffs("parker_et_al_2020-interface.csv")
COEFFS_SLAB = load_java_coeffs("parker_et_al_2020-slab.csv")
COEFFS_PRVI = load_java_coeffs("parker_et_al_2020-prvi.csv")
with open(os.path.join(COEFF_DIR, "nga_sub_ak_interface_adjustment.csv")) as fp:
    COEFFS_AK = {
        (
            {"PGA": 0.0, "PGV": -1.0}[row["T"]]
            if row["T"] in ("PGA", "PGV")
            else float(row["T"])
        ): float(row["adj_ak"])
        for row in csv.DictReader(line for line in fp if not line.startswith("#"))
    }

PGA, PGV = 0.0, -1.0
ZONES = {"global": "Global", "alaska": "Alaska", "cascadia": "Cascadia"}


def java_erf(x):
    # Maths.erf
    def base(x):
        t = 1 / (1 + 0.3275911 * x)
        tsq = t * t
        return 1 - (
            0.254829592 * t
            + -0.284496736 * tsq
            + 1.421413741 * tsq * t
            + -1.453152027 * tsq * tsq
            + 1.061405429 * tsq * tsq * t
        ) * math.exp(-x * x)

    return -base(-x) if x < 0.0 else base(x)


def java_coefficients(imt, zone, slab):
    coeffs = (COEFFS_SLAB if slab else COEFFS_INTERFACE)[imt]
    c = dict(imt=imt)
    if zone == "prvi":
        c["c0"] = COEFFS_PRVI[imt]["c0slab" if slab else "c0"]
        c["a0"] = coeffs["Global_a0"]
        c["s2"] = COEFFS_PRVI[imt]["s2"]
    else:
        c["c0"] = coeffs[ZONES[zone] + "_c0"]
        c["a0"] = coeffs[ZONES[zone] + "_a0"]
        c["s2"] = coeffs[ZONES[zone] + "_s2"]
    c["s1"] = c["s2"]
    if zone == "cascadia":
        c["mc"] = 7.20 if slab else 7.7
    elif zone == "alaska":
        c["mc"] = 7.20 if slab else 8.6
    else:
        c["mc"] = 7.6 if slab else 7.9
    for name in "c1 c4 c5 c6 V2 f4 f5 C_e1 C_e2 C_e3 Tau phi21 phi22 phi2V".split():
        c[name] = coeffs[name]
    c["d"] = coeffs["d"] if slab else 0.0
    c["m"] = coeffs["m"] if slab else 0.0
    if zone == "cascadia":
        epi = (0.35, 0.16, 0.20, 3.00) if slab else (0.43, 0.33, 0.20, 0.50)
    elif zone == "alaska":
        epi = (0.15, 0.12, 0.50, 1.00) if slab else (0.15, 0.10, 1.00, 4.00)
    else:
        epi = (0.35, 0.22, 0.15, 2.00) if slab else (0.40, 0.40, 0.20, 0.40)
    c["se1"], c["se2"], c["t1"], c["t2"] = epi
    c["ak_adjust"] = COEFFS_AK[imt]
    return c


def java_calc_mean(c, slab, ak_adjusted, Mw, rRup, zTor, site=None):
    # Equation 4
    if slab:
        mc = c["mc"]
        h = math.pow(10, (1.050 / (mc - 4.0)) * (Mw - mc) + 1.544) if Mw <= mc else 35
    else:
        h = math.pow(10.0, -0.82 + 0.252 * Mw)
    # Equation 2 and 3
    Rref = math.sqrt(1 + h * h)
    R = math.sqrt(rRup * rRup + h * h)
    Fp = c["c1"] * math.log(R) + 0.1 * Mw * math.log(R / Rref) + c["a0"] * R
    # Equation 5
    dM = Mw - c["mc"]
    Fm = c["c4"] * dM + c["c5"] * dM * dM if Mw <= c["mc"] else c["c6"] * dM
    # Equation 6
    Fd = 0.0
    if slab:
        zHyp = zTor + 0.48 * 6.5
        if zHyp < 20.0:
            Fd = c["m"] * (20.0 - 67.0) + c["d"]
        elif zHyp > 67.0:
            Fd = c["d"]
        else:
            Fd = c["m"] * (zHyp - 67.0) + c["d"]
    mu = c["c0"] + Fp + Fm + Fd
    if site is not None:
        mu += site
    if not slab and ak_adjusted:
        mu += c["ak_adjust"]
    return mu


def java_calc_fb(c, dz):
    if dz <= c["C_e1"] / c["C_e3"]:
        return c["C_e1"]
    elif dz >= c["C_e2"] / c["C_e3"]:
        return c["C_e2"]
    return c["C_e3"] * dz


def java_site(c, pgaRef, vs30, z2p5, slab, usgsBasin, m9):
    # Equation 8
    if vs30 <= 270.0:
        Flin = c["s1"] * math.log(vs30 / 270.0) + c["s2"] * math.log(270.0 / 760.0)
    elif vs30 > c["V2"]:
        Flin = c["s2"] * math.log(c["V2"] / 760.0)
    else:
        Flin = c["s2"] * math.log(vs30 / 760.0)
    # Equation 9
    f2 = c["f4"] * (
        math.exp(c["f5"] * (min(vs30, 760) - 200)) - math.exp(c["f5"] * (760 - 200))
    )
    Fsnl = f2 * math.log((pgaRef + 0.05) / 0.05)
    # Equations 11 to 13
    Fb = 0.0
    if not math.isnan(z2p5):
        z2p5m = z2p5 * 1000.0
        lnmu = math.log(10) * (
            -0.42
            * (
                1
                + java_erf((math.log10(vs30) - math.log10(200.0)) / 0.2 / math.sqrt(2))
            )
            + 3.94
        )
        Fb = java_calc_fb(c, math.log(z2p5m) - lnmu)
        period = c["imt"]
        if m9 and not slab and z2p5 > 6.0 and period > 1.9:
            Fb += math.log(2.0) - java_calc_fb(c, math.log(z2p5m) - math.log(1279.0))
        if usgsBasin:
            # GmmUtils.deltaZ25scale
            if z2p5 > 1.0 and period > 0.5:
                scale = (min(max(z2p5, 1.0), 3.0) - 1.0) / (3.0 - 1.0)
                Fb *= scale * 0.585 if period == 0.75 else scale
            else:
                Fb *= 0.0
    return Flin + Fsnl + Fb


def java_phi_total(c, rRup, vs30):
    rPrime = max(200.0, min(500.0, rRup))
    if vs30 <= 200:
        dVar = c["phi2V"] * (math.log(500 / rPrime) / math.log(500 / 200))
    elif vs30 >= 500:
        dVar = 0
    else:
        dVar = (
            c["phi2V"]
            * (math.log(500 / vs30) / math.log(500 / 200))
            * (math.log(500 / rPrime) / math.log(500 / 200))
        )
    if rRup <= 200:
        phi2 = c["phi21"]
    elif rRup >= 500:
        phi2 = c["phi22"]
    else:
        phi2 = (c["phi22"] - c["phi21"]) / math.log(500 / 200) * math.log(
            rRup / 200
        ) + c["phi21"]
    return math.sqrt(dVar + phi2)


def java_epistemic(c):
    period = c["imt"]
    if period in (PGA, PGV) or period < c["t1"]:
        return c["se1"]
    elif period < c["t2"]:
        return c["se1"] - (c["se1"] - c["se2"]) * math.log(period / c["t1"]) / math.log(
            c["t2"] / c["t1"]
        )
    return c["se2"]


def java_calc(imt, zone, slab, basin, m9, ak, epi, Mw, rRup, zTor, vs30, z2p5):
    """Return (ln median, sigma) as computed by nshmp-lib and GroundMotions.combine."""
    # Gmm variants only exist for these combinations
    usgs_basin = basin and zone == "cascadia"
    m9 = m9 and zone == "cascadia"
    ak = ak and zone == "global"
    c = java_coefficients(imt, zone, slab)
    c_pga = java_coefficients(PGA, zone, slab)
    pgaRef = math.exp(java_calc_mean(c_pga, slab, ak, Mw, rRup, zTor))
    site = java_site(c, pgaRef, vs30, z2p5, slab, usgs_basin, m9)
    mu = java_calc_mean(c, slab, ak, Mw, rRup, zTor, site)
    sigma = math.sqrt(c["Tau"] ** 2 + java_phi_total(c, rRup, vs30) ** 2)
    if not epi:
        return mu, sigma
    eps = java_epistemic(c)
    mus = [mu - eps * 1.645, mu, mu + eps * 1.645]
    wts = [0.185, 0.63, 0.185]
    # DoubleData.weightedSumLn and Maths.srssWeighted
    mean = math.log(sum(math.exp(m) * w for m, w in zip(mus, wts)))
    sigma = math.sqrt(sum(sigma * sigma * w for w in wts))
    return mean, sigma


# Grid with all regions and event types, magnitudes around the corner
# magnitudes, depths around the depth limits (Z_hyp of 20 and 67 km), distances
# around the phi limits, V_s30 values around the site and phi limits, and
# Z_2.5 values around the basin limits (including no basin)
EVENT_TYPES = ["interface", "intraslab"]
REGIONS = ["global", "alaska", "cascadia", "prvi"]
MAGS = [5.0, 7.2, 7.5, 8.0, 8.6, 9.2]
DISTS_RUP = [10.0, 300.0, 600.0]
DEPTHS_TOR = [10.0, 40.0, 70.0]
V_S30S = [180.0, 270.0, 400.0, 1000.0, 1600.0]
DEPTHS_2_5 = [np.nan, 0.3, 2.0, 7.0]
NAMES = [
    "mag",
    "dist_rup",
    "depth_tor",
    "v_s30",
    "depth_2_5",
    "event_type",
    "region",
]


def grid_scenario(rows):
    return Scenario(**{n: np.array(c) for n, c in zip(NAMES, zip(*rows))})


@pytest.fixture(scope="module")
def grid():
    rows = list(
        itertools.product(
            MAGS, DISTS_RUP, DEPTHS_TOR, V_S30S, DEPTHS_2_5, EVENT_TYPES, REGIONS
        )
    )
    return rows, grid_scenario(rows)


OPTIONS = [
    dict(),
    dict(epistemic=False),
    dict(basin=True),
    dict(basin=True, m9=True, epistemic=False),
    dict(ak_adjusted=True),
]


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("options", OPTIONS)
def test_matches_nshmp_transcription(grid, options):
    rows, scenario = grid
    m = PSBAH20(scenario, **options)
    kwds = dict(
        basin=options.get("basin", False),
        m9=options.get("m9", False),
        ak=options.get("ak_adjusted", False),
        epi=options.get("epistemic", True),
    )
    expected = np.array(
        [
            [
                java_calc(
                    imt,
                    region,
                    event_type == "intraslab",
                    Mw=mag,
                    rRup=dist,
                    zTor=depth,
                    vs30=vs,
                    z2p5=z2p5,
                    **kwds,
                )
                for imt in PSBAH20.PERIODS
            ]
            for mag, dist, depth, vs, z2p5, event_type, region in rows
        ]
    )
    np.testing.assert_allclose(m._ln_resp, expected[..., 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(m._ln_std, expected[..., 1], rtol=0, atol=1e-12)


def test_coefficient_files():
    for c in [PSBAH20.COEFF_INTERFACE, PSBAH20.COEFF_SLAB, PSBAH20.COEFF_PRVI]:
        np.testing.assert_array_equal(c["period"], PSBAH20.PERIODS)
    assert PSBAH20.PERIODS[PSBAH20.INDEX_PGA] == 0
    assert PSBAH20.PERIODS[PSBAH20.INDEX_PGV] == -1
    np.testing.assert_allclose(
        PSBAH20.PERIODS[PSBAH20.INDICES_PSA], PSBAH20.PERIODS_NGAWEST2_21
    )
    assert PSBAH20.ADJ_AK[PSBAH20.INDEX_PGA] == -0.26167
    assert PSBAH20.ADJ_AK[PSBAH20.INDEX_PGV] == 0.00833
    assert PSBAH20.ADJ_AK[-1] == -0.22860


# Options
#########


def test_m9_requires_basin():
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0, event_type="interface")
    with pytest.raises(ValueError, match="m9=True requires basin=True"):
        PSBAH20(s, m9=True)


def test_options_only_affect_nshmp_variants():
    base = dict(mag=9.0, dist_rup=100.0, depth_tor=40.0, v_s30=300.0, depth_2_5=7.0)
    for event_type, region in itertools.product(EVENT_TYPES, REGIONS):
        s = Scenario(event_type=event_type, region=region, **base)
        ref = PSBAH20(s).spec_accels
        affects = dict(
            basin=region == "cascadia",
            ak_adjusted=region == "global" and event_type == "interface",
        )
        for option, expected in affects.items():
            changed = np.any(PSBAH20(s, **{option: True}).spec_accels != ref)
            assert changed == expected, (event_type, region, option)
        basin = PSBAH20(s, basin=True).spec_accels
        changed = np.any(PSBAH20(s, basin=True, m9=True).spec_accels != basin)
        assert changed == (region == "cascadia" and event_type == "interface")


def test_m9_long_periods_only():
    s = Scenario(
        mag=9.0,
        dist_rup=100.0,
        v_s30=300.0,
        depth_2_5=np.array([5.0, 7.0]),
        event_type="interface",
        region="cascadia",
    )
    basin = PSBAH20(s, basin=True)
    m9 = PSBAH20(s, basin=True, m9=True)
    long = basin.periods > 1.9
    # Z_2.5 <= 6 km
    np.testing.assert_array_equal(m9.spec_accels[0], basin.spec_accels[0])
    np.testing.assert_array_equal(m9.spec_accels[1, ~long], basin.spec_accels[1, ~long])
    assert np.all(m9.spec_accels[1, long] != basin.spec_accels[1, long])


def test_basin_none_and_nan():
    params = dict(mag=8.0, dist_rup=100.0, v_s30=400.0, event_type="interface")
    for region in REGIONS:
        for basin in [False, True]:
            none = PSBAH20(Scenario(region=region, **params), basin=basin)
            nan = PSBAH20(Scenario(region=region, depth_2_5=np.nan, **params))
            np.testing.assert_array_equal(none.spec_accels, nan.spec_accels)


def test_usgs_basin_scaling():
    # With basin=True, the Cascadia basin term is removed at periods <= 0.5 s
    # and for Z_2.5 <= 1 km
    params = dict(
        mag=8.0, dist_rup=100.0, v_s30=400.0, event_type="interface", region="cascadia"
    )
    no_basin = PSBAH20(Scenario(**params))
    shallow = PSBAH20(Scenario(depth_2_5=0.8, **params), basin=True)
    np.testing.assert_array_equal(shallow.spec_accels, no_basin.spec_accels)
    deep = PSBAH20(Scenario(depth_2_5=5.0, **params), basin=True)
    short = deep.periods <= 0.5
    np.testing.assert_array_equal(deep.spec_accels[short], no_basin.spec_accels[short])
    assert np.all(deep.spec_accels[~short] > no_basin.spec_accels[~short])


def test_epistemic_branches():
    s = Scenario(
        mag=np.array([7.0, 8.5, 6.5, 7.5]),
        dist_rup=np.array([60.0, 30.0, 100.0, 200.0]),
        depth_tor=np.array([np.nan, np.nan, 50.0, 90.0]),
        v_s30=np.array([300.0, 760.0, 200.0, 1000.0]),
        event_type=np.array(["interface", "interface", "intraslab", "intraslab"]),
        region=np.array(["global", "alaska", "cascadia", "prvi"]),
    )
    epi = PSBAH20(s)
    center = PSBAH20(s, epistemic=False)
    assert epi.epistemic and not center.epistemic
    ln_epi = np.column_stack([epi.ln_pga, np.log(epi.spec_accels)])
    ln_center = np.column_stack([center.ln_pga, np.log(center.spec_accels)])
    assert np.all(ln_epi > ln_center)
    # Global interface: sigma_eps = 0.40 for all periods
    offset = math.log(
        0.185 * math.exp(-1.645 * 0.4) + 0.63 + 0.185 * math.exp(1.645 * 0.4)
    )
    np.testing.assert_allclose(ln_epi[0] - ln_center[0], offset, rtol=1e-12)
    # The standard deviation is the same for the branches
    np.testing.assert_array_equal(epi.ln_stds, center.ln_stds)


# Vectorized scenarios
######################

KEYS = ["pga", "ln_pga", "ln_std_pga", "pgv", "ln_std_pgv", "spec_accels", "ln_stds"]
PERIODS = [0.05, 0.3, 1.0, 2.5]


def scalar_results(*args, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = PSBAH20(Scenario(**dict(zip(NAMES, args))), **kwds)
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


@pytest.mark.parametrize("options", [OPTIONS[0], OPTIONS[3]])
def test_vectorized_matches_scalar(grid, options):
    rows, scenario = grid
    # Every 7th row covers all of the event types and regions
    rows = rows[::7]
    scenario = grid_scenario(rows)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = PSBAH20(scenario, **options)
    expected = [scalar_results(*row, **options) for row in rows]
    for key in KEYS + ["interp_spec_accels", "interp_ln_stds"]:
        if key.startswith("interp"):
            actual = getattr(m, key)(PERIODS)
        else:
            actual = getattr(m, key)
        desired = np.array([e[key] for e in expected])
        assert_matches(actual, desired, key)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_vectorized_shapes(grid):
    rows, scenario = grid
    m = PSBAH20(scenario)
    n = len(rows)
    assert m.pga.shape == (n,)
    assert m.ln_pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.spec_accels.shape == (n, 21)
    assert m.ln_stds.shape == (n, 21)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


def test_scalar_shapes_are_unchanged():
    for event_type, region in itertools.product(EVENT_TYPES, REGIONS):
        m = PSBAH20(
            Scenario(
                mag=7.0,
                dist_rup=80.0,
                depth_tor=50.0,
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


def test_default_region_is_global():
    params = dict(mag=8.0, dist_rup=100.0, v_s30=400.0, event_type="interface")
    np.testing.assert_array_equal(
        PSBAH20(Scenario(**params)).spec_accels,
        PSBAH20(Scenario(region="global", **params)).spec_accels,
    )


def test_scalars_broadcast_with_arrays():
    # The standard deviation does not depend on the magnitude
    mag = np.array([6.0, 7.0, 8.0])
    m = PSBAH20(
        Scenario(
            mag=mag,
            dist_rup=300.0,
            depth_tor=60.0,
            v_s30=300.0,
            event_type="intraslab",
            region="alaska",
        )
    )
    assert m.ln_stds.shape == (3, 21)
    for i in range(3):
        expected = scalar_results(
            mag[i], 300.0, 60.0, 300.0, None, "intraslab", "alaska"
        )
        assert m.pga[i] == expected["pga"]
        np.testing.assert_array_equal(m.spec_accels[i], expected["spec_accels"])
        np.testing.assert_array_equal(m.ln_stds[i], expected["ln_stds"])


def test_mixed_event_types_and_regions():
    event_type = np.array(["interface", "intraslab", "intraslab", "interface"])
    region = np.array(["cascadia", "prvi", "global", "alaska"])
    mag = np.array([9.0, 7.0, 6.0, 8.0])
    depth_tor = np.array([np.nan, 50.0, 120.0, np.nan])
    depth_2_5 = np.array([4.0, np.nan, 2.0, np.nan])
    m = PSBAH20(
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
        m9=True,
    )
    for i in range(4):
        expected = scalar_results(
            mag[i],
            100.0,
            depth_tor[i],
            360.0,
            depth_2_5[i],
            event_type[i],
            region[i],
            basin=True,
            m9=True,
        )
        np.testing.assert_array_equal(m.spec_accels[i], expected["spec_accels"])
        assert m.pga[i] == expected["pga"]
    assert np.all(np.isfinite(m.spec_accels))


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [6.0, 7.0, 8.0], [20.0, 50.0, 100.0, 300.0], indexing="ij"
    )
    m = PSBAH20(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            depth_tor=60.0,
            v_s30=400.0,
            event_type="intraslab",
            region="cascadia",
        )
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    assert m.ln_stds.shape == (3, 4, 21)
    expected = scalar_results(8.0, 50.0, 60.0, 400.0, None, "intraslab", "cascadia")
    assert m.pga[2, 1] == expected["pga"]


def test_interface_depth_can_be_none():
    params = dict(
        mag=np.array([7.0, 8.5]), dist_rup=np.array([30.0, 90.0]), v_s30=500.0
    )
    m = PSBAH20(Scenario(event_type="interface", **params))
    expected = PSBAH20(
        Scenario(event_type="interface", depth_tor=np.array([10.0, 500.0]), **params)
    )
    np.testing.assert_array_equal(m.spec_accels, expected.spec_accels)


@pytest.mark.parametrize(
    "event_type", ["intraslab", np.array(["interface", "intraslab"])]
)
def test_missing_depth_raises(event_type):
    with pytest.raises(ValueError, match="depth_tor is required"):
        PSBAH20(Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type=event_type))


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_event_types():
    # The event type has no default, so results for invalid entries are NaN
    with pytest.warns(UserWarning, match="2 of 4 values"):
        m = PSBAH20(
            Scenario(
                mag=7.0,
                dist_rup=50.0,
                depth_tor=60.0,
                v_s30=400.0,
                event_type=np.array(["interface", "deep", "intraslab", "crustal"]),
            )
        )
    assert np.all(np.isnan(m.spec_accels[[1, 3]]))
    assert np.all(np.isnan(m.ln_stds[[1, 3]]))
    assert np.all(np.isnan(m.pga[[1, 3]]))
    assert np.all(np.isfinite(m.spec_accels[[0, 2]]))
    args = (7.0, 50.0, 60.0, 400.0, None)
    assert m.pga[0] == scalar_results(*args, "interface")["pga"]
    assert m.pga[2] == scalar_results(*args, "intraslab")["pga"]
    assert np.isnan(scalar_results(*args, "deep")["pga"])


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_regions_use_global():
    with pytest.warns(UserWarning, match="1 of 3 values"):
        m = PSBAH20(
            Scenario(
                mag=8.0,
                dist_rup=100.0,
                v_s30=400.0,
                event_type="interface",
                region=np.array(["alaska", "japan", "global"]),
            )
        )
    np.testing.assert_array_equal(m.spec_accels[1], m.spec_accels[2])
    assert np.all(m.spec_accels[0] != m.spec_accels[2])


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(5.0, 9.5, n),
        dist_rup=rng.uniform(10.0, 700.0, n),
        depth_tor=rng.uniform(10.0, 150.0, n),
        depth_2_5=np.where(rng.uniform(size=n) < 0.2, np.nan, rng.uniform(0, 8, n)),
        v_s30=rng.uniform(150.0, 1500.0, n),
        event_type=rng.choice(EVENT_TYPES, n),
        region=rng.choice(REGIONS, n),
    )


@pytest.mark.parametrize("options", [OPTIONS[0], OPTIONS[3], OPTIONS[4]])
@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(ims, options):
    s = vector_scenario()
    full = PSBAH20(s, **options)
    m = PSBAH20(s, ims=ims, **options)
    for key in ["pga", "ln_pga", "ln_std_pga"]:
        np.testing.assert_array_equal(getattr(m, key), getattr(full, key), err_msg=key)
    # Only one period is computed
    assert m._ln_resp.shape == (500, 1)
    assert m._ln_std.shape == (500, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(
        mag=7.0, dist_rup=50.0, depth_tor=70.0, v_s30=300.0, event_type="intraslab"
    )
    m = PSBAH20(s, ims=["pga"])
    assert m.pga == PSBAH20(s).pga
    assert isinstance(m.pga, float)


@pytest.mark.parametrize(
    "ims",
    [["psa_all"], ["psa_ngawest2_21"], ["pga", "psa_1p000", "psa_0p200"]],
)
def test_psa_matches_default(ims):
    s = vector_scenario()
    full = PSBAH20(s, basin=True, m9=True)
    m = PSBAH20(s, basin=True, m9=True, ims=ims)
    cols = np.isin(full.periods, m.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])


def test_pgv_only_matches_default():
    s = vector_scenario()
    full = PSBAH20(s)
    m = PSBAH20(s, ims=["pgv"])
    np.testing.assert_array_equal(m.pgv, full.pgv)
    np.testing.assert_array_equal(m.ln_std_pgv, full.ln_std_pgv)
    assert m._ln_resp.shape == (500, 1)


def test_psa_ngawest2_21_is_all_periods():
    m = PSBAH20(vector_scenario(), ims=["psa_ngawest2_21"])
    np.testing.assert_allclose(m.periods, PSBAH20.PERIODS_NGAWEST2_21)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    m = PSBAH20(vector_scenario(), ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
