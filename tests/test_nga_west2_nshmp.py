"""Test the nshmp-lib options of the NGA-West2 models.

The nshmp-lib variants (``GMM_IDS``) of the NGA-West2 models are compared with
the nshmp-lib reference results (``NgaWest2``, ``DeepBasin``,
``TotalTreeActiveCrust``, and ``UsgsPrviActiveCrust`` tests) and with a
deliberately literal, scalar transcription of the nshmp-lib (Java) classes for
scenarios and variants that the reference results do not cover. BSSA14 uses
the published coefficients, which differ from the rounded nshmp-lib
coefficients, so the BSSA14 results are compared with the reference results
with bounds on the differences.
"""

import inspect
import itertools
import math
import os
import time
import warnings

import numpy as np
import pandas as pd
import pytest

import pygmm
from pygmm.ground_motion import nshmp as usgs
from pygmm.model import Scenario

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

pytestmark = pytest.mark.filterwarnings("ignore:.*recommended limit:UserWarning")

ASK = pygmm.AbrahamsonSilvaKamai2014
BSSA = pygmm.BooreStewartSeyhanAtkinson2014
CB = pygmm.CampbellBozorgnia2014
CY = pygmm.ChiouYoungs2014
MODELS = {"ASK": ASK, "BSSA": BSSA, "CB": CB, "CY": CY}

#: Periods of nshmp-lib (PGV, PGA, and the 21 NGA-West2 periods)
NSHMP_PERIODS = np.r_[-1.0, 0.0, ASK.PERIODS_NGAWEST2_21]


def split_gmm_id(cls, gmm):
    """Model options and scenario values of an nshmp-lib Gmm id."""
    params = set(inspect.signature(cls.__init__).parameters)
    values = cls.GMM_IDS[gmm]
    options = {k: v for k, v in values.items() if k in params}
    scenario = {k: v for k, v in values.items() if k not in params}
    return options, scenario


def make(cls, scenario, gmm, **kwds):
    """Model of an nshmp-lib Gmm id, e.g. "ASK_14_BASIN"."""
    options, values = split_gmm_id(cls, gmm)
    if not isinstance(scenario, Scenario):
        scenario = Scenario(**scenario)
    if values:
        scenario = scenario.copy_with(**values)
    return cls(scenario, **options, **kwds)


def variant(cls, scenario, suffix="", **kwds):
    """Model of the nshmp-lib variant with a suffix, e.g. "_BASIN"."""
    name = {ASK: "ASK", BSSA: "BSSA", CB: "CB", CY: "CY"}[cls]
    return make(cls, scenario, f"{name}_14{suffix}", **kwds)


def period_columns(cls, periods):
    """Columns of the periods (0 for PGA and -1 for PGV) in the model results."""
    return np.array(
        [
            np.flatnonzero(np.isclose(cls.PERIODS, p, rtol=1e-9, atol=0))[0]
            for p in periods
        ]
    )


# Active crust logic trees of the 2023 CONUS NSHM (nshmp-lib GmmTotalTree) as
# (model, Gmm id suffix, weight)
TOTAL_TREES = {
    "TOTAL_TREE_CONUS_ACTIVE_CRUST_2023": [
        (cls, "_BASIN", 0.25) for cls in MODELS.values()
    ],
    "TOTAL_TREE_CONUS_ACTIVE_CRUST_2023_LOS_ANGELES": [
        (cls, suffix, weight)
        for suffix, weight in [("", 0.125), ("_BASIN", 0.0625), ("_CYBERSHAKE", 0.0625)]
        for cls in MODELS.values()
    ],
    "TOTAL_TREE_CONUS_ACTIVE_CRUST_2023_SAN_FRANCISCO": [
        (cls, suffix, 0.125) for suffix in ["", "_BASIN"] for cls in MODELS.values()
    ],
}


def combine_tree(scenario, components, periods):
    """Collapsed median and standard deviation of a logic tree at periods.

    As nshmp-lib ``GroundMotions.combine``: ln(sum(w * exp(mu))) and
    sqrt(sum(w * sigma**2)). `periods` has one period per scenario. Also
    returns the fraction of the collapsed median from BSSA14.
    """
    rows = np.arange(len(periods))
    resp = 0.0
    resp_bssa = 0.0
    var = 0.0
    for cls, suffix, w in components:
        m = variant(cls, scenario, suffix)
        cols = period_columns(cls, periods)
        resp = resp + w * np.exp(m._ln_resp[rows, cols])
        if cls is BSSA:
            resp_bssa = resp_bssa + w * np.exp(m._ln_resp[rows, cols])
        var = var + w * m._ln_std[rows, cols] ** 2
    return np.log(resp), np.sqrt(var), resp_bssa / resp


# Atol of the nshmp-lib results, which are written with 10 decimals
ATOL = 5.1e-11

# nshmp-lib reference results
#############################

NSHMP = pd.read_csv(os.path.join(DATA_DIR, "nga_west2_nshmp-nshmp.csv.gz"), comment="#")


def reference_scenario(df):
    """Scenario of the nshmp-lib inputs."""
    return Scenario(
        mag=df["mag"].values,
        dist_jb=df["dist_jb"].values,
        dist_rup=df["dist_rup"].values,
        dist_x=df["dist_x"].values,
        dip=df["dip"].values,
        width=df["width"].values,
        depth_tor=df["depth_tor"].values,
        depth_hyp=df["depth_hyp"].values,
        mechanism=usgs.mechanism_from_rake(df["rake"].values),
        v_s30=df["v_s30"].values,
        depth_1_0=df["depth_1_0"].values,
        depth_2_5=df["depth_2_5"].values,
    )


REFERENCE_GMMS = sorted(g for g in set(NSHMP["gmm"]) if not g.startswith("PRVI"))

# Bounds of the differences of the natural log of the BSSA14 medians from the
# nshmp-lib results, which are due to the coefficients: nshmp-lib has a
# different f_4 at 0.05, 0.075, 0.5, and 1.5 s (e.g., -0.0679 instead of
# -0.062 at 1.5 s, which changes the nonlinear site term by about 10%), and
# coefficients rounded to four digits at the other periods (e.g., c of -0.6037
# instead of -0.60372 at 0.01 s, which changes the linear site term by about
# 3e-5 per unit of ln(V_s30 / 760)). The PGA and PGV coefficients are the same.
BSSA_F4_PERIODS = [0.05, 0.075, 0.5, 1.5]
BSSA_ATOL_F4 = 0.05
BSSA_ATOL = 1e-3


def bssa_bounds(periods):
    """Bounds of the BSSA14 differences of the natural log of the median."""
    periods = np.asarray(periods)
    return np.select(
        [periods <= 0, np.isin(periods, BSSA_F4_PERIODS)],
        [0.0, BSSA_ATOL_F4],
        BSSA_ATOL,
    )


def test_reference_results_coverage():
    counts = NSHMP.groupby("gmm").size()
    for name in ["ASK", "BSSA", "CB", "CY"]:
        # 107 inputs and PGA, 0.02, 0.2, 1, and 3 s
        assert counts[name + "_14_BASE"] == 107 * 5
    for name in ["ASK", "BSSA", "CB", "CY"]:
        # 20 inputs and PGA, 0.02, 0.2, 0.5, 1, 3, 5, and 10 s
        assert counts[name + "_14_BASIN"] == 20 * 8
        assert counts[name + "_14_CYBERSHAKE"] == 20 * 8
    # 107 inputs and PGA, PGV, and the 21 periods
    for gmm in TOTAL_TREES:
        assert counts[gmm] == 107 * 23
    # 107 inputs and PGA and the 21 periods
    assert counts["PRVI_2025_ACTIVE_CRUST"] == 107 * 22
    assert counts["PRVI_2025_ACTIVE_CRUST_ADJUSTED"] == 107 * 22
    assert len(REFERENCE_GMMS) == 15


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("gmm", [g for g in REFERENCE_GMMS if not g.startswith("BSSA")])
def test_nshmp_reference_results(gmm):
    df = NSHMP[NSHMP["gmm"] == gmm]
    s = reference_scenario(df)
    periods = df["period"].values
    median = df["median"].values
    if gmm in TOTAL_TREES:
        ln_resp, ln_std, share = combine_tree(s, TOTAL_TREES[gmm], periods)
        # BSSA14 uses the published coefficients, so the collapsed median
        # differs by its fraction of the median times the BSSA14 differences
        diff = ln_resp - np.log(median)
        tol = ATOL / median + 1e-10
        assert np.all(np.abs(diff) <= 1.03 * share * bssa_bounds(periods) + tol)
        assert np.all(np.abs(diff[periods <= 0]) <= tol[periods <= 0])
        assert np.max(np.abs(diff)) > 1e-4
    else:
        cls = MODELS[gmm.split("_")[0]]
        m = make(cls, s, gmm)
        rows = np.arange(len(df))
        cols = period_columns(cls, periods)
        ln_resp, ln_std = m._ln_resp[rows, cols], m._ln_std[rows, cols]
        np.testing.assert_allclose(np.exp(ln_resp), median, rtol=0, atol=ATOL)
    # BSSA14 has the same standard deviations as nshmp-lib
    np.testing.assert_allclose(ln_std, df["sigma"], rtol=0, atol=ATOL)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("gmm", [g for g in REFERENCE_GMMS if g.startswith("BSSA")])
def test_bssa_nshmp_reference_results(gmm):
    """BSSA14 with the published coefficients and the nshmp-lib results."""
    df = NSHMP[NSHMP["gmm"] == gmm]
    m = make(BSSA, reference_scenario(df), gmm)
    periods = df["period"].values
    rows = np.arange(len(df))
    cols = period_columns(BSSA, periods)
    ln_resp, ln_std = m._ln_resp[rows, cols], m._ln_std[rows, cols]
    diff = ln_resp - np.log(df["median"].values)
    # Rounding of the medians to 10 decimals
    tol = ATOL / df["median"].values + 1e-10
    assert np.all(np.abs(diff) <= bssa_bounds(periods) + tol)
    # The coefficients of PGA and PGV are the same
    assert np.all(np.abs(diff[periods <= 0]) <= tol[periods <= 0])
    assert np.max(np.abs(diff)) > 1e-4
    np.testing.assert_allclose(ln_std, df["sigma"], rtol=0, atol=ATOL)


def test_reference_variants_differ():
    # The reference inputs exercise the basin and CyberShake terms
    def results(gmm):
        return NSHMP[NSHMP["gmm"] == gmm]["median"].values

    for name in ["ASK", "BSSA", "CB", "CY"]:
        assert np.any(
            np.abs(results(name + "_14_BASIN") - results(name + "_14_CYBERSHAKE"))
            > 1e-3
        )
    conus, la, sf = (results(g) for g in TOTAL_TREES)
    assert np.any(np.abs(conus - la) > 1e-3)
    assert np.any(np.abs(conus - sf) > 1e-3)


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_prvi_reference_results():
    """PRVI coefficients from the nshmp-lib PRVI backbone model results.

    The PRVI 2025 active crust backbone medians are the weighted means (0.25
    each) of the natural log medians of ``ASK_14_BASE``, ``BSSA_14_BASE``,
    ``CB_14_BASE``, and ``CY_14_BASE`` (``PRVI_2025_ACTIVE_CRUST``) and of the
    ``_PRVI`` variants (``PRVI_2025_ACTIVE_CRUST_ADJUSTED``). Both have the
    same epistemic branches and standard deviations, so the ratio of their
    medians is the difference of the backbone medians. BSSA14 uses the
    published coefficients, which change the difference of its PRVI and base
    medians by up to a fraction of the bounds of the BSSA14 differences.
    """
    base = NSHMP[NSHMP["gmm"] == "PRVI_2025_ACTIVE_CRUST"].reset_index(drop=True)
    adj = NSHMP[NSHMP["gmm"] == "PRVI_2025_ACTIVE_CRUST_ADJUSTED"].reset_index(
        drop=True
    )
    np.testing.assert_array_equal(base[["index", "period"]], adj[["index", "period"]])
    np.testing.assert_array_equal(base["sigma"], adj["sigma"])
    s = reference_scenario(base)
    periods = base["period"].values
    rows = np.arange(len(base))
    diff = 0.0
    for name in ["ASK", "BSSA", "CB", "CY"]:
        cls = MODELS[name]
        cols = period_columns(cls, periods)
        ln_prvi = make(cls, s, name + "_14_PRVI")._ln_resp
        ln_base = make(cls, s, name + "_14_BASE")._ln_resp
        diff = diff + 0.25 * (ln_prvi[rows, cols] - ln_base[rows, cols])
    observed = np.log(adj["median"].values) - np.log(base["median"].values)
    # Rounding of the medians to 10 decimals
    tol = ATOL / adj["median"].values + ATOL / base["median"].values
    assert np.all(np.abs(observed - diff) <= tol + 0.25 * bssa_bounds(periods))
    assert np.all(np.abs(observed - diff)[periods <= 0] <= tol[periods <= 0])
    assert np.max(np.abs(diff)) > 0.5


# Literal transcription of nshmp-lib
####################################


def java_name(name):
    """Coefficient name without underscores or case (``v_lin`` is ``Vlin``)."""
    return name.replace("_", "").lower()


class JavaCoeffs(dict):
    """Coefficients of a period, accessed with the Java (nshmp-lib) names."""

    ALIASES = {"rho": "rholnpgalny"}

    def __getitem__(self, key):
        key = java_name(key)
        return super().__getitem__(self.ALIASES.get(key, key))


# nshmp-lib names of the CyberShake coefficients
CYBERSHAKE_NAMES = {
    "abrahamson_silva_kamai_2014": dict(a44="a44cy", a45="a45cy"),
    "boore_stewart_seyhan_atkinson_2014": dict(f_6="f6cy", f_7="f7cy", dz1_max="dz1cy"),
    "campbell_bozorgnia_2014": dict(slope="slope_cy"),
    "chiou_youngs_2014": dict(phi_5="phi5cy", phi_6="phi6cy"),
}
MODULES = {
    "abrahamson_silva_kamai_2014": ASK,
    "boore_stewart_seyhan_atkinson_2014": BSSA,
    "campbell_bozorgnia_2014": CB,
    "chiou_youngs_2014": CY,
}


def java_coeffs(name, imt, prvi, prvi_names):
    """Coefficients of the model (the published coefficients for BSSA14).

    The nshmp-lib coefficients of ASK14, CB14, and CY14 are the same as the
    published coefficients, except for the CY14 c_8 (of the directivity
    term, which nshmp-lib does not have).
    """
    cls = MODULES[name]
    (j,) = np.flatnonzero(np.isclose(cls.PERIODS, imt, rtol=1e-9, atol=0))
    c = JavaCoeffs(
        {java_name(k): float(cls.COEFF[k][j]) for k in cls.COEFF.dtype.names}
    )
    for key, java in CYBERSHAKE_NAMES[name].items():
        c[java_name(java)] = float(cls.COEFF_CYBERSHAKE[key][j])
    if prvi:
        names = {java_name(k): k for k in cls.COEFF_PRVI.dtype.names}
        for key in prvi_names:
            c[java_name(key)] = float(cls.COEFF_PRVI[names[java_name(key)]][j])
    return c


PGA, PGV = 0.0, -1.0


def is_sa(imt):
    return imt > 0


def check_basin(imt, z, z_min):
    # GmmUtils.checkBasin: imt.ordinal() > SA0P5.ordinal()
    return not math.isnan(z) and z > z_min and is_sa(imt) and imt > 0.5


def basin_scale(z, upper, lower):
    constrained = min(max(z, upper), lower)
    return (constrained - upper) / (lower - upper)


def delta_z1_scale(imt, z1p0):
    if check_basin(imt, z1p0, 0.3):
        z_scale = basin_scale(z1p0, 0.3, 0.5)
        return z_scale * 0.585 if imt == 0.75 else z_scale
    return 0.0


def cybershake_imt(imt):
    return is_sa(imt) and imt > 1.9


def fault_style(mechanism):
    return {"U": "UNKNOWN", "SS": "STRIKE_SLIP", "NS": "NORMAL", "RS": "REVERSE"}[
        mechanism
    ]


def epi_tree(mu, Mw, rJB):
    # GroundMotions.createNgaTree and combine
    v = [0.37, 0.25, 0.4, 0.22, 0.23, 0.36, 0.22, 0.23, 0.33]
    values = [[v[0], v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]]
    mi = 0 if Mw < 6 else 1 if Mw < 7 else 2
    ri = 0 if rJB < 10 else 1 if rJB < 30 else 2
    e = values[ri][mi]
    mus = [mu - e, mu, mu + e]
    wts = [0.185, 0.63, 0.185]
    return math.log(sum(w * math.exp(m) for w, m in zip(wts, mus)))


def java_find_y(xs, ys, x):
    # Interpolator.findY with binaryIndex
    import bisect

    if x in xs:
        i = xs.index(x)
    else:
        ip = bisect.bisect_left(xs, x)
        i = 0 if ip == 0 else ip - 1
    if i >= len(xs) - 1:
        i -= 1
    x1, y1, x2, y2 = xs[i], ys[i], xs[i + 1], ys[i + 1]
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1)


# AbrahamsonEtAl_2014


def java_ask(imt, inp, usgs_basin, vs_inferred, cybershake, prvi):
    c = java_coeffs("abrahamson_silva_kamai_2014", imt, prvi, ["a1", "a10"])
    Mw, rJB, rRup, rX = inp["Mw"], inp["rJB"], inp["rRup"], inp["rX"]
    dip, zTor, vs30 = inp["dip"], inp["zTor"], inp["vs30"]
    C4, A3, A4, A5, M2, N = 4.5, 0.275, -0.1, -0.41, 5.0, 1.5
    A = math.pow(610, 4)
    B = math.pow(1360, 4) + A

    c4mag = C4 if Mw > 5 else C4 - (C4 - 1.0) * (5.0 - Mw) if Mw > 4 else 1.0
    R = math.sqrt(rRup * rRup + c4mag * c4mag)
    MaxMwSq = (8.5 - Mw) * (8.5 - Mw)
    MwM1 = Mw - c["M1"]
    f1 = c["a1"] + c["a17"] * rRup
    if Mw > c["M1"]:
        f1 += A5 * MwM1 + c["a8"] * MaxMwSq + (c["a2"] + A3 * MwM1) * math.log(R)
    elif Mw >= M2:
        f1 += A4 * MwM1 + c["a8"] * MaxMwSq + (c["a2"] + A3 * MwM1) * math.log(R)
    else:
        M2M1 = M2 - c["M1"]
        MaxM2Sq = (8.5 - M2) * (8.5 - M2)
        MwM2 = Mw - M2
        f1 += (
            A4 * M2M1
            + c["a8"] * MaxM2Sq
            + c["a6"] * MwM2
            + (c["a2"] + A3 * M2M1) * math.log(R)
        )

    f4 = 0.0
    if rJB < 30 and rX >= 0.0 and Mw > 5.5 and zTor <= 10.0:
        T1 = (90.0 - dip) / 45 if dip > 30.0 else 1.33333333
        dM = Mw - 6.5
        T2 = 1 + 0.2 * dM if Mw >= 6.5 else 1 + 0.2 * dM - (1 - 0.2) * dM * dM
        T3 = 0.0
        r1 = inp["width"] * math.cos(math.radians(dip))
        r2 = 3 * r1
        if rX <= r1:
            rXr1 = rX / r1
            T3 = 0.25 + 1.5 * rXr1 + -0.75 * rXr1 * rXr1
        elif rX <= r2:
            T3 = 1 - (rX - r1) / (r2 - r1)
        T4 = 1 - (zTor * zTor) / 100.0
        T5 = 1.0 if rJB == 0.0 else 1 - rJB / 30.0
        f4 = c["a13"] * T1 * T2 * T3 * T4 * T5

    f6 = c["a15"]
    if zTor < 20.0:
        f6 *= zTor / 20.0

    style = fault_style(inp["mechanism"])
    f78 = 0.0
    if style == "NORMAL":
        f78 = c["a12"] if Mw > 5.0 else c["a12"] * (Mw - 4) if Mw >= 4.0 else 0.0

    # calcSoilTerm
    z1p0 = inp["z1p0"]
    if math.isnan(z1p0):
        f10 = 0.0
    else:
        vsPow4 = vs30 * vs30 * vs30 * vs30
        z1ref = math.exp(-7.67 / 4.0 * math.log((vsPow4 + A) / B)) / 1000.0
        if cybershake:
            vsCoeff = [c["a43"], c["a44cy"], c["a45cy"], c["a46"], c["a46"]]
        else:
            vsCoeff = [c["a43"], c["a44"], c["a45"], c["a46"], c["a46"]]
        z1c = java_find_y([150.0, 250.0, 400.0, 700.0, 1000.0], vsCoeff, vs30)
        z1c *= math.log((z1p0 + 0.01) / (z1ref + 0.01))
        if cybershake:
            z1c += 0.1
        f10 = z1c
    if usgs_basin:
        f10 *= delta_z1_scale(imt, z1p0)

    # getV1
    v1 = 1500.0
    if is_sa(imt):
        if imt >= 3.0:
            v1 = 800.0
        elif imt > 0.5:
            v1 = math.exp(-0.35 * math.log(imt / 0.5) + math.log(1500.0))
    vs30s = vs30 if vs30 < v1 else v1

    saRock = 0.0
    if vs30 < c["Vlin"]:
        vs30s_rk = 1180.0 if 1180.0 < v1 else v1
        f5_rk = (c["a10"] + c["b"] * N) * math.log(vs30s_rk / c["Vlin"])
        saRock = math.exp(f1 + f78 + f5_rk + f4 + f6)
        f5 = (
            c["a10"] * math.log(vs30s / c["Vlin"])
            - c["b"] * math.log(saRock + c["c"])
            + c["b"] * math.log(saRock + c["c"] * math.pow(vs30s / c["Vlin"], N))
        )
    else:
        f5 = (c["a10"] + c["b"] * N) * math.log(vs30s / c["Vlin"])

    mu = f1 + f78 + f5 + f4 + f6 + f10

    def get_phi_a(Mw, s1, s2):
        return s1 if Mw < 4.0 else s2 if Mw > 6.0 else s1 + ((s2 - s1) / 2) * (Mw - 4.0)

    def get_tau_a(Mw, s3, s4):
        return s3 if Mw < 5.0 else s4 if Mw > 7.0 else s3 + ((s4 - s3) / 2) * (Mw - 5.0)

    if vs_inferred:
        phiAsq = get_phi_a(Mw, c["s1e"], c["s2e"])
    else:
        phiAsq = get_phi_a(Mw, c["s1m"], c["s2m"])
    phiAsq *= phiAsq
    tauB = get_tau_a(Mw, c["s3"], c["s4"])
    phiBsq = phiAsq - 0.16
    if vs30 >= c["Vlin"]:
        dAmp = 0.0
    else:
        dAmp = (-c["b"] * saRock) / (saRock + c["c"]) + (c["b"] * saRock) / (
            saRock + c["c"] * math.pow(vs30 / c["Vlin"], N)
        )
    dAmp_p1 = dAmp + 1.0
    phiSq = phiBsq * dAmp_p1 * dAmp_p1 + 0.16
    tau = tauB * dAmp_p1
    sigma = math.sqrt(phiSq + tau * tau)
    return mu, sigma


# BooreEtAl_2014


def java_bssa(imt, inp, usgs_basin, cybershake, prvi):
    names = ["e0", "e1", "e2", "e3", "c"]
    c = java_coeffs("boore_stewart_seyhan_atkinson_2014", imt, prvi, names)
    cPGA = java_coeffs("boore_stewart_seyhan_atkinson_2014", PGA, prvi, names)
    Mw, rJB, vs30 = inp["Mw"], inp["rJB"], inp["vs30"]
    A = math.pow(570.94, 4)
    B = math.pow(1360, 4) + A
    style = fault_style(inp["mechanism"])

    def source(c):
        Fe = (
            c["e1"]
            if style == "STRIKE_SLIP"
            else c["e3"]
            if style == "REVERSE"
            else c["e2"]
            if style == "NORMAL"
            else c["e0"]
        )
        MwMh = Mw - c["Mh"]
        Fe += (
            c["e4"] * MwMh + c["e5"] * MwMh * MwMh if Mw <= c["Mh"] else c["e6"] * MwMh
        )
        return Fe

    def path(c, R):
        return (c["c1"] + c["c2"] * (Mw - 4.5)) * math.log(R / 1.0) + (
            c["c3"] + 0.0
        ) * (R - 1.0)

    R = math.sqrt(rJB * rJB + cPGA["h"] * cPGA["h"])
    pgaRock = math.exp(source(cPGA) + path(cPGA, R))

    Fe = source(c)
    R = math.sqrt(rJB * rJB + c["h"] * c["h"])
    Fp = path(c, R)
    vsLin = vs30 if vs30 <= c["Vc"] else c["Vc"]
    lnFlin = c["c"] * math.log(vsLin / 760.0)
    f2 = c["f4"] * (
        math.exp(c["f5"] * (min(vs30, 760.0) - 360.0))
        - math.exp(c["f5"] * (760.0 - 360.0))
    )
    lnFnl = 0.0 + f2 * math.log((pgaRock + 0.1) / 0.1)
    Fdz1 = 0.0
    if is_sa(imt) and imt >= 0.65:
        z1p0 = inp["z1p0"]
        if math.isnan(z1p0):
            DZ1 = 0.0
        else:
            vsPow4 = vs30 * vs30 * vs30 * vs30
            z1ref = math.exp(-7.15 / 4.0 * math.log((vsPow4 + A) / B)) / 1000.0
            DZ1 = z1p0 - z1ref
        if cybershake:
            Fdz1 = (c["f6cy"] * DZ1 if DZ1 <= c["dz1cy"] else c["f7cy"]) + 0.1
        else:
            Fdz1 = c["f6"] * DZ1 if DZ1 <= c["f7"] / c["f6"] else c["f7"]
    if usgs_basin:
        Fdz1 *= delta_z1_scale(imt, inp["z1p0"])
    Fs = lnFlin + lnFnl + Fdz1
    mu = Fe + Fp + Fs

    # calcStdDev
    tau = (
        c["tau2"]
        if Mw >= 5.5
        else c["tau1"]
        if Mw <= 4.5
        else c["tau1"] + (c["tau2"] - c["tau1"]) * (Mw - 4.5)
    )
    phi_m = (
        c["phi2"]
        if Mw >= 5.5
        else c["phi1"]
        if Mw <= 4.5
        else c["phi1"] + (c["phi2"] - c["phi1"]) * (Mw - 4.5)
    )
    phi_mr = phi_m
    if rJB > c["R2"]:
        phi_mr += c["dPhiR"]
    elif rJB > c["R1"]:
        phi_mr += c["dPhiR"] * (math.log(rJB / c["R1"]) / math.log(c["R2"] / c["R1"]))
    phi_mrv = phi_mr
    if vs30 <= 225:
        phi_mrv -= c["dPhiV"]
    elif vs30 < 300:
        phi_mrv -= c["dPhiV"] * (math.log(300 / vs30) / math.log(300 / 225))
    sigma = math.sqrt(phi_mrv * phi_mrv + tau * tau)
    return mu, sigma


# CampbellBozorgnia_2014


def java_cb(imt, inp, usgs_basin, cybershake, prvi):
    c = java_coeffs("campbell_bozorgnia_2014", imt, prvi, ["c0", "c11"])
    cPGA = java_coeffs("campbell_bozorgnia_2014", PGA, prvi, ["c0", "c11"])
    c["imt"], cPGA["imt"] = imt, PGA
    C, N = 1.88, 1.18
    style = fault_style(inp["mechanism"])

    def calc_basin_term(c, z2p5, cybershake):
        if cybershake:
            if z2p5 <= 1.0:
                return c["c14"] * (z2p5 - 1.0)
            elif z2p5 > 3.0:
                return (
                    c["slope_cy"]
                    * math.exp(-0.75)
                    * (1.0 - math.exp(-0.25 * (z2p5 - 3.0)))
                )
            return 0.0
        if z2p5 <= 1.0:
            return c["c14"] * (z2p5 - 1.0)
        elif z2p5 > 3.0:
            return (
                c["c16"]
                * c["k3"]
                * math.exp(-0.75)
                * (1.0 - math.exp(-0.25 * (z2p5 - 3.0)))
            )
        return 0.0

    def basin_response_term(c, vs30, z2p5, usgs_basin, cybershake):
        zRef = math.exp(7.089 - 1.144 * math.log(vs30))
        zRefTerm = calc_basin_term(c, zRef, False)
        if math.isnan(z2p5):
            return zRefTerm
        z2p5Term = calc_basin_term(c, z2p5, cybershake)
        if usgs_basin:
            if check_basin(c["imt"], z2p5, 1.0):
                zScale = basin_scale(z2p5, 1.0, 3.0)
                zScaled = zRefTerm * (1.0 - zScale) + z2p5Term * zScale
                if cybershake:
                    zScaled += 0.1
                return zScaled * 0.585 if c["imt"] == 0.75 else zScaled
            else:
                return zRefTerm
        return z2p5Term

    def calc_mean(c, vs30, z2p5, pgaRock, usgs_basin, cybershake):
        Mw, rRup, rX, dip = inp["Mw"], inp["rRup"], inp["rX"], inp["dip"]
        Fmag = c["c0"] + c["c1"] * Mw
        if Mw > 6.5:
            Fmag += c["c2"] * (Mw - 4.5) + c["c3"] * (Mw - 5.5) + c["c4"] * (Mw - 6.5)
        elif Mw > 5.5:
            Fmag += c["c2"] * (Mw - 4.5) + c["c3"] * (Mw - 5.5)
        elif Mw > 4.5:
            Fmag += c["c2"] * (Mw - 4.5)
        r = math.sqrt(rRup * rRup + c["c7"] * c["c7"])
        Fr = (c["c5"] + c["c6"] * Mw) * math.log(r)
        Fflt = 0.0
        if style == "NORMAL" and Mw > 4.5:
            Fflt = c["c9"]
            if Mw <= 5.5:
                Fflt *= Mw - 4.5
        Fhw = 0.0
        if rX >= 0.0 and Mw > 5.5 and inp["zTor"] <= 16.66:
            r1 = inp["width"] * math.cos(math.radians(dip))
            r2 = 62.0 * Mw - 350.0
            rXr1 = rX / r1
            rXr2r1 = (rX - r1) / (r2 - r1)
            f1_rX = c["h1"] + c["h2"] * rXr1 + c["h3"] * (rXr1 * rXr1)
            f2_rX = 1.0 + c["h5"] * rXr2r1 + c["h6"] * rXr2r1 * rXr2r1
            Fhw_rX = max(f2_rX, 0.0) if rX >= r1 else f1_rX
            Fhw_rRup = 1.0 if rRup == 0.0 else (rRup - inp["rJB"]) / rRup
            Fhw_m = 1.0 + c["a2"] * (Mw - 6.5)
            if Mw <= 6.5:
                Fhw_m *= Mw - 5.5
            Fhw_z = 1.0 - 0.06 * inp["zTor"]
            Fhw_d = (90.0 - dip) / 45.0
            Fhw = c["c10"] * Fhw_rX * Fhw_rRup * Fhw_m * Fhw_z * Fhw_d
        vsk1 = vs30 / c["k1"]
        if vs30 <= c["k1"]:
            Fsite = c["c11"] * math.log(vsk1) + c["k2"] * (
                math.log(pgaRock + C * math.pow(vsk1, N)) - math.log(pgaRock + C)
            )
        else:
            Fsite = (c["c11"] + c["k2"] * N) * math.log(vsk1)
        Fsed = basin_response_term(c, vs30, z2p5, usgs_basin, cybershake)
        zHyp = inp["zHyp"]
        Fhyp = 0.0 if zHyp <= 7.0 else zHyp - 7.0 if zHyp <= 20.0 else 13.0
        if Mw <= 5.5:
            Fhyp *= c["c17"]
        elif Mw <= 6.5:
            Fhyp *= c["c17"] + (c["c18"] - c["c17"]) * (Mw - 5.5)
        else:
            Fhyp *= c["c18"]
        Fdip = (
            0.0
            if Mw > 5.5
            else c["c19"] * (5.5 - Mw) * dip
            if Mw > 4.5
            else c["c19"] * dip
        )
        Fatn = c["c20"] * (rRup - 80.0) if rRup > 80.0 else 0.0
        return Fmag + Fr + Fflt + Fhw + Fsite + Fsed + Fhyp + Fdip + Fatn

    vs30, z2p5, Mw = inp["vs30"], inp["z2p5"], inp["Mw"]
    pgaRock = (
        math.exp(calc_mean(cPGA, 1100.0, 0.398, 0.0, usgs_basin, False))
        if vs30 < c["k1"]
        else 0.0
    )
    mu = calc_mean(c, vs30, z2p5, pgaRock, usgs_basin, cybershake)
    if is_sa(imt) and 0.01 <= imt <= 0.25:
        pgaMean = calc_mean(cPGA, vs30, z2p5, pgaRock, usgs_basin, False)
        mu = max(mu, pgaMean)

    # calcStdDev
    vsk1 = vs30 / c["k1"]
    alpha = (
        c["k2"] * pgaRock * (1 / (pgaRock + C * math.pow(vsk1, N)) - 1 / (pgaRock + C))
        if vs30 < c["k1"]
        else 0.0
    )

    def std_mag_dep(lo, hi):
        return hi + (lo - hi) * (5.5 - Mw)

    if Mw <= 4.5:
        tau_lnYB, phi_lnY = c["tau1"], c["phi1"]
        tau_lnPGAB, phi_lnPGAB = cPGA["tau1"], cPGA["phi1"]
    elif Mw < 5.5:
        tau_lnYB = std_mag_dep(c["tau1"], c["tau2"])
        phi_lnY = std_mag_dep(c["phi1"], c["phi2"])
        tau_lnPGAB = std_mag_dep(cPGA["tau1"], cPGA["tau2"])
        phi_lnPGAB = std_mag_dep(cPGA["phi1"], cPGA["phi2"])
    else:
        tau_lnYB, phi_lnY = c["tau2"], c["phi2"]
        tau_lnPGAB, phi_lnPGAB = cPGA["tau2"], cPGA["phi2"]
    alphaTau = alpha * tau_lnPGAB
    tauSq = (
        tau_lnYB * tau_lnYB
        + alphaTau * alphaTau
        + 2.0 * alpha * c["rho"] * tau_lnYB * tau_lnPGAB
    )
    phi_lnYB = math.sqrt(phi_lnY * phi_lnY - 0.09)
    phi_lnPGAB = math.sqrt(phi_lnPGAB * phi_lnPGAB - 0.09)
    aPhi_lnPGAB = alpha * phi_lnPGAB
    phiSq = (
        phi_lnY * phi_lnY
        + aPhi_lnPGAB * aPhi_lnPGAB
        + 2.0 * c["rho"] * phi_lnYB * aPhi_lnPGAB
    )
    sigma = math.sqrt(phiSq + tauSq)
    return mu, sigma


# ChiouYoungs_2014


def java_cy(imt, inp, usgs_basin, vs_inferred, cybershake, prvi):
    c = java_coeffs("chiou_youngs_2014", imt, prvi, ["c1", "phi1"])
    Mw, rJB, rRup, zTor, vs30 = (
        inp["Mw"],
        inp["rJB"],
        inp["rRup"],
        inp["zTor"],
        inp["vs30"],
    )
    C2, C4, C4A, CRB = 1.06, -2.1, -0.5, 50.0
    dC4 = C4A - C4
    A = math.pow(571, 4)
    B = math.pow(1360, 4) + A
    style = fault_style(inp["mechanism"])

    # calcSAref
    r1 = (
        c["c1"]
        + C2 * (Mw - 6.0)
        + ((C2 - c["c3"]) / c["cn"])
        * math.log(1.0 + math.exp(c["cn"] * (c["cM"] - Mw)))
    )
    r2 = C4 * math.log(rRup + c["c5"] * math.cosh(c["c6"] * max(Mw - c["cHM"], 0.0)))
    gamma = c["cgamma1"] + c["cgamma2"] / math.cosh(max(Mw - c["cgamma3"], 0.0))
    r3 = dC4 * math.log(math.sqrt(rRup * rRup + CRB * CRB)) + rRup * gamma
    coshM = math.cosh(2 * max(Mw - 4.5, 0))
    cosDelta = math.cos(math.radians(inp["dip"]))
    if style == "REVERSE":
        mzTor = 2.704 if Mw <= 5.849 else max(2.704 - 1.226 * (Mw - 5.849), 0)
    else:
        mzTor = 2.673 if Mw <= 4.970 else max(2.673 - 1.136 * (Mw - 4.970), 0)
    dzTor = zTor - mzTor * mzTor
    r4 = (c["c7"] + c["c7b"] / coshM) * dzTor + (
        0.0 + c["c11b"] / coshM
    ) * cosDelta * cosDelta
    if style == "REVERSE":
        r4 += c["c1a"] + c["c1c"] / coshM
    elif style == "NORMAL":
        r4 += c["c1b"] + c["c1d"] / coshM
    r5 = 0.0
    if inp["rX"] >= 0.0:
        r5 = (
            c["c9"]
            * math.cos(math.radians(inp["dip"]))
            * (c["c9a"] + (1.0 - c["c9a"]) * math.tanh(inp["rX"] / c["c9b"]))
            * (1 - math.sqrt(rJB * rJB + zTor * zTor) / (rRup + 1.0))
        )
    saRef = math.exp(r1 + r2 + r3 + r4 + r5)

    # calcSoilNonLin
    exp1 = math.exp(c["phi3"] * (min(vs30, 1130.0) - 360.0))
    exp2 = math.exp(c["phi3"] * (1130.0 - 360.0))
    snl = c["phi2"] * (exp1 - exp2)

    # calcMean
    sl = c["phi1"] * min(math.log(vs30 / 1130.0), 0.0)
    snl_mod = snl * math.log((saRef + c["phi4"]) / c["phi4"])
    z1p0 = inp["z1p0"]
    if math.isnan(z1p0):
        dZ1 = 0.0
    else:
        vsPow4 = vs30 * vs30 * vs30 * vs30
        z1ref = math.exp(-7.15 / 4 * math.log((vsPow4 + A) / B))
        dZ1 = z1p0 * 1000.0 - z1ref
    if cybershake:
        rkdepth = c["phi5cy"] * (1.0 - math.exp(-dZ1 / c["phi6cy"])) + 0.1
    else:
        rkdepth = c["phi5"] * (1.0 - math.exp(-dZ1 / 300.0))
    if usgs_basin:
        rkdepth *= delta_z1_scale(imt, z1p0)
    mu = math.log(saRef) + sl + snl_mod + rkdepth

    # calcStdDev
    NL0 = snl * saRef / (saRef + c["phi4"])
    mTest = min(max(Mw, 5.0), 6.5) - 5.0
    tau = c["tau1"] + (c["tau2"] - c["tau1"]) / 1.5 * mTest
    sigmaNL0 = c["sigma1"] + (c["sigma2"] - c["sigma1"]) / 1.5 * mTest
    vsTerm = c["sigma3"] if vs_inferred else 0.7
    NL0sq = (1 + NL0) * (1 + NL0)
    sigmaNL0 *= math.sqrt(vsTerm + NL0sq)
    sigma = math.sqrt(tau * tau * NL0sq + sigmaNL0 * sigmaNL0)
    return mu, sigma


def java_calc(name, imt, inp, gmm):
    """nshmp-lib ground motion (mean, sigma) of a model and Gmm id."""
    options, values = split_gmm_id(MODELS[name], gmm)
    assert options["nshmp"]
    epistemic = options.get("epistemic", False)
    basin = options.get("basin", False)
    cybershake = options.get("cybershake", False) and cybershake_imt(imt)
    prvi = options.get("prvi", False)
    vs_inferred = values.get("vs_source", "inferred") == "inferred"
    if name == "ASK":
        mu, sigma = java_ask(imt, inp, basin, vs_inferred, cybershake, prvi)
    elif name == "BSSA":
        mu, sigma = java_bssa(imt, inp, basin, cybershake, prvi)
    elif name == "CB":
        mu, sigma = java_cb(imt, inp, basin, cybershake, prvi)
    else:
        mu, sigma = java_cy(imt, inp, basin, vs_inferred, cybershake, prvi)
    if epistemic:
        mu = epi_tree(mu, inp["Mw"], inp["rJB"])
    return mu, sigma


# Scenario grid of the transcription tests
MECHS = ["SS", "NS", "RS", "U"]
GRID_NAMES = [
    "mag",
    "dist_jb",
    "dist_rup",
    "dist_x",
    "dip",
    "width",
    "depth_tor",
    "depth_hyp",
    "mechanism",
    "v_s30",
    "depth_1_0",
    "depth_2_5",
]
JAVA_NAMES = [
    "Mw",
    "rJB",
    "rRup",
    "rX",
    "dip",
    "width",
    "zTor",
    "zHyp",
    "mechanism",
    "vs30",
    "z1p0",
    "z2p5",
]


def make_grid():
    rng = np.random.default_rng(1)
    rows = []
    geometries = [
        # rJB, rRup, rX, dip, width, zTor, zHyp
        (0.0, 0.0, 0.0, 90.0, 10.0, 0.0, 5.0),
        (0.0, 2.0, 3.0, 45.0, 15.0, 2.0, 7.5),
        (5.0, 7.0, 12.0, 30.0, 20.0, 3.0, 8.0),
        (2.0, 6.0, -2.0, 60.0, 12.0, 5.0, 12.0),
        (15.0, 16.0, 40.0, 20.0, 25.0, 9.0, 14.0),
        (25.0, 26.0, 30.0, 70.0, 12.0, 12.0, 18.0),
        (40.0, 41.0, 50.0, 45.0, 10.0, 17.0, 22.0),
        (95.0, 95.2, -95.0, 90.0, 14.0, 21.0, 25.0),
        (250.0, 250.1, 300.0, 50.0, 18.0, 1.0, 6.0),
    ]
    sites = [
        # vs30, z1p0, z2p5
        (120.0, 0.8, 4.0),
        (180.0, np.nan, np.nan),
        (200.0, 0.05, 0.5),
        (260.0, 0.35, 1.5),
        (300.0, 0.45, 2.5),
        (400.0, 0.6, 3.5),
        (550.0, 2.0, 7.0),
        (760.0, np.nan, 1.0),
        (900.0, 0.3, np.nan),
        (1100.0, 0.02, 0.2),
        (1600.0, 0.1, 0.4),
    ]
    mags = [3.5, 4.5, 5.2, 5.8, 6.3, 6.75, 6.9, 7.4, 8.2]
    for mag, geom, site in itertools.product(mags, geometries, sites):
        mech = MECHS[rng.integers(4)]
        rjb, rrup, rx, dip, width, ztor, zhyp = geom
        vs30, z1, z25 = site
        rows.append((mag, rjb, rrup, rx, dip, width, ztor, zhyp, mech, vs30, z1, z25))
    return rows


def grid_scenario(rows):
    columns = list(zip(*rows))
    values = {
        name: np.array(col) if name == "mechanism" else np.array(col, dtype=float)
        for name, col in zip(GRID_NAMES, columns)
    }
    return Scenario(**values)


@pytest.fixture(scope="module")
def grid():
    rows = make_grid()
    return rows, grid_scenario(rows)


VARIANTS = [(name, gmm) for name, cls in MODELS.items() for gmm in cls.GMM_IDS]


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("name,gmm", VARIANTS, ids=[v[1] for v in VARIANTS])
def test_matches_nshmp_transcription(grid, name, gmm):
    rows, scenario = grid
    cls = MODELS[name]
    m = make(cls, scenario, gmm)
    # Every 3rd row (with the nshmp-lib IMTs) keeps the test fast
    index = np.arange(0, len(rows), 3)
    cols = period_columns(cls, NSHMP_PERIODS)
    expected = np.array(
        [
            [
                java_calc(name, imt, dict(zip(JAVA_NAMES, rows[i])), gmm)
                for imt in NSHMP_PERIODS
            ]
            for i in index
        ]
    )
    np.testing.assert_allclose(
        m._ln_resp[index][:, cols], expected[..., 0], rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(
        m._ln_std[index][:, cols], expected[..., 1], rtol=1e-12, atol=1e-12
    )


def test_grid_exercises_branches(grid):
    rows, scenario = grid
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        # Variants differ on the grid
        for name in ["ASK", "BSSA", "CB", "CY"]:
            cls = MODELS[name]
            cols = period_columns(cls, NSHMP_PERIODS)
            results = {
                gmm: np.concatenate([m._ln_resp[:, cols], m._ln_std[:, cols]])
                for gmm in cls.GMM_IDS
                for m in [make(cls, scenario, gmm)]
            }
            values = list(results.values())
            for a, b in itertools.combinations(range(len(values)), 2):
                assert np.any(np.abs(values[a] - values[b]) > 1e-6)
        # Hanging-wall terms and the CB14 PGA limit at short periods
        assert np.any(np.asarray(scenario.dist_x) < 0)


# Model details
###############


def test_coefficient_files():
    # The nshmp-lib PRVI and CyberShake coefficients are given at the nshmp-lib
    # periods, and are NaN at the other periods of the models
    for cls, prvi_names, cybershake_names in [
        (ASK, {"a1", "a10"}, {"a44", "a45"}),
        (BSSA, {"e_0", "e_1", "e_2", "e_3", "c"}, {"f_6", "f_7", "dz1_max"}),
        (CB, {"c_0", "c_11"}, {"slope"}),
        (CY, {"c_1", "phi_1"}, {"phi_5", "phi_6"}),
    ]:
        nshmp = np.zeros(len(cls.PERIODS), dtype=bool)
        nshmp[period_columns(cls, NSHMP_PERIODS)] = True
        for coeff, names, rows in [
            (cls.COEFF_PRVI, prvi_names, nshmp),
            (cls.COEFF_CYBERSHAKE, cybershake_names, nshmp & (cls.PERIODS > 1.9)),
        ]:
            np.testing.assert_array_equal(coeff["period"], cls.PERIODS)
            assert set(coeff.dtype.names) == names | {"period"}
            # The PRVI coefficients replace the model coefficients
            if coeff is cls.COEFF_PRVI:
                assert names <= set(cls.COEFF.dtype.names)
            for name in names:
                assert np.all(np.isfinite(coeff[name][rows]))
                assert np.all(np.isnan(coeff[name][~rows]))
        # The PRVI coefficients differ from the published coefficients
        for name in prvi_names:
            assert np.all(cls.COEFF_PRVI[name][nshmp] != cls.COEFF[name][nshmp])


def test_epistemic_epsilon():
    eps = usgs.epistemic_epsilon(
        [5.9, 6.0, 6.99, 7.0, 5.0, 6.5, 7.5, 5.0, 6.5, 7.5],
        [9.9, 9.9, 9.9, 9.9, 10.0, 29.9, 10.0, 30.0, 30.0, 300.0],
    )
    np.testing.assert_array_equal(
        eps, [0.37, 0.25, 0.25, 0.40, 0.22, 0.23, 0.36, 0.22, 0.23, 0.33]
    )
    factor = usgs.epistemic_ln_factor(6.5, 20.0)
    expected = math.log(0.185 * math.exp(-0.23) + 0.63 + 0.185 * math.exp(0.23))
    assert math.isclose(factor, expected, rel_tol=1e-15)


def test_mechanism_from_rake():
    rake = [np.nan, -180, -135, -134, -90, -46, -45, -44, 0, 44, 45, 90, 135, 136, 180]
    expected = ["U", "SS", "NS", "NS", "NS", "NS", "NS", "SS", "SS"]
    expected += ["SS", "RS", "RS", "RS", "SS", "SS"]
    assert usgs.mechanism_from_rake(rake).tolist() == expected
    assert usgs.mechanism_from_rake(90.0) == "RS"


def test_basin_scale():
    periods = np.array([-1, 0, 0.5, 0.75, 1.0, 3.0])
    scale = usgs.basin_scale(
        periods, np.array([[np.nan], [0.3], [0.4], [0.6]]), 0.3, 0.5
    )
    np.testing.assert_allclose(
        scale,
        [
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0.5 * 0.585, 0.5, 0.5],
            [0, 0, 0, 0.585, 1, 1],
        ],
        rtol=1e-14,
    )


SCENARIO = dict(
    mag=6.8,
    dist_jb=4.0,
    dist_rup=6.0,
    dist_x=5.0,
    dip=50.0,
    width=16.0,
    depth_tor=1.5,
    depth_hyp=8.0,
    mechanism="RS",
    v_s30=350.0,
    depth_1_0=0.45,
    depth_2_5=2.5,
)


def scenario(**kwds):
    return Scenario(**{**SCENARIO, **kwds})


@pytest.mark.parametrize("cls", [ASK, BSSA, CB, CY])
def test_basin_variant_removes_short_period_basin_term(cls):
    # Z_2.5 from 1 to 3 km has no CB14 basin term
    basin = variant(cls, scenario(depth_2_5=5.0), "_BASIN")
    # No basin depth
    none = variant(cls, scenario(depth_1_0=None, depth_2_5=None), "_BASIN")
    short = basin.periods <= 0.5
    np.testing.assert_array_equal(basin.spec_accels[short], none.spec_accels[short])
    assert basin.pga == none.pga
    assert basin.pgv == none.pgv
    assert np.all(
        basin.spec_accels[basin.periods >= 1] != none.spec_accels[basin.periods >= 1]
    )


@pytest.mark.parametrize("cls", [ASK, BSSA, CY])
def test_basin_none_is_nan(cls):
    for gmm in cls.GMM_IDS:
        none = make(cls, scenario(depth_1_0=None), gmm)
        nan = make(cls, scenario(depth_1_0=np.nan), gmm)
        np.testing.assert_array_equal(none.spec_accels, nan.spec_accels)
        np.testing.assert_array_equal(none.ln_stds, nan.ln_stds)


def test_cb_basin_none_is_nan():
    for gmm in CB.GMM_IDS:
        none = make(CB, scenario(depth_2_5=None), gmm)
        nan = make(CB, scenario(depth_2_5=np.nan), gmm)
        np.testing.assert_array_equal(none.spec_accels, nan.spec_accels)


@pytest.mark.parametrize("cls", [ASK, BSSA, CB, CY])
def test_cybershake_long_periods_only(cls):
    ims = ["pga", "pgv", "psa_ngawest2_21"]
    basin = variant(cls, scenario(), "_BASIN", ims=ims)
    cy = variant(cls, scenario(), "_CYBERSHAKE", ims=ims)
    long = basin.periods > 1.9
    np.testing.assert_array_equal(cy.spec_accels[~long], basin.spec_accels[~long])
    assert cy.pga == basin.pga
    assert np.all(cy.spec_accels[long] != basin.spec_accels[long])


@pytest.mark.parametrize("cls", [ASK, BSSA, CB, CY])
def test_nan_without_nshmp_coefficients(cls):
    # The CyberShake and PRVI coefficients are only given at the nshmp-lib
    # periods
    nshmp = np.zeros(len(cls.PERIODS), dtype=bool)
    nshmp[period_columns(cls, NSHMP_PERIODS)] = True
    for suffix, rows in [
        ("_PRVI", nshmp),
        ("_CYBERSHAKE", nshmp | (cls.PERIODS <= 1.9)),
    ]:
        m = variant(cls, scenario(), suffix)
        psa = rows[cls.INDICES_PSA]
        assert np.all(np.isfinite(m.spec_accels[psa]))
        assert np.all(np.isnan(m.spec_accels[~psa]))
        assert np.isfinite(m.pga) and np.isfinite(m.pgv)
    assert np.all(np.isfinite(variant(cls, scenario(), "_BASIN").spec_accels))


@pytest.mark.parametrize("cls", [ASK, CY])
def test_vs30_measured_only_affects_sigma(cls):
    inferred = variant(cls, scenario(), "_BASE")
    measured = variant(cls, scenario(), "_VS30_MEASURED")
    np.testing.assert_array_equal(measured.spec_accels, inferred.spec_accels)
    assert np.all(measured.ln_stds <= inferred.ln_stds)
    assert np.all(
        measured.ln_stds[measured.periods < 1] < inferred.ln_stds[measured.periods < 1]
    )
    # nshmp-lib uses the variability of an inferred V_s30
    np.testing.assert_array_equal(cls(scenario(), nshmp=True).ln_stds, inferred.ln_stds)


@pytest.mark.parametrize("cls", list(MODELS.values()))
def test_epistemic_branches(cls):
    epi = variant(cls, scenario())
    center = variant(cls, scenario(), "_BASE")
    assert epi.epistemic and not center.epistemic
    # M 6.8, R_JB of 4 km: epsilon = 0.25
    offset = math.log(0.185 * math.exp(-0.25) + 0.63 + 0.185 * math.exp(0.25))
    np.testing.assert_allclose(
        np.log(epi.spec_accels) - np.log(center.spec_accels), offset, rtol=1e-12
    )
    np.testing.assert_array_equal(epi.ln_stds, center.ln_stds)
    # The branches used for hazard (nshmp-lib GroundMotions.createNgaTree)
    for im in ["pga", "psa", "pgv"]:
        branches = epi.ln_branches(im)
        [(w_c, ln_c, std_c)] = center.ln_branches(im)
        assert w_c == 1.0
        assert [w for w, _, _ in branches] == [0.185, 0.63, 0.185]
        for (_, ln, std), sign in zip(branches, [-1, 0, 1]):
            np.testing.assert_allclose(ln, ln_c + sign * 0.25, rtol=0, atol=1e-12)
            np.testing.assert_array_equal(std, std_c)
    w = np.array([b[0] for b in epi.ln_branches()])
    ln = np.array([b[1] for b in epi.ln_branches()])
    np.testing.assert_allclose(np.log(w @ np.exp(ln)), epi.ln_pga, rtol=1e-12)
    # The epistemic branches can also be used with the published model
    published = cls(scenario(), epistemic=True)
    np.testing.assert_allclose(
        np.log(published.spec_accels) - np.log(cls(scenario()).spec_accels),
        offset,
        rtol=1e-12,
    )


def test_width_from_depth_bor_and_depth_hyp():
    params = dict(SCENARIO)
    params.pop("width")
    params.pop("depth_hyp")
    depth_bor = 1.5 + 16.0 * math.sin(math.radians(50.0))
    for cls in [ASK, CB]:
        expected = cls(
            scenario(depth_hyp=1.5 + 8.0 * math.sin(math.radians(50.0))), nshmp=True
        )
        m = cls(Scenario(depth_bor=depth_bor, **params), nshmp=True)
        np.testing.assert_allclose(m.spec_accels, expected.spec_accels, rtol=1e-13)
        # Without depth_bor, the width is estimated as in the published model
        m = cls(Scenario(**params), nshmp=True)
        published = cls(Scenario(**params))
        np.testing.assert_allclose(m.scenario.width, published.scenario.width)
    # The hypocentral depth is the middle of the rupture
    m = CB(Scenario(**params), nshmp=True)
    assert math.isclose(
        m.scenario.depth_hyp,
        1.5 + m.scenario.width * math.sin(math.radians(50.0)) / 2,
        rel_tol=1e-15,
    )
    # The width is only needed by CB14 for the hanging-wall term
    m = CB(scenario(width=None, dist_x=-5.0), nshmp=True)
    np.testing.assert_allclose(
        m.spec_accels, CB(scenario(dist_x=-5.0), nshmp=True).spec_accels, rtol=1e-14
    )


def test_unspecified_mechanism():
    # "U" uses e0 in BSSA14 and is the same as strike-slip in the other models
    for cls in [ASK, CB, CY]:
        np.testing.assert_array_equal(
            cls(scenario(mechanism="U"), nshmp=True).spec_accels,
            cls(scenario(mechanism="SS"), nshmp=True).spec_accels,
        )
        mechanism = np.array(["U", "SS", "NS"])
        m = cls(scenario(mechanism=mechanism), nshmp=True)
        np.testing.assert_array_equal(m.spec_accels[0], m.spec_accels[1])
        assert np.any(m.spec_accels[2] != m.spec_accels[1])
    assert np.all(
        BSSA(scenario(mechanism="U"), nshmp=True).spec_accels
        != BSSA(scenario(mechanism="SS"), nshmp=True).spec_accels
    )


@pytest.mark.parametrize("cls", [ASK, CY])
def test_hanging_wall_from_dist_x(cls):
    hw = cls(scenario(), nshmp=True)
    np.testing.assert_array_equal(
        hw.spec_accels, cls(scenario(on_hanging_wall=True), nshmp=True).spec_accels
    )
    fw = cls(scenario(on_hanging_wall=False), nshmp=True)
    assert np.any(hw.spec_accels != fw.spec_accels)
    # The published models use on_hanging_wall=False if it is not given
    np.testing.assert_array_equal(
        cls(scenario()).spec_accels,
        cls(scenario(on_hanging_wall=False)).spec_accels,
    )


def test_options():
    s = scenario()
    for cls in MODELS.values():
        with pytest.raises(ValueError, match="cybershake=True requires basin=True"):
            cls(s, nshmp=True, cybershake=True)
        with pytest.raises(TypeError):
            cls(s, vs30_measured=True)
        m = cls(s)
        assert not (m.nshmp or m.epistemic or m.basin or m.cybershake or m.prvi)
    m = variant(CY, s, "_CYBERSHAKE")
    assert m.nshmp and m.basin and m.cybershake and m.epistemic and not m.prvi
    assert m.scenario.vs_source == "inferred"


def test_missing_dist_x():
    # Without dist_x, the hanging-wall term is not applied
    params = dict(SCENARIO)
    params.pop("dist_x")
    np.testing.assert_array_equal(
        ASK(Scenario(**params), nshmp=True).spec_accels,
        ASK(scenario(dist_x=-1.0), nshmp=True).spec_accels,
    )


def test_gmm_ids():
    assert set(ASK.GMM_IDS) == {
        "ASK_14",
        "ASK_14_BASE",
        "ASK_14_BASIN",
        "ASK_14_CYBERSHAKE",
        "ASK_14_VS30_MEASURED",
        "ASK_14_PRVI",
    }
    assert set(CY.GMM_IDS) == {g.replace("ASK", "CY") for g in ASK.GMM_IDS}
    for cls in [BSSA, CB]:
        prefix = "BSSA" if cls is BSSA else "CB"
        assert set(cls.GMM_IDS) == {
            prefix + s
            for s in ["_14", "_14_BASE", "_14_BASIN", "_14_CYBERSHAKE", "_14_PRVI"]
        }
    for cls in MODELS.values():
        for gmm in cls.GMM_IDS:
            options, values = split_gmm_id(cls, gmm)
            assert options["nshmp"]
            assert options.get("epistemic", False) == (
                not gmm.endswith(("_BASE", "_VS30_MEASURED", "_PRVI"))
            )
            # The scenario values are vs_source of ASK14 and CY14
            assert set(values) == ({"vs_source"} if cls in [ASK, CY] else set())
            make(cls, scenario(), gmm)


# Comparison with the published models
#######################################


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.filterwarnings("ignore::Warning")
def test_published_models():
    """Compare the _BASE reference results with the published models.

    The differences are explained by the published models' treatment of
    some inputs, coefficient revisions, and nshmp-lib details:

    - ASK14: nshmp-lib uses a hanging-wall dip taper of 1.33333333 instead of
      4 / 3.
    - BSSA14: nshmp-lib uses rounded coefficients (e.g., c_2 at 1.5, 2, and
      3 s, c, f_4, f_6, and f_7) and different f_4 values at 0.05, 0.075,
      0.5, and 1.5 s than the published model (revised 2014-07-15
      coefficients).
    - CB14: the reference rock Z_2.5 is 0.398 km in nshmp-lib, and the PGA
      limit of the short periods uses the reference rock PGA of the period.
    - CY14: none (both use 571 m/sec in the reference Z_1.0, as in the paper).
    """
    inp = NSHMP[NSHMP["gmm"] == "ASK_14_BASE"].drop_duplicates("index")
    inp = inp.sort_values("index").reset_index(drop=True)
    mech = usgs.mechanism_from_rake(inp["rake"].values)
    v_s30 = inp["v_s30"].values.astype(float)
    z1, z25 = inp["depth_1_0"].values, inp["depth_2_5"].values
    geom = dict(
        mag=inp["mag"].values,
        dist_jb=inp["dist_jb"].values,
        dist_rup=inp["dist_rup"].values,
        dist_x=inp["dist_x"].values,
        dip=inp["dip"].values.astype(float),
        depth_tor=inp["depth_tor"].values.astype(float),
        v_s30=v_s30,
    )
    on_hw = inp["dist_x"].values >= 0
    # nshmp-lib has no basin term for NaN Z_1.0, the same as the reference
    # depth of the published models
    bssa_z1_ref = (
        np.exp(-7.15 / 4 * np.log((v_s30**4 + 570.94**4) / (1360.0**4 + 570.94**4)))
        / 1000
    )
    published = {
        "ASK": pygmm.AbrahamsonSilvaKamai2014(
            Scenario(
                **geom,
                width=inp["width"].values,
                mechanism=np.where(mech == "U", "SS", mech),
                depth_1_0=np.where(
                    np.isnan(z1),
                    pygmm.AbrahamsonSilvaKamai2014.calc_depth_1_0(v_s30),
                    z1,
                ),
                region="california",
                vs_source="inferred",
                on_hanging_wall=on_hw,
            )
        ),
        "BSSA": pygmm.BooreStewartSeyhanAtkinson2014(
            Scenario(
                mag=geom["mag"],
                dist_jb=geom["dist_jb"],
                v_s30=v_s30,
                mechanism=mech,
                depth_1_0=np.where(np.isnan(z1), bssa_z1_ref, z1),
                region="california",
            )
        ),
        "CB": pygmm.CampbellBozorgnia2014(
            Scenario(
                **geom,
                width=inp["width"].values,
                depth_hyp=inp["depth_hyp"].values,
                mechanism=np.where(mech == "U", "SS", mech),
                depth_2_5=np.where(
                    np.isnan(z25),
                    pygmm.CampbellBozorgnia2014.calc_depth_2_5(v_s30, "california"),
                    z25,
                ),
                region="california",
            )
        ),
        "CY": pygmm.ChiouYoungs2014(
            Scenario(
                **geom,
                mechanism=mech,
                depth_1_0=np.where(
                    np.isnan(z1),
                    pygmm.ChiouYoungs2014.calc_depth_1_0(v_s30, "california"),
                    z1,
                ),
                region="california",
                vs_source="inferred",
                on_hanging_wall=on_hw,
            )
        ),
    }
    periods = [0.0, 0.02, 0.2, 1.0, 3.0]
    for name, m in published.items():
        df = NSHMP[NSHMP["gmm"] == name + "_14_BASE"]
        ref = df.pivot(index="index", columns="period", values="median")[periods]
        ref_sigma = df.pivot(index="index", columns="period", values="sigma")[periods]
        cols = [np.flatnonzero(np.isclose(m.periods, p))[0] for p in periods[1:]]
        ln_resp = np.column_stack([np.log(m.pga), np.log(m.spec_accels[:, cols])])
        ln_std = np.column_stack([m.ln_std_pga, m.ln_stds[:, cols]])
        diff = ln_resp - np.log(ref.values)
        diff_sigma = ln_std - ref_sigma.values
        if name == "ASK":
            # Hanging-wall dip taper of 1.33333333 instead of 4 / 3
            other = np.repeat(
                (inp["dip"].values <= 30)[:, np.newaxis], len(periods), axis=1
            )
            assert np.max(np.abs(diff[other])) < 2e-9
        elif name == "BSSA":
            # Rounded coefficients (e.g., c_2 at 3 s and c)
            other = np.ones(diff.shape, dtype=bool)
            assert np.max(np.abs(diff)) < 2e-4
        elif name == "CB":
            # Short-period PGA limit at T = 0.2 s, V_s30 = 760 m/sec (> k_1)
            other = np.zeros(diff.shape, dtype=bool)
            other[(v_s30 == 760) & (inp["mag"].values == 8.0), 2] = True
            assert np.max(np.abs(diff[other])) < 1e-3
            # Reference rock Z_2.5
            assert np.max(np.abs(diff[~other])) < 3e-6
            other[:] = True
        else:
            # CY14: no differences
            other = np.zeros(diff.shape, dtype=bool)
        # Rounding of the reference medians to 10 decimals
        tol = ATOL / ref.values + 1e-10
        assert np.all(np.abs(diff[~other]) <= tol[~other]), name
        np.testing.assert_allclose(diff_sigma, 0.0, atol=2e-7)


# Vectorized scenarios
######################

KEYS = ["pga", "ln_pga", "ln_std_pga", "spec_accels", "ln_stds", "pgv", "ln_std_pgv"]
PERIODS = [0.05, 0.3, 1.0, 2.5]


def scalar_results(cls, row, gmm):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = make(cls, dict(zip(GRID_NAMES, row)), gmm)
    results = {key: getattr(m, key) for key in KEYS}
    results["interp_spec_accels"] = m.interp_spec_accels(PERIODS)
    results["interp_ln_stds"] = m.interp_ln_stds(PERIODS)
    return results


def assert_matches(actual, desired, key):
    assert actual.shape == desired.shape, key
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different round-off
        np.testing.assert_allclose(actual, desired, rtol=1e-12, err_msg=key)
    else:
        np.testing.assert_array_equal(actual, desired, err_msg=key)


ALL_CLASSES = list(MODELS.values())


@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_vectorized_matches_scalar(grid, cls):
    rows, _ = grid
    rows = rows[::11]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = variant(cls, grid_scenario(rows), "_BASIN")
    gmm = next(g for g in cls.GMM_IDS if g.endswith("_14_BASIN"))
    expected = [scalar_results(cls, row, gmm) for row in rows]
    for key in expected[0]:
        if key.startswith("interp"):
            actual = getattr(m, key)(PERIODS)
        else:
            actual = getattr(m, key)
        desired = np.array([e[key] for e in expected])
        assert_matches(actual, desired, key)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_vectorized_shapes(grid, cls):
    rows, s = grid
    m = variant(cls, s)
    n = len(rows)
    n_psa = len(cls.INDICES_PSA)
    assert m.pga.shape == (n,)
    assert m.ln_pga.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.spec_accels.shape == (n, n_psa)
    assert m.ln_stds.shape == (n, n_psa)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_scalar_shapes(cls):
    m = variant(cls, scenario())
    assert np.ndim(m.pga) == 0
    assert m.spec_accels.shape == (len(cls.INDICES_PSA),)
    assert m.ln_stds.shape == (len(cls.INDICES_PSA),)


@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_scalars_broadcast_with_arrays(cls):
    mag = np.array([5.5, 6.5, 7.5])
    m = variant(cls, scenario(mag=mag))
    n_psa = len(cls.INDICES_PSA)
    assert m.spec_accels.shape == (3, n_psa)
    assert m.ln_stds.shape == (3, n_psa)
    for i in range(3):
        expected = variant(cls, scenario(mag=mag[i]))
        np.testing.assert_array_equal(m.spec_accels[i], expected.spec_accels)
        np.testing.assert_array_equal(m.ln_stds[i], expected.ln_stds)
        assert m.pga[i] == expected.pga
    # Values that only affect the median
    m = variant(cls, scenario(v_s30=np.array([300.0, 760.0])), ims=["pga"])
    assert m.ln_std_pga.shape == (2,)


@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_multidimensional_inputs(cls):
    mag, dist = np.meshgrid([5.5, 6.5, 7.5], [5.0, 20.0, 60.0, 150.0], indexing="ij")
    m = variant(cls, scenario(mag=mag, dist_jb=dist, dist_rup=dist + 1.0))
    n_psa = len(cls.INDICES_PSA)
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, n_psa)
    assert m.ln_stds.shape == (3, 4, n_psa)
    expected = variant(cls, scenario(mag=7.5, dist_jb=20.0, dist_rup=21.0))
    assert m.pga[2, 1] == expected.pga
    np.testing.assert_array_equal(m.spec_accels[2, 1], expected.spec_accels)


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    dist_jb = rng.uniform(0.0, 200.0, n)
    dip = rng.choice([30.0, 45.0, 60.0, 90.0], n)
    depth_tor = rng.uniform(0.0, 10.0, n)
    return Scenario(
        mag=rng.uniform(5.0, 8.0, n),
        dist_jb=dist_jb,
        dist_rup=np.hypot(dist_jb, depth_tor + rng.uniform(0.0, 5.0, n)),
        dist_x=rng.choice([-1.0, 1.0], n) * (dist_jb + rng.uniform(0.0, 10.0, n)),
        dip=dip,
        depth_tor=depth_tor,
        depth_bor=depth_tor + rng.uniform(5.0, 15.0, n),
        mechanism=rng.choice(MECHS, n),
        v_s30=rng.uniform(180.0, 1000.0, n),
        depth_1_0=np.where(rng.random(n) < 0.3, np.nan, rng.uniform(0.0, 1.0, n)),
        depth_2_5=np.where(rng.random(n) < 0.3, np.nan, rng.uniform(0.0, 6.0, n)),
    )


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(cls, ims):
    s = vector_scenario()
    full = variant(cls, s)
    m = variant(cls, s, ims=ims)
    assert m._ln_resp.shape == (500, 1)
    np.testing.assert_array_equal(m.pga, full.pga)
    np.testing.assert_array_equal(m.ln_pga, full.ln_pga)
    np.testing.assert_array_equal(m.ln_std_pga, full.ln_std_pga)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
@pytest.mark.parametrize(
    "ims", [["psa_0p200"], ["pga", "psa_0p100", "psa_1p000", "psa_5p000"]]
)
@pytest.mark.parametrize("suffix", ["_BASIN", "_CYBERSHAKE", "_PRVI"])
def test_psa_matches_default(cls, ims, suffix):
    s = vector_scenario()
    full = variant(cls, s, suffix)
    m = variant(cls, s, suffix, ims=ims)
    cols = [
        np.flatnonzero(np.isclose(full.periods, m.psa_period(im)))[0]
        for im in ims
        if im != "pga"
    ]
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_psa_ngawest2_21(cls):
    s = vector_scenario(50)
    m = variant(cls, s, "_CYBERSHAKE", ims=["psa_ngawest2_21"])
    full = variant(cls, s, "_CYBERSHAKE")
    cols = period_columns(cls, cls.PERIODS_NGAWEST2_21) - cls.INDICES_PSA[0]
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_allclose(m.periods, cls.PERIODS_NGAWEST2_21)
    assert np.all(np.isfinite(m.spec_accels))


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
def test_pgv_only_matches_default(cls):
    s = vector_scenario(50)
    np.testing.assert_array_equal(variant(cls, s, ims=["pgv"]).pgv, variant(cls, s).pgv)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("cls", ALL_CLASSES, ids=lambda c: c.__name__)
@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(cls, ims):
    m = variant(cls, vector_scenario(50), ims=ims)
    np.testing.assert_allclose(m.ln_pga, np.log(m.pga), rtol=1e-15)


@pytest.mark.slow
@pytest.mark.filterwarnings("ignore::UserWarning")
def test_timing():
    s = vector_scenario(100_000)
    for cls in ALL_CLASSES:
        start = time.perf_counter()
        variant(cls, s, "_BASIN")
        all_ims = time.perf_counter() - start
        start = time.perf_counter()
        variant(cls, s, "_BASIN", ims=["pga"])
        pga = time.perf_counter() - start
        print(f"{cls.__name__}: {all_ims:.2f} s (all), {pga:.2f} s (pga)")
