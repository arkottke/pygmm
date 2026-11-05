"""Test the nshmp-lib Atkinson and Macias (2009) model.

The model is compared with the nshmp-lib reference results (``interface`` and
``deep-basin`` tests), the shared conditional PGV model with the nshmp-lib
``UsgsPgvSupport`` reference results, and both with a deliberately literal,
scalar transcription of ``AtkinsonMacias_2009``, ``BooreAtkinson_2008``,
``CampbellBozorgnia_2014.deepBasinScaling``, ``InterpolatedGmm``, and
``UsgsPgvSupport`` in nshmp-lib (Java) for scenarios and intensity measures
that the reference results do not cover. The transcription is also used by
the Zhao et al. (2006) tests.
"""

import csv
import itertools
import math
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from pygmm.ground_motion.atkinson_macias_2009 import (
    PERIODS_SA,
    conditional_pgv,
    deep_basin_term,
    interpolate_periods,
)
from pygmm.ground_motion.atkinson_macias_2009 import AtkinsonMacias2009 as AM09
from pygmm.ground_motion.campbell_bozorgnia_2014 import CampbellBozorgnia2014
from pygmm.ground_motion.zhao_et_al_2006 import bchydro_2012
from pygmm.model import Scenario

DATA_DIR = os.path.join(os.path.dirname(__file__), "data")
COEFF_DIR = os.path.join(os.path.dirname(__file__), "..", "src", "pygmm", "data")

PGA, PGV = 0.0, -1.0
SA_IMTS = [float(p) for p in PERIODS_SA]

# nshmp-lib reference results
#############################

NSHMP = pd.read_csv(
    os.path.join(DATA_DIR, "atkinson_macias_2009-nshmp.csv.gz"), comment="#"
)


def period_columns(cls, periods):
    return np.array(
        [np.flatnonzero(np.isclose(cls.PERIODS, p, atol=0))[0] for p in periods]
    )


def test_reference_results_coverage():
    assert set(NSHMP["gmm"]) == set(AM09.GMM_IDS)
    assert len(NSHMP) == 37 * 4 + 4 * 20 * 8


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("gmm", sorted(AM09.GMM_IDS))
def test_nshmp_reference_results(gmm):
    df = NSHMP[NSHMP["gmm"] == gmm]
    scenario = Scenario(
        mag=df["mag"].values,
        dist_rup=df["dist_rup"].values,
        v_s30=df["v_s30"].values,
        depth_2_5=df["depth_2_5"].values,
    )
    m = AM09(scenario, **AM09.GMM_IDS[gmm])
    rows = np.arange(len(df))
    cols = period_columns(AM09, df["period"].values)
    median = np.exp(m._ln_resp[rows, cols])
    sigma = m._ln_std[rows, cols]
    # nshmp-lib writes the values with 10 decimals
    np.testing.assert_allclose(median, df["median"], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(sigma, df["sigma"], rtol=0, atol=5.1e-11)


def test_reference_variants_differ():
    def results(gmm):
        return NSHMP[NSHMP["gmm"] == gmm]["median"].values

    for a, b in [
        ("AM_09_INTERFACE_BASIN", "AM_09_INTERFACE_BASIN_M9"),
        ("AM_09_INTERFACE_BASIN", "AM_09_INTERFACE_BASIN_SITE_FIX"),
        ("AM_09_INTERFACE_BASIN_M9", "AM_09_INTERFACE_BASIN_M9_SITE_FIX"),
    ]:
        assert np.any(np.abs(results(a) - results(b)) > 1e-4), (a, b)


def test_nshmp_conditional_pgv_reference_results():
    # UsgsPgvSupport.calcAB20Pgv for BCHYDRO_12_INTERFACE
    df = pd.read_csv(
        os.path.join(DATA_DIR, "usgs_pgv_support-nshmp.csv.gz"), comment="#"
    )
    assert len(df) == 90 and set(df["gmm"]) == {"BCHYDRO_12_INTERFACE"}

    def col(name):
        return df[name].values[:, np.newaxis]

    v = dict(
        mag=col("mag"),
        dist_rup=col("dist_rup"),
        depth_tor=col("depth_tor"),
        v_s30=col("v_s30"),
        depth_2_5=col("depth_2_5"),
        is_slab=np.zeros((len(df), 1), dtype=bool),
    )
    ln_resp, ln_std = interpolate_periods(
        PERIODS_SA, {0.03: (0.02, 0.05)}, lambda base: bchydro_2012(base, v)
    )
    ln_pgv, ln_std_pgv = conditional_pgv(
        PERIODS_SA, ln_resp, ln_std, v["mag"], v["dist_rup"], v["v_s30"]
    )
    np.testing.assert_allclose(np.exp(ln_pgv[:, 0]), df["median"], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(ln_std_pgv[:, 0], df["sigma"], rtol=0, atol=5.1e-11)


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


COEFFS_AM09 = load_java_coeffs("atkinson_macias_2009-nshmp.csv")
COEFFS_BA08 = load_java_coeffs("boore_atkinson_2008-nshmp.csv")
# The CB14 basin coefficients of nshmp-lib (CB14.csv) are the same as the
# published coefficients
COEFFS_CB14 = {
    float(row["period"]): {
        "c14": float(row["c_14"]),
        "c16": float(row["c_16"]),
        "k3": float(row["k_3"]),
    }
    for row in CampbellBozorgnia2014.COEFF
}
COEFFS_BCHYDRO = load_java_coeffs("bchydro_2012-nshmp.csv")

BASE_10_TO_E = math.log(10.0)
LN_G_CM_TO_M = math.log(980.0)


def find_y(x1, y1, x2, y2, x):
    # Interpolator.findY
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1)


def is_sa(imt):
    return imt > 0


def java_deep_basin_scaling(imt, z2p5):
    # CampbellBozorgnia_2014.deepBasinScaling; GmmUtils.checkBasin requires an
    # SA ordinal > SA0P5 (PGA and PGV have lower ordinals)
    if not math.isnan(z2p5) and z2p5 > 1.0 and is_sa(imt) and imt > 0.5:
        constrained = min(max(z2p5, 1.0), 3.0)
        z_scale = (constrained - 1.0) / (3.0 - 1.0)
        c = COEFFS_CB14[imt]
        # calcBasinTerm(c, z2p5, false)
        if z2p5 <= 1.0:
            term = c["c14"] * (z2p5 - 1.0)
        elif z2p5 > 3.0:
            term = (
                c["c16"]
                * c["k3"]
                * math.exp(-0.75)
                * (1.0 - math.exp(-0.25 * (z2p5 - 3.0)))
            )
        else:
            term = 0.0
        basin_amp = term * z_scale
        return basin_amp * 0.585 if imt == 0.75 else basin_amp
    return 0.0


def java_usgs_basin(imt, z2p5, m9):
    fb = java_deep_basin_scaling(imt, z2p5)
    if m9 and z2p5 > 6.0 and (is_sa(imt) and imt > 1.9):
        fb = math.log(2.0)
    return fb


def java_ba08_site(imt, pgaRef, vs30):
    # BooreAtkinson_2008.calcSite
    c = COEFFS_BA08[imt]
    PGAlo, A2, A1, V1, V2, Vref = 0.06, 0.09, 0.03, 180.0, 300.0, 760
    Flin = c["b_lin"] * math.log(vs30 / Vref)
    bnl = 0.0
    if vs30 < Vref:
        if vs30 > V2:
            bnl = c["b2"] * math.log(vs30 / Vref) / math.log(V2 / Vref)
        elif vs30 > V1:
            bnl = (c["b1"] - c["b2"]) * math.log(vs30 / V2) / math.log(V1 / V2) + c[
                "b2"
            ]
        else:
            bnl = c["b1"]
    if pgaRef <= A1:
        Fnl = bnl * math.log(PGAlo / 0.1)
    elif pgaRef <= A2:
        dX = math.log(A2 / A1)
        dY = bnl * math.log(A2 / PGAlo)
        _c = (3.0 * dY - bnl * dX) / (dX * dX)
        d = -(2.0 * dY - bnl * dX) / (dX * dX * dX)
        p = math.log(pgaRef / A1)
        Fnl = bnl * math.log(PGAlo / 0.1) + (_c * p * p) + (d * p * p * p)
    else:
        Fnl = bnl * math.log(pgaRef / 0.1)
    return Flin + Fnl


def java_ab20_pgv(calc, sa_imts, Mw, rRup, vs30):
    """UsgsPgvSupport.calcAB20Pgv; calc(imt) returns (mean, sigma)."""
    target = math.exp(-4.09 + 0.66 * Mw)
    lower = upper = sa_imts[0]
    for sa in sa_imts[1:]:
        upper = sa
        if upper > target:
            break
        lower = upper
    lo, hi = calc(lower), calc(upper)
    try:
        mu = find_y(math.log(lower), lo[0], math.log(upper), hi[0], math.log(target))
        sd = find_y(math.log(lower), lo[1], math.log(upper), hi[1], math.log(target))
    except ZeroDivisionError:
        # Java division of 0 by 0 is NaN
        mu = sd = math.nan
    a1, a2, a3, a4 = 5.39, 0.799, 0.654, 0.479
    a5, a6, a7, a8 = -0.062, -0.359, -0.134, 0.023
    if Mw < 5.0:
        f1 = a2
    elif Mw <= 7.5:
        f1 = a2 + (a3 - a2) * (Mw - 5.0) / 2.5
    else:
        f1 = a3
    lnPgv = (
        a1
        + f1 * mu
        + a4 * (Mw - 6.0)
        + a5 * math.pow(8.5 - Mw, 2)
        + a6 * math.log(rRup + 5.0 * math.exp(0.4 * (Mw - 6.0)))
        + (a7 + a8 * (Mw - 5.0)) * math.log(vs30 / 425)
    )
    sigma = math.sqrt(f1 * f1 * sd * sd + 0.33 * 0.33)
    return lnPgv, sigma


AM09_INTERPOLATED = {
    0.02: (0.01, 0.05),
    0.03: (0.01, 0.05),
    0.075: (0.05, 0.1),
    0.15: (0.1, 0.2),
    0.25: (0.2, 0.3),
    1.5: (1.0, 2.0),
}


def java_am09_mean(c, Mw, rRup):
    h = (Mw * Mw) - (3.1 * Mw) - 14.55
    dM = Mw - 8.0
    gnd = c["c0"] + (c["c3"] * dM) + (c["c4"] * dM * dM)
    r = math.sqrt(rRup * rRup + h * h)
    gnd += c["c1"] * math.log10(r) + c["c2"] * r
    return gnd * BASE_10_TO_E - LN_G_CM_TO_M


def java_am09(imt, basin, m9, site_fix, Mw, rRup, vs30, z2p5):
    """AtkinsonMacias_2009.calc: (mean, sigma)."""

    def calc(t):
        return java_am09(t, basin, m9, site_fix, Mw, rRup, vs30, z2p5)

    if imt in AM09_INTERPOLATED:
        lo, hi = AM09_INTERPOLATED[imt]
        g_lo, g_hi = calc(lo), calc(hi)
        return (
            find_y(lo, g_lo[0], hi, g_hi[0], imt),
            find_y(lo, g_lo[1], hi, g_hi[1], imt),
        )
    if imt == PGV:
        return java_ab20_pgv(calc, SA_IMTS, Mw, rRup, vs30)
    c = COEFFS_AM09[imt]
    sigma = c["sig"] * BASE_10_TO_E
    mu_ref = java_am09_mean(c, Mw, rRup)
    mu_pga = java_am09_mean(COEFFS_AM09[PGA], Mw, rRup)
    site = java_ba08_site(imt, math.exp(mu_pga) if site_fix else mu_pga, vs30)
    mu = mu_ref + site
    fb = java_usgs_basin(imt, z2p5, m9) if basin else 0.0
    mu += fb
    return mu, sigma


BCHYDRO_INTERPOLATED = {0.03: (0.02, 0.05)}


def java_bchydro_mean(c, slab, pgaRock, Mw, rRup, zTor, vs30):
    T3, T4, T5, T9 = 0.1, 0.9, 0.0, 0.4
    C1, C4, C, N = 7.8, 10.0, 1.88, 1.18
    dC1 = -0.3 if slab else c["dC1mid"]
    mCut = C1 + dC1
    t13m = c["t13"] * (10 - Mw) * (10 - Mw)
    fMag = (T4 if Mw <= mCut else T5) * (Mw - mCut) + t13m
    fDepth = c["t11"] * (min(zTor, 120.0) - 60.0) if slab else 0.0
    vsS = min(vs30, 1000.0)
    fFaba = 0.0
    fSite = c["t12"] * math.log(vsS / c["vlin"])
    if vs30 < c["vlin"]:
        fSite += -c["b"] * math.log(pgaRock + C) + c["b"] * math.log(
            pgaRock + C * math.pow((vsS / c["vlin"]), N)
        )
    else:
        fSite += c["b"] * N * math.log(vsS / c["vlin"])
    return (
        c["t1"]
        + T4 * dC1
        + (c["t2"] + (c["t14"] if slab else 0.0) + T3 * (Mw - 7.8))
        * math.log(rRup + C4 * math.exp((Mw - 6.0) * T9))
        + c["t6"] * rRup
        + (c["t10"] if slab else 0.0)
        + fMag
        + fDepth
        + fFaba
        + fSite
    )


def java_bchydro(imt, slab, basin, Mw, rRup, zTor, vs30, z2p5):
    """BcHydro_2012.calc (forearc): (mean, sigma)."""

    def calc(t):
        return java_bchydro(t, slab, basin, Mw, rRup, zTor, vs30, z2p5)

    if imt == PGV:
        return java_ab20_pgv(calc, SA_IMTS, Mw, rRup, vs30)
    if imt in BCHYDRO_INTERPOLATED:
        lo, hi = BCHYDRO_INTERPOLATED[imt]
        g_lo, g_hi = calc(lo), calc(hi)
        return (
            find_y(lo, g_lo[0], hi, g_hi[0], imt),
            find_y(lo, g_lo[1], hi, g_hi[1], imt),
        )
    pgaRock = math.exp(
        java_bchydro_mean(COEFFS_BCHYDRO[PGA], slab, 0.0, Mw, rRup, zTor, 1000.0)
    )
    mu = java_bchydro_mean(COEFFS_BCHYDRO[imt], slab, pgaRock, Mw, rRup, zTor, vs30)
    if basin:
        mu += java_deep_basin_scaling(imt, z2p5)
    return mu, 0.74


# Grid of magnitudes (including T_PGV between the longest periods), distances
# (including 0 km), V_s30 values around the BA08 limits, and Z_2.5 values
# around the basin limits (including no basin)
MAGS = [5.0, 6.8, 8.0, 9.0, 9.4]
DISTS_RUP = [0.0, 15.0, 120.0, 700.0]
V_S30S = [150.0, 180.0, 250.0, 300.0, 500.0, 760.0, 1000.0, 1500.0]
DEPTHS_2_5 = [np.nan, 0.5, 2.0, 4.0, 7.0]
NAMES = ["mag", "dist_rup", "v_s30", "depth_2_5"]

OPTIONS = [
    dict(),
    dict(site_fix=True),
    dict(basin=True),
    dict(basin=True, m9=True),
    dict(basin=True, m9=True, site_fix=True),
]


def grid_scenario(rows):
    return Scenario(**{n: np.array(c) for n, c in zip(NAMES, zip(*rows))})


@pytest.fixture(scope="module")
def grid():
    rows = list(itertools.product(MAGS, DISTS_RUP, V_S30S, DEPTHS_2_5))
    return rows, grid_scenario(rows)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("options", OPTIONS)
def test_matches_nshmp_transcription(grid, options):
    rows, scenario = grid
    m = AM09(scenario, **options)
    kwds = dict(
        basin=options.get("basin", False),
        m9=options.get("m9", False),
        site_fix=options.get("site_fix", False),
    )
    expected = np.array(
        [
            [
                java_am09(imt, Mw=mag, rRup=dist, vs30=vs, z2p5=z2p5, **kwds)
                for imt in AM09.PERIODS
            ]
            for mag, dist, vs, z2p5 in rows
        ]
    )
    np.testing.assert_allclose(m._ln_resp, expected[..., 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(m._ln_std, expected[..., 1], rtol=0, atol=1e-12)


@pytest.mark.parametrize("basin", [False, True])
@pytest.mark.parametrize("slab", [False, True])
def test_bchydro_matches_nshmp_transcription(slab, basin):
    rows = list(
        itertools.product(
            MAGS, DISTS_RUP, [10.0, 70.0, 150.0], V_S30S + [1050.0], [np.nan, 7.0]
        )
    )
    mag, dist_rup, depth_tor, v_s30, depth_2_5 = (np.array(c) for c in zip(*rows))
    v = dict(
        mag=mag[:, None],
        dist_rup=dist_rup[:, None],
        depth_tor=depth_tor[:, None],
        v_s30=v_s30[:, None],
        depth_2_5=depth_2_5[:, None],
        is_slab=np.full((len(rows), 1), slab),
    )
    ln_resp, ln_std = interpolate_periods(
        PERIODS_SA, BCHYDRO_INTERPOLATED, lambda base: bchydro_2012(base, v, basin)
    )
    ln_pgv, ln_std_pgv = conditional_pgv(
        PERIODS_SA, ln_resp, ln_std, v["mag"], v["dist_rup"], v["v_s30"]
    )
    ln_resp = np.column_stack([ln_pgv, ln_resp])
    ln_std = np.column_stack(
        [ln_std_pgv, np.broadcast_to(ln_std, ln_resp[:, 1:].shape)]
    )
    expected = np.array(
        [
            [java_bchydro(imt, slab, basin, *row) for imt in [PGV] + SA_IMTS]
            for row in rows
        ]
    )
    np.testing.assert_allclose(ln_resp, expected[..., 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(ln_std, expected[..., 1], rtol=0, atol=1e-12)


def test_coefficient_files():
    assert list(AM09.COEFF.period) == [0.0] + sorted(
        set(SA_IMTS) - set(AM09.INTERPOLATED)
    )
    np.testing.assert_array_equal(AM09.COEFF_BA08.period, AM09.PERIODS)
    assert AM09.PERIODS[AM09.INDEX_PGA] == 0
    assert AM09.PERIODS[AM09.INDEX_PGV] == -1
    np.testing.assert_array_equal(AM09.PERIODS[AM09.INDICES_PSA], PERIODS_SA)
    # PGA coefficients are used for 0.01 s
    pga, sa01 = AM09.COEFF[0], AM09.COEFF[1]
    assert sa01.period == 0.01
    assert all(pga[k] == sa01[k] for k in ["c0", "c1", "c2", "c3", "c4", "sig"])


# Model behavior
################


def test_gmm_ids():
    assert set(AM09.GMM_IDS) == {
        "AM_09_INTERFACE",
        "AM_09_INTERFACE_BASIN",
        "AM_09_INTERFACE_BASIN_M9",
        "AM_09_INTERFACE_BASIN_SITE_FIX",
        "AM_09_INTERFACE_BASIN_M9_SITE_FIX",
    }
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0, depth_2_5=7.0)
    for gmm, options in AM09.GMM_IDS.items():
        m = AM09(s, **options)
        assert m.basin == ("BASIN" in gmm)
        assert m.m9 == ("M9" in gmm)
        assert m.site_fix == ("SITE_FIX" in gmm)


def test_m9_requires_basin():
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0)
    with pytest.raises(ValueError, match="m9=True requires basin=True"):
        AM09(s, m9=True)


def test_invalid_options_raise():
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0)
    with pytest.raises(TypeError):
        AM09(s, epistemic=False)


def test_site_fix_nonlinear_term():
    # Without the site fix, the nonlinear term uses ln(PGA) < 0.03, so it is
    # the same as for a weak motion reference PGA
    s = Scenario(mag=9.0, dist_rup=np.array([20.0, 400.0]), v_s30=250.0)
    m = AM09(s)
    fixed = AM09(s, site_fix=True)
    # Strong shaking: the nonlinear term reduces the motions
    assert fixed.pga[0] < m.pga[0]
    # Weak shaking (PGA < 0.03 g): no difference
    assert fixed.pga[1] < 0.03
    np.testing.assert_allclose(fixed.pga[1], m.pga[1], rtol=1e-14)
    # No nonlinear term at V_s30 >= 760 m/sec
    s = Scenario(mag=9.0, dist_rup=20.0, v_s30=760.0)
    np.testing.assert_array_equal(
        AM09(s).spec_accels, AM09(s, site_fix=True).spec_accels
    )


def test_basin_term():
    params = dict(mag=9.0, dist_rup=100.0, v_s30=400.0)
    none = AM09(Scenario(**params), basin=True)
    nan = AM09(Scenario(depth_2_5=np.nan, **params), basin=True)
    no_basin = AM09(Scenario(depth_2_5=5.0, **params))
    np.testing.assert_array_equal(none.spec_accels, nan.spec_accels)
    np.testing.assert_array_equal(none.spec_accels, no_basin.spec_accels)
    # Zero for Z_2.5 <= 3 km
    shallow = AM09(Scenario(depth_2_5=3.0, **params), basin=True)
    np.testing.assert_array_equal(shallow.spec_accels, no_basin.spec_accels)
    deep = AM09(Scenario(depth_2_5=5.0, **params), basin=True)
    short = deep.periods <= 0.5
    np.testing.assert_array_equal(deep.spec_accels[short], no_basin.spec_accels[short])
    assert np.all(deep.spec_accels[~short] > no_basin.spec_accels[~short])
    assert deep.pga == no_basin.pga
    # 1.5 s is interpolated between 1 and 2 s
    np.testing.assert_allclose(
        np.log(deep.spec_accels[deep.periods == 1.5]),
        np.mean(np.log(deep.spec_accels[np.isin(deep.periods, [1.0, 2.0])])),
        rtol=1e-14,
    )


def test_m9_adjustment():
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0, depth_2_5=np.array([5.0, 7.0]))
    basin = AM09(s, basin=True)
    m9 = AM09(s, basin=True, m9=True)
    long = basin.periods > 1.9
    # 1.5 s is interpolated between 1 and 2 s
    affected = basin.periods >= 1.5
    np.testing.assert_array_equal(m9.spec_accels[0], basin.spec_accels[0])
    np.testing.assert_array_equal(
        m9.spec_accels[1, ~affected], basin.spec_accels[1, ~affected]
    )
    assert np.all(m9.spec_accels[1, affected] != basin.spec_accels[1, affected])
    no_basin = AM09(s)
    np.testing.assert_allclose(
        m9.spec_accels[1, long] / no_basin.spec_accels[1, long], 2.0, rtol=1e-14
    )


def test_deep_basin_term():
    periods = np.array([0.0, 0.5, 0.75, 1.0, 3.0, 10.0])
    depths = np.array([np.nan, 1.0, 2.0, 3.0, 3.5, 6.5, 9.0])[:, None]
    for m9 in [False, True]:
        expected = [[java_usgs_basin(p, z, m9) for p in periods] for z in depths[:, 0]]
        np.testing.assert_allclose(
            deep_basin_term(periods, depths, m9), expected, rtol=0, atol=1e-15
        )


def test_pgv_uses_conditional_model():
    s = Scenario(mag=8.0, dist_rup=50.0, v_s30=400.0)
    m = AM09(s)
    # T_PGV = exp(-4.09 + 0.66 * 8) = 3.29 s, between 3 and 4 s
    ln_pgv, ln_std_pgv = conditional_pgv(
        PERIODS_SA, np.log(m.spec_accels), m.ln_stds, 8.0, 50.0, 400.0
    )
    np.testing.assert_allclose(np.log(m.pgv), ln_pgv, rtol=1e-14)
    np.testing.assert_allclose(m.ln_std_pgv, ln_std_pgv, rtol=1e-14)
    assert 20 < m.pgv < 60


def test_conditional_pgv_period_limits():
    # T_PGV is extrapolated below 0.01 s, and NaN at 10 s and longer
    periods = np.array([0.01, 0.1, 1.0, 10.0])
    ln_resp = np.log([[0.3, 0.6, 0.2, 0.01]])
    ln_std = np.array([[0.6, 0.7, 0.7, 0.8]])
    mag = np.array([[-1.0], [7.0], [9.7]])
    ln_pgv, ln_std_pgv = conditional_pgv(periods, ln_resp, ln_std, mag, 50.0, 400.0)
    for i, mw in enumerate(mag[:, 0]):

        def calc(t):
            j = list(periods).index(t)
            return float(ln_resp[0, j]), float(ln_std[0, j])

        with np.errstate(invalid="ignore"):
            expected = java_ab20_pgv(calc, list(periods), mw, 50.0, 400.0)
        np.testing.assert_allclose(
            [ln_pgv[i, 0], ln_std_pgv[i, 0]], expected, rtol=1e-14
        )
    assert np.isnan(ln_pgv[2, 0]) and np.isfinite(ln_pgv[:2]).all()


# Vectorized scenarios
######################

KEYS = ["pga", "ln_pga", "ln_std_pga", "pgv", "ln_std_pgv", "spec_accels", "ln_stds"]


def scalar_results(**kwds):
    options = {k: kwds.pop(k) for k in ["basin", "m9", "site_fix"] if k in kwds}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = AM09(Scenario(**kwds), **options)
    return {key: getattr(m, key) for key in KEYS}


@pytest.mark.parametrize("options", [OPTIONS[0], OPTIONS[4]])
def test_vectorized_matches_scalar(grid, options):
    rows, _ = grid
    rows = rows[::7]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = AM09(grid_scenario(rows), **options)
    expected = [scalar_results(**dict(zip(NAMES, row)), **options) for row in rows]
    for key in KEYS:
        desired = np.array([e[key] for e in expected])
        actual = getattr(m, key)
        assert actual.shape == desired.shape, key
        np.testing.assert_array_equal(actual, desired, err_msg=key)


def test_scalar_shapes():
    m = AM09(Scenario(mag=8.0, dist_rup=80.0, v_s30=400.0))
    for key in ["pga", "ln_pga", "ln_std_pga", "pgv", "ln_std_pgv"]:
        assert isinstance(getattr(m, key), float), key
    assert m.spec_accels.shape == (21,)
    assert m.ln_stds.shape == (21,)


def test_scalars_broadcast_with_arrays():
    mag = np.array([7.0, 8.0, 9.0])
    m = AM09(Scenario(mag=mag, dist_rup=100.0, v_s30=300.0))
    # The standard deviation does not depend on the scenario
    assert m.ln_stds.shape == (3, 21)
    for i in range(3):
        expected = scalar_results(mag=mag[i], dist_rup=100.0, v_s30=300.0)
        for key in KEYS:
            np.testing.assert_array_equal(
                getattr(m, key)[i], expected[key], err_msg=key
            )


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [7.0, 8.0, 9.0], [20.0, 50.0, 100.0, 300.0], indexing="ij"
    )
    m = AM09(
        Scenario(mag=mag, dist_rup=dist_rup, v_s30=400.0, depth_2_5=7.0), basin=True
    )
    assert m.pga.shape == (3, 4)
    assert m.pgv.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    expected = scalar_results(
        mag=8.0, dist_rup=50.0, v_s30=400.0, depth_2_5=7.0, basin=True
    )
    assert m.pga[1, 1] == expected["pga"]
    assert m.pgv[1, 1] == expected["pgv"]
    np.testing.assert_array_equal(m.spec_accels[1, 1], expected["spec_accels"])


def test_event_types():
    params = dict(mag=8.0, dist_rup=100.0, v_s30=400.0)
    default = AM09(Scenario(**params))
    interface = AM09(Scenario(event_type="interface", **params))
    np.testing.assert_array_equal(default.spec_accels, interface.spec_accels)
    # The model is for interface events only
    m = AM09(Scenario(event_type=np.array(["interface", "intraslab"]), **params))
    np.testing.assert_array_equal(m.spec_accels[0], default.spec_accels)
    assert np.all(np.isnan(m.spec_accels[1]))
    assert np.all(np.isnan(m.ln_stds[1]))
    assert np.isnan(m.pga[1]) and np.isnan(m.pgv[1])
    assert np.isnan(AM09(Scenario(event_type="intraslab", **params)).pga)


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(5.0, 9.5, n),
        dist_rup=rng.uniform(0.0, 700.0, n),
        depth_2_5=np.where(rng.uniform(size=n) < 0.2, np.nan, rng.uniform(0, 9, n)),
        v_s30=rng.uniform(150.0, 1500.0, n),
    )


@pytest.mark.parametrize(
    "ims",
    [["pga"], "pga", ["pgv"], ["psa_all"], ["pga", "psa_1p500", "psa_0p200", "pgv"]],
)
@pytest.mark.parametrize("options", [OPTIONS[0], OPTIONS[4]])
def test_ims_match_default(ims, options):
    s = vector_scenario()
    full = AM09(s, **options)
    m = AM09(s, ims=ims, **options)
    computed = [ims] if isinstance(ims, str) else ims
    if "pga" in computed:
        for key in ["pga", "ln_pga", "ln_std_pga"]:
            np.testing.assert_array_equal(getattr(m, key), getattr(full, key))
    if "pgv" in computed:
        np.testing.assert_array_equal(m.pgv, full.pgv)
        np.testing.assert_array_equal(m.ln_std_pgv, full.ln_std_pgv)
    if any(im.startswith("psa") for im in computed):
        cols = np.isin(full.periods, m.periods)
        np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
        np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])
    if computed in (["pga"], ["pgv"]):
        # Only one column is computed
        assert m._ln_resp.shape == (500, 1)
        assert m._ln_std.shape == (500, 1)


def test_ln_pga():
    m = AM09(vector_scenario(), ims=["pga"])
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
