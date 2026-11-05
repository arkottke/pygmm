"""Test the NGA-East model of USGS nshmp-lib."""

import itertools
import math
import os
import warnings

import numpy as np
import pandas as pd
import pytest

import pygmm
from pygmm.ground_motion import nga_east
from pygmm.ground_motion.nga_east import (
    TABLE_DISTS,
    TABLE_MAGS,
    NgaEast,
    _load_tables,
)
from pygmm.model import Scenario

PERIODS = [0.05, 0.3, 1.0, 2.5]

# Results of nshmp-lib for the 432 inputs of nga-east-inputs.csv and 8 intensity
# measures, given with 10 decimals
NSHMP = pd.read_csv(
    os.path.join(os.path.dirname(__file__), "data", "nga_east-nshmp.csv.gz"),
    comment="#",
)

VARIANTS = list(NgaEast.GMM_IDS.items())


def compute(version="2026", adjusted=False, cpa=False, ims=None, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NgaEast(
            Scenario(**kwds), version=version, adjusted=adjusted, cpa=cpa, ims=ims
        )


@pytest.fixture(scope="module")
def nshmp_scenario():
    return dict(
        mag=NSHMP["mag"].values,
        dist_jb=NSHMP["dist_jb"].values,
        dist_rup=NSHMP["dist_rup"].values,
        v_s30=NSHMP["v_s30"].values,
        depth_sed=NSHMP["depth_sed"].values,
    )


def nshmp_columns():
    period = [0.0 if s == "PGA" else float(s[2:].replace("P", ".")) for s in NSHMP.imt]
    cols = np.searchsorted(NgaEast.PERIODS, period)
    np.testing.assert_array_equal(NgaEast.PERIODS[cols], period)
    return cols


def test_nshmp_inputs():
    # 432 inputs and 8 intensity measures, with 11 inputs on the coastal plain
    assert len(NSHMP) == 432 * 8
    assert NSHMP.loc[NSHMP.depth_sed.notna(), "index"].nunique() == 11


@pytest.mark.parametrize("gmm,options", VARIANTS)
def test_nshmp(nshmp_scenario, gmm, options):
    m = compute(**options, **nshmp_scenario)
    rows = np.arange(len(NSHMP))
    cols = nshmp_columns()
    median = np.exp(m._ln_resp[rows, cols])
    sigma = m._ln_std[rows, cols]
    # The reference values are rounded to 10 decimals
    np.testing.assert_allclose(median, NSHMP["median_" + gmm], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(sigma, NSHMP["sigma"], rtol=0, atol=5.1e-11)


def test_nshmp_variants_differ():
    # The reference results distinguish the versions and options
    on_cpa = NSHMP.depth_sed.notna()
    for a, b in [
        ("2018", "2023"),
        ("2023", "2026"),
        ("2023", "2023_ADJUSTED"),
        ("2026", "2026_ADJUSTED"),
    ]:
        diff = NSHMP[f"median_NGA_EAST_{a}"] - NSHMP[f"median_NGA_EAST_{b}"]
        assert np.max(np.abs(diff)) > 1e-3
    for v in ["2023", "2026"]:
        diff = NSHMP[f"median_NGA_EAST_{v}_CPA"] - NSHMP[f"median_NGA_EAST_{v}"]
        assert np.all(diff[~on_cpa] == 0)
        assert np.all(np.abs(diff[on_cpa & (NSHMP.v_s30 < 3000)]) > 0)


# Scalar transliteration of nshmp-lib NgaEast.NgaEast_2018 and NgaEast_2023, used
# to check the vectorized model at cases that are not in the nshmp-lib results


def ref_position(keys, value):
    i = int(np.searchsorted(keys, value, side="right")) - 1
    i = min(max(i, 0), len(keys) - 2)
    lo, hi = keys[i], keys[i + 1]
    frac = 0.0 if value < lo else 1.0 if value > hi else (value - lo) / (hi - lo)
    return i, frac


def ref_table(table, mag, dist_rup):
    # table has shape (distances, magnitudes)
    ir, fr = ref_position(
        np.log10(TABLE_DISTS),
        math.log10(dist_rup) if dist_rup > 0 else -math.inf,
    )
    im, fm = ref_position(TABLE_MAGS, mag)
    i1 = table[ir, im] + fm * (table[ir, im + 1] - table[ir, im])
    i2 = table[ir + 1, im] + fm * (table[ir + 1, im + 1] - table[ir + 1, im])
    return i1 + fr * (i2 - i1)


def ref_cpa_ratio(j, z, m, r):
    data = nga_east._load_cpa_tables()[..., j]
    i, zf = ref_position(nga_east.CPA_DEPTHS, z)
    jm, mf = ref_position(nga_east.CPA_MAGS, m)
    k, rf = ref_position(nga_east.CPA_DISTS, r)

    def interp(lo, hi, f):
        return lo + f * (hi - lo)

    z1m1 = interp(data[i, jm, k], data[i, jm, k + 1], rf)
    z1m2 = interp(data[i, jm + 1, k], data[i, jm + 1, k + 1], rf)
    z2m1 = interp(data[i + 1, jm, k], data[i + 1, jm, k + 1], rf)
    z2m2 = interp(data[i + 1, jm + 1, k], data[i + 1, jm + 1, k + 1], rf)
    return interp(interp(z1m1, z1m2, mf), interp(z2m1, z2m2, mf), zf)


def ref_site_amp(c, version, period, pga_rock, v_s30):
    if v_s30 >= 3000.0:
        return 0.0, 0.0
    v_s30 = max(v_s30, 150.0)
    wt_scale = (0.767 - 0.1) / (math.log(600.0) - math.log(400.0))
    if v_s30 < 400.0:
        wti = 0.1
    elif v_s30 < 600.0:
        wti = wt_scale * math.log(v_s30 / 400.0) + 0.1
    else:
        wti = 0.767
    wtg = 1.0 - wti
    f760 = c["f760i"] * wti + c["f760g"] * wtg
    f760s = c["f760is"] * wti + c["f760gs"] * wtg
    if v_s30 <= c["V1"]:
        fv = c["c"] * math.log(c["V1"] / 760.0)
    elif v_s30 <= c["V2"]:
        fv = c["c"] * math.log(v_s30 / 760.0)
    elif v_s30 <= 2000.0:
        fv = c["c"] * math.log(c["V2"] / 760.0)
    else:
        f2000 = c["c"] * math.log(c["V2"] / 760.0)
        x1, x2 = math.log(2000.0), math.log(3000.0)
        fv = f2000 + (math.log(v_s30) - x1) * (-f760 - f2000) / (x2 - x1)
    if v_s30 < c["Vf"]:
        st = c["sig_l"] - c["sig_vc"]
        vt = (v_s30 - 200.0) / (c["Vf"] - 200.0)
        fvs = c["sig_l"] - 2.0 * st * vt + st * vt * vt
    elif v_s30 <= c["V2"]:
        fvs = c["sig_vc"]
    elif v_s30 <= 2000.0:
        vt = (v_s30 - c["V2"]) / (2000.0 - c["V2"])
        fvs = c["sig_vc"] + (c["sig_u"] - c["sig_vc"]) * vt * vt
    else:
        fvs = c["sig_u"] * (1.0 - math.log(v_s30 / 2000.0) / math.log(3000.0 / 2000.0))
    f_lin = fv + f760
    s_lin = math.sqrt(fvs**2 + f760s**2)
    v_ref_nl = 3000.0 if period >= 0.4 else 760.0
    v_nl = max(v_s30, 200.0 if version == "2026" else 150.0)
    if version == "2026":
        pga_rock = min(pga_rock, 1.0)
    f4 = c["f4"] if version == "2018" else c["f4"] * 0.5 + c["f4mod"] * 0.5
    rk = math.log((pga_rock + c["f3"]) / c["f3"])
    f_nl = 0.0
    if v_nl < c["Vc"]:
        f2 = f4 * (
            math.exp(c["f5"] * (min(v_nl, v_ref_nl) - 360.0))
            - math.exp(c["f5"] * (v_ref_nl - 360.0))
        )
        f_nl = f2 * rk
    sf2 = 0.0
    if v_nl < 300.0:
        sf2 = c["sig_c"]
    elif v_nl < 1000.0:
        sf2 = c["sig_c"] - c["sig_c"] / math.log(1000.0 / 300.0) * math.log(
            v_nl / 300.0
        )
    return f_lin + f_nl, math.sqrt(s_lin**2 + (sf2 * rk) ** 2)


def ref_sigma(c, mag, v_s30):
    if mag <= 4.5:
        tau = c["t1"]
    elif mag <= 5.0:
        tau = c["t1"] + (c["t2"] - c["t1"]) * (mag - 4.5) / 0.5
    elif mag <= 5.5:
        tau = c["t2"] + (c["t3"] - c["t2"]) * (mag - 5.0) / 0.5
    elif mag <= 6.5:
        tau = c["t3"] + (c["t4"] - c["t3"]) * (mag - 5.5)
    else:
        tau = c["t4"]
    if mag <= 5.0:
        phi_ss = c["ss_a"]
    elif mag <= 6.5:
        phi_ss = c["ss_a"] + (mag - 5.0) * (c["ss_b"] - c["ss_a"]) / 1.5
    else:
        phi_ss = c["ss_b"]
    if v_s30 < 1200.0:
        phi_s2s = c["s2s1"]
    elif v_s30 < 1500.0:
        phi_s2s = c["s2s1"] - (c["s2s1"] - c["s2s2"]) / 300.0 * (v_s30 - 1200.0)
    else:
        phi_s2s = c["s2s2"]
    panel = math.sqrt(tau**2 + phi_ss**2 + phi_s2s**2)

    def interp(m5, m6, m7):
        if mag <= 5.0:
            return c[m5]
        if mag <= 6.0:
            return c[m5] + (mag - 5.0) * (c[m6] - c[m5])
        if mag <= 7.0:
            return c[m6] + (mag - 6.0) * (c[m7] - c[m6])
        return c[m7]

    epri = math.hypot(
        interp("phi_M5", "phi_M6", "phi_M7"), interp("tau_M5", "tau_M6", "tau_M7")
    )
    return epri, panel


def ref_model(version, adjusted, cpa, mag, dist_rup, dist_jb, v_s30, depth_sed):
    """ln median and sigma at all periods of the nshmp-lib logic tree."""
    tables = _load_tables()
    n_r, n_m = len(TABLE_DISTS), len(TABLE_MAGS)
    z_score = 1.0 if version == "2018" else 1.645
    if depth_sed is None or math.isnan(depth_sed):
        z_scale = 0.0
    else:
        z_scale = (1.0 - math.exp(-depth_sed / 0.2)) ** 4
    ln_resps, ln_stds = [], []
    coeff = NgaEast.COEFF
    for j, period in enumerate(NgaEast.PERIODS):
        c = {name: coeff[name][j] for name in coeff.dtype.names}
        nga_adj = c["nga_adj"]
        if v_s30 > 1000.0:
            nga_adj += c["vs30_b"] * math.log(min(v_s30, 2000.0) / 1000.0)
        mu_adj = (1.0 - z_scale) * nga_adj if adjusted else 0.0
        apply_cpa = cpa and z_scale > 0.0
        f_cpa = (
            math.log(ref_cpa_ratio(j, depth_sed, mag, dist_jb)) if apply_cpa else 0.0
        )
        # Branches of the logic tree: 17 median models times 2 sigma models
        mus, wts = [], []
        for i in range(17):
            mu = ref_table(tables[i, :, j].reshape(n_r, n_m), mag, dist_rup) + mu_adj
            pga = math.exp(
                ref_table(
                    tables[i, :, NgaEast.INDEX_PGA].reshape(n_r, n_m), mag, dist_rup
                )
                + mu_adj
            )
            f_s, s_s = ref_site_amp(c, version, period, pga, v_s30)
            if f_s == 0.0:
                mu_site = mu
            else:
                mu_amp = mu + f_s
                mu_site = math.log(
                    0.185 * math.exp(mu_amp + s_s * z_score)
                    + 0.63 * math.exp(mu_amp)
                    + 0.185 * math.exp(mu_amp - s_s * z_score)
                )
            if apply_cpa:
                f_ref, _ = ref_site_amp(c, version, period, pga, 1000.0)
                mu_site += f_cpa - f_ref * z_scale
            mus.append(mu_site)
            wts.append(NgaEast.WEIGHTS[j, i])
        sigmas = ref_sigma(c, mag, v_s30)
        total = 0.0
        var = 0.0
        for mu, w in zip(mus, wts):
            for sigma, ws in zip(sigmas, (0.8, 0.2)):
                total += math.exp(mu) * w * ws
                var += sigma**2 * w * ws
        ln_resps.append(math.log(total))
        ln_stds.append(math.sqrt(var))
    return np.array(ln_resps), np.array(ln_stds)


# Includes V_S30 below 200 and 150 m/s, at the phi_S2S, f760, and fv breaks,
# between 1000 and 2000 m/s and 2000 and 3000 m/s, at and above 3000 m/s,
# magnitudes and distances at and beyond the table limits, a rock PGA above 1 g,
# and sediment thicknesses of 0, in the taper, and beyond the tables
REF_CASES = [
    (3.5, 0.0, 0.0, 100.0, np.nan),
    (4.0, 0.5, 0.5, 150.0, 0.05),
    (4.7, 1500.0, 1499.0, 185.0, 0.3),
    (5.2, 2000.0, 1999.0, 250.0, 0.0),
    (5.75, 7.5, 6.0, 450.0, 0.7),
    (6.3, 33.0, 30.0, 600.0, 1.2),
    (6.5, 1e-6, 0.0, 760.0, 30.0),
    (7.0, 120.0, 118.0, 1300.0, 2.0),
    (7.6, 15.0, 10.0, 1500.0, np.nan),
    (8.2, 3.0, 0.0, 1999.0, 4.0),
    (8.5, 75.0, 70.0, 2500.0, 0.15),
    (9.0, 3000.0, 3000.0, 3000.0, 9.0),
    (6.0, 10.0, 8.0, 3500.0, 0.4),
    (8.0, 0.5, 0.0, 190.0, 0.2),
]


@pytest.mark.parametrize("gmm,options", VARIANTS)
@pytest.mark.parametrize("mag,dist_rup,dist_jb,v_s30,depth_sed", REF_CASES)
def test_reference_implementation(
    gmm, options, mag, dist_rup, dist_jb, v_s30, depth_sed
):
    m = compute(
        **options,
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=dist_jb,
        v_s30=v_s30,
        depth_sed=depth_sed,
    )
    kwds = {"adjusted": False, "cpa": False, **options}
    ln_resp, ln_std = ref_model(
        kwds["version"],
        kwds["adjusted"],
        kwds["cpa"],
        mag,
        dist_rup,
        dist_jb,
        v_s30,
        depth_sed,
    )
    np.testing.assert_allclose(m._ln_resp, ln_resp, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(m._ln_std, ln_std, rtol=1e-12)


def test_rock_pga_cap_is_exercised():
    # The rock PGA of the 2026 nonlinear site term is capped at 1 g
    a = compute(version="2023", mag=8.0, dist_rup=0.5, v_s30=190.0)
    b = compute(version="2026", mag=8.0, dist_rup=0.5, v_s30=190.0)
    assert np.exp(_load_tables()[:, 0 * 11 + 9, NgaEast.INDEX_PGA]).max() > 1.0
    assert not np.allclose(a._ln_resp, b._ln_resp)


def test_z_site_scale():
    np.testing.assert_array_equal(nga_east.z_site_scale([np.nan, 0.0]), [0.0, 0.0])
    np.testing.assert_allclose(
        nga_east.z_site_scale([0.2, 1.0, 25.0]),
        [(1 - math.exp(-1)) ** 4, (1 - math.exp(-5)) ** 4, 1.0],
    )


def test_depth_sed_none_nan_and_zero():
    # depth_sed of None, NaN, and 0 are not on the coastal plain
    kwds = dict(mag=6.0, dist_rup=20.0, dist_jb=18.0, v_s30=400.0)
    for options in [dict(adjusted=True), dict(cpa=True), dict()]:
        none = compute(**options, **kwds)
        for depth_sed in [np.nan, 0.0]:
            m = compute(**options, depth_sed=depth_sed, **kwds)
            np.testing.assert_array_equal(m._ln_resp, none._ln_resp)
            np.testing.assert_array_equal(m._ln_std, none._ln_std)
    # The CPA model without a coastal plain site is the base model
    np.testing.assert_array_equal(
        compute(cpa=True, **kwds)._ln_resp, compute(**kwds)._ln_resp
    )
    # dist_jb is only needed on the coastal plain
    kwds.pop("dist_jb")
    compute(cpa=True, depth_sed=np.nan, **kwds)
    with pytest.raises(ValueError, match="dist_jb is required"):
        compute(cpa=True, depth_sed=np.array([np.nan, 0.5]), **kwds)


@pytest.mark.parametrize("version", ["2023", "2026"])
def test_adjustment_taper(version):
    # Without site amplification, the adjustment is added to the median
    depth_sed = np.array([np.nan, 0.0, 0.05, 0.2, 0.5, 1.0, 3.0])
    kwds = dict(mag=6.0, dist_rup=20.0, v_s30=3000.0, depth_sed=depth_sed)
    base = compute(version=version, **kwds)
    adj = compute(version=version, adjusted=True, **kwds)
    scale = nga_east.z_site_scale(depth_sed)[:, np.newaxis]
    c = NgaEast.COEFF
    # V_S30 above 2000 m/s is used as 2000 m/s
    expected = (1 - scale) * (c.nga_adj + c.vs30_b * np.log(2.0))
    np.testing.assert_allclose(adj._ln_resp - base._ln_resp, expected, atol=1e-14)
    # The adjustment decreases with sediment thickness
    diff = np.abs(adj._ln_resp - base._ln_resp)[:, NgaEast.INDEX_PGA]
    assert np.all(np.diff(diff) <= 0)


@pytest.mark.parametrize("version", ["2023", "2026"])
def test_cpa_taper(version):
    depth_sed = np.array([np.nan, 0.0, 0.05, 0.2, 0.5, 1.0, 3.0])
    kwds = dict(mag=6.0, dist_rup=50.0, dist_jb=48.0, v_s30=760.0, depth_sed=depth_sed)
    base = compute(version=version, **kwds)
    cpa = compute(version=version, cpa=True, **kwds)
    diff = cpa._ln_resp - base._ln_resp
    np.testing.assert_array_equal(diff[:2], 0.0)
    assert np.all(np.abs(diff[2:]) > 0)
    for i, z in enumerate(depth_sed):
        ln_resp, _ = ref_model(version, False, True, 6.0, 50.0, 48.0, 760.0, z)
        np.testing.assert_allclose(cpa._ln_resp[i], ln_resp, rtol=1e-12, atol=1e-12)


def test_cpa_ratio_at_table_nodes():
    data = nga_east._load_cpa_tables()
    np.testing.assert_allclose(
        nga_east.cpa_ln_ratio(2.5, 6.0, 300.0), np.log(data[9, 4, 7]), rtol=1e-15
    )
    # Clamped beyond the tables, and NaN sediment thickness is used as 0
    np.testing.assert_array_equal(
        nga_east.cpa_ln_ratio(30.0, 9.0, 2000.0), np.log(data[-1, -1, -1])
    )
    np.testing.assert_array_equal(
        nga_east.cpa_ln_ratio(np.nan, 6.0, 300.0), np.log(data[0, 4, 7])
    )


def test_adjusted_and_cpa_combined():
    kwds = dict(mag=6.0, dist_rup=50.0, dist_jb=48.0, v_s30=400.0, depth_sed=0.3)
    m = compute(adjusted=True, cpa=True, **kwds)
    ln_resp, _ = ref_model("2026", True, True, 6.0, 50.0, 48.0, 400.0, 0.3)
    np.testing.assert_allclose(m._ln_resp, ln_resp, rtol=1e-12, atol=1e-12)


def test_sigma_is_the_same_for_all_versions():
    kwds = dict(mag=6.0, dist_rup=50.0, dist_jb=48.0, v_s30=1300.0, depth_sed=0.3)
    expected = compute(version="2018", **kwds)._ln_std
    for _, options in VARIANTS:
        np.testing.assert_array_equal(compute(**options, **kwds)._ln_std, expected)


@pytest.mark.parametrize(
    "options,match",
    [
        (dict(version="2020"), "version must be one of"),
        (dict(version="2018", adjusted=True), "not available for version 2018"),
        (dict(version="2018", cpa=True), "not available for version 2018"),
    ],
)
def test_invalid_options(options, match):
    with pytest.raises(ValueError, match=match):
        NgaEast(Scenario(mag=6.0, dist_rup=20.0, v_s30=760.0), **options)


def test_version_as_int():
    s = Scenario(mag=6.0, dist_rup=20.0, v_s30=760.0)
    np.testing.assert_array_equal(
        NgaEast(s, version=2023)._ln_resp, NgaEast(s, version="2023")._ln_resp
    )
    assert NgaEast(s).version == "2026"


def test_out_of_range_values_warn():
    with pytest.warns(UserWarning, match="v_s30"):
        NgaEast(Scenario(mag=6.0, dist_rup=20.0, v_s30=np.array([100.0, 760.0])))
    with pytest.warns(UserWarning, match="depth_sed"):
        NgaEast(Scenario(mag=6.0, dist_rup=20.0, v_s30=760.0, depth_sed=30.0))


# Vectorized model


MAGS = [4.0, 6.75, 8.5]
DISTS = [0.0, 75.0, 1500.0]
V_S30S = [140.0, 260.0, 1300.0, 3000.0]
DEPTHS = [np.nan, 0.3, 5.0]


def scalar_results(options, mag, dist_rup, v_s30, depth_sed):
    m = compute(
        **options,
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=0.9 * dist_rup,
        v_s30=v_s30,
        depth_sed=depth_sed,
    )
    return {
        "pga": m.pga,
        "ln_pga": m.ln_pga,
        "ln_std_pga": m.ln_std_pga,
        "pgv": m.pgv,
        "ln_std_pgv": m.ln_std_pgv,
        "spec_accels": m.spec_accels,
        "ln_stds": m.ln_stds,
        "interp_spec_accels": m.interp_spec_accels(PERIODS),
        "interp_ln_stds": m.interp_ln_stds(PERIODS),
    }


GRID_ROWS = list(itertools.product(MAGS, DISTS, V_S30S, DEPTHS))


@pytest.fixture(scope="module", params=[v[0] for v in VARIANTS])
def grid(request):
    options = NgaEast.GMM_IDS[request.param]
    mag, dist_rup, v_s30, depth_sed = (np.array(c) for c in zip(*GRID_ROWS))
    m = compute(
        **options,
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=0.9 * dist_rup,
        v_s30=v_s30,
        depth_sed=depth_sed,
    )
    expected = [scalar_results(options, *row) for row in GRID_ROWS]
    return m, expected


@pytest.mark.parametrize(
    "key",
    [
        "pga",
        "ln_pga",
        "ln_std_pga",
        "pgv",
        "ln_std_pgv",
        "spec_accels",
        "ln_stds",
        "interp_spec_accels",
        "interp_ln_stds",
    ],
)
def test_vectorized_matches_scalar(grid, key):
    m, expected = grid
    if key.startswith("interp"):
        actual = getattr(m, key)(PERIODS)
    else:
        actual = getattr(m, key)
    expected = np.array([e[key] for e in expected])
    assert actual.shape == expected.shape
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, expected, rtol=1e-12)
    else:
        np.testing.assert_array_equal(actual, expected)


def test_vectorized_shapes(grid):
    m, _ = grid
    n = len(GRID_ROWS)
    assert m.pga.shape == (n,)
    assert m.pgv.shape == (n,)
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 21)
    assert m.ln_stds.shape == (n, 21)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


@pytest.mark.parametrize("gmm,options", VARIANTS)
def test_scalar_shapes_are_unchanged(gmm, options):
    m = compute(
        **options, mag=6.5, dist_rup=20.0, dist_jb=19.0, v_s30=760.0, depth_sed=2.0
    )
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert isinstance(m.pgv, float)
    assert m.spec_accels.shape == (21,)
    assert m.ln_stds.shape == (21,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


@pytest.mark.parametrize("options", [dict(adjusted=True), dict(cpa=True)])
def test_scalars_broadcast_with_arrays(options):
    mag = np.array([5.5, 6.5, 7.5])
    dist_rup = np.array([5.0, 20.0, 80.0])

    def scalar(mg, d, v, z):
        return scalar_results(options, mg, d, v, z)["spec_accels"]

    m = compute(
        **options,
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=0.9 * dist_rup,
        v_s30=400.0,
        depth_sed=0.3,
    )
    expected = [scalar(mg, d, 400.0, 0.3) for mg, d in zip(mag, dist_rup)]
    np.testing.assert_array_equal(m.spec_accels, expected)

    # Only V_S30 is an array
    v_s30 = np.array([200.0, 760.0, 2500.0])
    m = compute(
        **options, mag=6.0, dist_rup=20.0, dist_jb=18.0, v_s30=v_s30, depth_sed=0.3
    )
    assert m.spec_accels.shape == (3, 21)
    assert m.ln_stds.shape == (3, 21)
    np.testing.assert_array_equal(
        m.spec_accels, [scalar(6.0, 20.0, v, 0.3) for v in v_s30]
    )

    # Only the sediment thickness is an array
    depth_sed = np.array([np.nan, 0.1, 0.5, 2.0])
    m = compute(
        **options,
        mag=6.0,
        dist_rup=20.0,
        dist_jb=18.0,
        v_s30=400.0,
        depth_sed=depth_sed,
    )
    assert m.spec_accels.shape == (4, 21)
    assert m.ln_stds.shape == (4, 21)
    np.testing.assert_array_equal(
        m.spec_accels, [scalar(6.0, 20.0, 400.0, z) for z in depth_sed]
    )

    # Only the magnitude is an array, the standard deviation has the full shape
    m = compute(**options, mag=mag, dist_rup=20.0, dist_jb=18.0, v_s30=760.0)
    assert m.ln_stds.shape == (3, 21)


@pytest.mark.parametrize("options", [dict(), dict(adjusted=True), dict(cpa=True)])
def test_multidimensional_inputs(options):
    mag, depth_sed = np.meshgrid(
        [5.5, 6.5, 7.5], [np.nan, 0.1, 0.5, 2.0], indexing="ij"
    )
    m = compute(
        **options,
        mag=mag,
        dist_rup=20.0,
        dist_jb=18.0,
        v_s30=400.0,
        depth_sed=depth_sed,
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    assert m.ln_stds.shape == (3, 4, 21)
    assert m.pga[2, 1] == scalar_results(options, 7.5, 20.0, 400.0, 0.1)["pga"]


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    depth_sed = rng.uniform(-1.0, 3.0, n)
    return Scenario(
        mag=rng.uniform(4.0, 8.2, n),
        dist_rup=rng.uniform(0.0, 1000.0, n),
        dist_jb=rng.uniform(0.0, 1000.0, n),
        v_s30=rng.uniform(200.0, 3000.0, n),
        depth_sed=np.where(depth_sed < 0, np.nan, depth_sed),
    )


@pytest.mark.parametrize("gmm,options", VARIANTS)
@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(gmm, options, ims):
    s = vector_scenario()
    full = NgaEast(s, **options)
    pga_only = NgaEast(s, ims=ims, **options)
    np.testing.assert_array_equal(pga_only.pga, full.pga)
    np.testing.assert_array_equal(pga_only.ln_std_pga, full.ln_std_pga)
    # Only one period is computed
    assert pga_only._ln_resp.shape == (500, 1)
    assert pga_only._ln_std.shape == (500, 1)


def test_pga_only_scalar_matches_default():
    s = Scenario(mag=6.5, dist_rup=20.0, dist_jb=19.0, v_s30=760.0, depth_sed=0.4)
    for _, options in VARIANTS:
        full = NgaEast(s, **options)
        pga_only = NgaEast(s, ims=["pga"], **options)
        assert pga_only.pga == full.pga
        assert pga_only.ln_std_pga == full.ln_std_pga
        assert isinstance(pga_only.pga, float)


@pytest.mark.parametrize("options", [dict(), dict(adjusted=True), dict(cpa=True)])
@pytest.mark.parametrize(
    "ims", [["psa_all"], ["psa_ngawest2_21"], ["pga", "pgv", "psa_all"]]
)
def test_psa_matches_default(options, ims):
    s = vector_scenario()
    full = NgaEast(s, **options)
    m = NgaEast(s, ims=ims, **options)
    np.testing.assert_array_equal(m.periods, full.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels)
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds)
    np.testing.assert_allclose(
        m.interp_spec_accels(PERIODS), full.interp_spec_accels(PERIODS), rtol=1e-12
    )


@pytest.mark.parametrize("options", [dict(), dict(adjusted=True), dict(cpa=True)])
def test_selected_periods_match_default(options):
    s = vector_scenario()
    full = NgaEast(s, **options)
    m = NgaEast(s, ims=["pgv", "psa_0p200", "psa_1p000"], **options)
    assert m.psa_ims == ["psa_0p200", "psa_1p000"]
    cols = np.searchsorted(full.periods, [0.2, 1.0])
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])
    np.testing.assert_array_equal(m.pgv, full.pgv)
    np.testing.assert_array_equal(m.ln_std_pgv, full.ln_std_pgv)
    with pytest.raises(ValueError, match="pga was not computed"):
        m.pga


@pytest.mark.parametrize(
    "ims,match",
    [
        (["pgd"], "does not provide 'pgd'"),
        (["psa_0p040"], "does not provide 'psa_0p040'"),
        (["sa"], "not a valid intensity measure"),
    ],
)
def test_invalid_ims_raise(ims, match):
    with pytest.raises(ValueError, match=match):
        NgaEast(vector_scenario(), ims=ims)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    s = vector_scenario()
    m = NgaEast(s, cpa=True, ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_array_equal(m.ln_pga, NgaEast(s, cpa=True).ln_pga)


@pytest.mark.parametrize(
    "cls, options",
    [
        (pygmm.NgaEast, dict()),
        (pygmm.NgaEast, dict(version="2026", cpa=True)),
        (pygmm.NgaEast, dict(version="2023", adjusted=True)),
        (pygmm.NgaEastSeeds, dict(version="2026", cpa=True)),
        (pygmm.NgaEastSeeds, dict(version="2026", adjusted=True)),
    ],
)
def test_ln_branches(cls, options):
    # nshmp-lib logic tree: median models x (EPRI 0.8, panel 0.2) standard deviations;
    # collapsing it as GroundMotions.combine gives the model's ln_pga and ln_std_pga
    s = pygmm.Scenario(
        mag=np.array([7.5, 7.0, 5.5]),
        dist_rup=np.array([11.0, 90.0, 10.0]),
        dist_jb=np.array([10.0, 90.0, 9.0]),
        v_s30=760.0,
        depth_sed=0.9,
    )
    m = cls(s, ims=["pga"], **options)
    branches = m.ln_branches("pga")
    n_means = 17 if cls is pygmm.NgaEast else 14
    assert len(branches) == 2 * n_means
    weights = np.array([b[0] for b in branches])
    np.testing.assert_allclose(weights.sum(), 1.0)
    np.testing.assert_allclose(weights[1::2] / weights[0::2], 0.25)
    for _, ln_mean, ln_std in branches:
        assert ln_mean.shape == (3,) and ln_std.shape == (3,)
    np.testing.assert_allclose(
        np.log(sum(w * np.exp(mu) for w, mu, _ in branches)), m.ln_pga
    )
    np.testing.assert_allclose(
        np.sqrt(sum(w * sd**2 for w, _, sd in branches)), m.ln_std_pga
    )
    # the median models span about 2 (NgaEast) and 0.7 to 1.6 (seeds) natural log units
    means = np.array([mu for _, mu, _ in branches])
    assert np.all(means.max(axis=0) - means.min(axis=0) > 0.5)


def test_ln_branches_period_weights():
    # the weights of the 17 median models depend on the period
    m = pygmm.NgaEast(pygmm.Scenario(mag=6.5, dist_rup=20.0, v_s30=760.0))
    psa = m.ln_branches("psa")
    assert np.shape(psa[0][0]) == np.shape(m.periods)
    np.testing.assert_allclose(sum(w for w, _, _ in psa), 1.0)
    pga = m.ln_branches("pga")
    assert np.ndim(pga[0][0]) == 0
    assert not np.allclose([w for w, _, _ in psa][0], pga[0][0])
