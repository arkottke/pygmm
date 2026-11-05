"""Test the nshmp-lib Zhao et al. (2006) model.

The model is compared with the nshmp-lib reference results (``interface``,
``slab``, and ``deep-basin`` tests) and with a deliberately literal, scalar
transcription of ``ZhaoEtAl_2006`` and ``ExtrapolatedGmm`` in nshmp-lib (Java)
for scenarios and intensity measures that the reference results do not cover
(the transcription of the models and helpers used by it is in
``test_atkinson_macias_2009``).
"""

import csv
import itertools
import math
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from pygmm.ground_motion.zhao_et_al_2006 import ZhaoEtAl2006 as Zea06
from pygmm.model import Scenario

from .test_atkinson_macias_2009 import (
    COEFF_DIR,
    DATA_DIR,
    PGV,
    SA_IMTS,
    find_y,
    java_ab20_pgv,
    java_am09,
    java_bchydro,
    java_usgs_basin,
    period_columns,
)

# nshmp-lib reference results
#############################

NSHMP = pd.read_csv(os.path.join(DATA_DIR, "zhao_et_al_2006-nshmp.csv.gz"), comment="#")


def split_gmm_id(gmm):
    options = dict(Zea06.GMM_IDS[gmm])
    event_type = options.pop("event_type")
    return event_type, options


def test_reference_results_coverage():
    assert set(NSHMP["gmm"]) == set(Zea06.GMM_IDS)
    assert len(NSHMP) == 37 * 4 + 75 * 4 + 3 * 20 * 8


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("gmm", sorted(Zea06.GMM_IDS))
def test_nshmp_reference_results(gmm):
    df = NSHMP[NSHMP["gmm"] == gmm]
    event_type, options = split_gmm_id(gmm)
    assert event_type == ("intraslab" if "_SLAB" in gmm else "interface")
    scenario = Scenario(
        mag=df["mag"].values,
        dist_rup=df["dist_rup"].values,
        depth_tor=df["depth_tor"].values,
        v_s30=df["v_s30"].values,
        depth_2_5=df["depth_2_5"].values,
        event_type=event_type,
    )
    m = Zea06(scenario, **options)
    rows = np.arange(len(df))
    cols = period_columns(Zea06, df["period"].values)
    median = np.exp(m._ln_resp[rows, cols])
    sigma = m._ln_std[rows, cols]
    # nshmp-lib writes the values with 10 decimals
    np.testing.assert_allclose(median, df["median"], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(sigma, df["sigma"], rtol=0, atol=5.1e-11)


def test_reference_variants_differ():
    def results(gmm):
        return NSHMP[NSHMP["gmm"] == gmm]["median"].values

    a, b = "ZHAO_06_INTERFACE_BASIN", "ZHAO_06_INTERFACE_BASIN_M9"
    assert np.any(np.abs(results(a) - results(b)) > 1e-4)
    # The 10 s results are extrapolated
    assert 10.0 in set(NSHMP[NSHMP["gmm"] == "ZHAO_06_SLAB_BASIN"]["period"])


# Literal transcription of nshmp-lib
####################################


def load_java_coeffs(name):
    with open(os.path.join(COEFF_DIR, name)) as fp:
        lines = [line for line in fp if not line.startswith("#")]
    coeffs = {}
    for row in csv.DictReader(lines):
        coeffs[float(row.pop("period"))] = {k: float(v) for k, v in row.items()}
    return coeffs


COEFFS = load_java_coeffs("zhao_et_al_2006-nshmp.csv")
LN_G_CM_TO_M = math.log(980.0)

INTERPOLATED = {
    0.02: (0.01, 0.05),
    0.03: (0.01, 0.05),
    0.075: (0.05, 0.1),
    0.75: (0.5, 1.0),
}
EXTRAPOLATED = (7.5, 10.0)


def java_calc_mean(c, slab, fSite, Mw, rRup, zTor):
    rRup = max(rRup, 1.0)
    zTor = min(zTor, 125.0) if slab else 20.0
    hfac = 0.0 if zTor < 15.0 else zTor - 15.0
    m2 = Mw - (6.5 if slab else 6.3)
    if slab:
        afac = c["Ssl"] * math.log(rRup) + c["Ss"]
        xmcor = c["Ps"] * m2 + c["Qs"] * m2 * m2 + c["Ws"]
    else:
        afac = c["Si"]
        xmcor = c["Qi"] * m2 * m2 + c["Wi"]
    r = rRup + c["c"] * math.exp(c["d"] * Mw)
    return (
        c["a"] * Mw
        + c["b"] * rRup
        - math.log(r)
        + c["e"] * hfac
        + afac
        + fSite
        + xmcor
        - LN_G_CM_TO_M
    )


def java_site_term_step(c, vs30):
    return c["C1"] if vs30 >= 600.0 else c["C2"] if vs30 >= 300.0 else c["C3"]


def java_extrapolation_refs(slab, basin, m9):
    """extrapolationRefs: [(calc(imt) -> (mean, sigma), weight)]."""
    if slab:
        # BCHYDRO_12_SLAB or BCHYDRO_12_SLAB_BASIN
        return [
            (lambda imt, a: java_bchydro(imt, True, basin, *a), 1.0),
        ]
    # BCHYDRO_12_INTERFACE(_BASIN(_M9)) and AM_09_INTERFACE(_BASIN(_M9)); the
    # BC Hydro M9 variant does not apply the M9 adjustment
    return [
        (lambda imt, a: java_bchydro(imt, False, basin, *a), 0.5),
        (
            lambda imt, a: java_am09(imt, basin, m9, False, a[0], a[1], a[3], a[4]),
            0.5,
        ),
    ]


def java_calc(imt, slab, basin, m9, Mw, rRup, zTor, vs30, z2p5):
    """ZhaoEtAl_2006.calc: (mean, sigma)."""
    args = (Mw, rRup, zTor, vs30, z2p5)

    def calc(t):
        return java_calc(t, slab, basin, m9, *args)

    if imt == PGV:
        return java_ab20_pgv(calc, SA_IMTS, Mw, rRup, vs30)
    if imt in INTERPOLATED:
        lo, hi = INTERPOLATED[imt]
        g_lo, g_hi = calc(lo), calc(hi)
        return (
            find_y(lo, g_lo[0], hi, g_hi[0], imt),
            find_y(lo, g_lo[1], hi, g_hi[1], imt),
        )
    if imt in EXTRAPOLATED:
        # ExtrapolatedGmm with the common IMT SA5P0
        refs = java_extrapolation_refs(slab, basin, m9)
        mu_ref_common = sigma_ref_common = 0.0
        for ref, weight in refs:
            mu, sigma = ref(5.0, args)
            mu_ref_common += mu * weight
            sigma_ref_common += sigma * weight
        mu_ref_target = sigma_ref_target = 0.0
        for ref, weight in refs:
            mu, sigma = ref(imt, args)
            mu_ref_target += mu * weight
            sigma_ref_target += sigma * weight
        mu, sigma = calc(5.0)
        mu_scale = mu / mu_ref_common
        sigma_scale = sigma / sigma_ref_common
        return mu_ref_target * mu_scale, sigma_ref_target * sigma_scale
    c = COEFFS[imt]
    tau = c["tauS"] if slab else c["tau"]
    sigma = math.sqrt(c["sigma"] * c["sigma"] + tau * tau)
    mu = java_calc_mean(c, slab, java_site_term_step(c, vs30), Mw, rRup, zTor)
    if basin:
        mu += java_usgs_basin(imt, z2p5, m9)
    return mu, sigma


# Grid of magnitudes (including T_PGV beyond 7.5 s), distances (including
# less than 1 km), depths around the depth limits, V_s30 values around the
# site class and BC Hydro limits, and Z_2.5 values around the basin limits
MAGS = [5.0, 6.4, 7.8, 9.0, 9.4]
DISTS_RUP = [0.0, 0.5, 30.0, 300.0]
DEPTHS_TOR = [10.0, 40.0, 130.0]
V_S30S = [200.0, 300.0, 599.0, 600.0, 1050.0]
DEPTHS_2_5 = [np.nan, 2.0, 4.0, 7.0]
EVENT_TYPES = ["interface", "intraslab"]
NAMES = ["mag", "dist_rup", "depth_tor", "v_s30", "depth_2_5", "event_type"]

OPTIONS = [dict(), dict(basin=True), dict(basin=True, m9=True)]


def grid_scenario(rows):
    return Scenario(**{n: np.array(c) for n, c in zip(NAMES, zip(*rows))})


@pytest.fixture(scope="module")
def grid():
    rows = list(
        itertools.product(MAGS, DISTS_RUP, DEPTHS_TOR, V_S30S, DEPTHS_2_5, EVENT_TYPES)
    )
    return rows, grid_scenario(rows)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize("options", OPTIONS)
def test_matches_nshmp_transcription(grid, options):
    rows, scenario = grid
    m = Zea06(scenario, **options)
    basin = options.get("basin", False)
    m9 = options.get("m9", False)
    expected = np.array(
        [
            [
                java_calc(
                    imt,
                    event_type == "intraslab",
                    basin,
                    # nshmp-lib has no intraslab M9 variant
                    m9 and event_type == "interface",
                    mag,
                    dist,
                    depth,
                    vs,
                    z2p5,
                )
                for imt in Zea06.PERIODS
            ]
            for mag, dist, depth, vs, z2p5, event_type in rows
        ]
    )
    np.testing.assert_allclose(m._ln_resp, expected[..., 0], rtol=0, atol=1e-12)
    np.testing.assert_allclose(m._ln_std, expected[..., 1], rtol=0, atol=1e-12)


def test_coefficient_files():
    assert list(Zea06.COEFF.period) == [0.0] + sorted(
        set(SA_IMTS) - set(Zea06.INTERPOLATED) - set(Zea06.EXTRAPOLATED)
    )
    assert Zea06.PERIODS[Zea06.INDEX_PGA] == 0
    assert Zea06.PERIODS[Zea06.INDEX_PGV] == -1
    np.testing.assert_array_equal(Zea06.PERIODS[Zea06.INDICES_PSA], SA_IMTS)
    pga, sa01 = Zea06.COEFF[0], Zea06.COEFF[1]
    assert sa01.period == 0.01
    assert all(pga[k] == sa01[k] for k in Zea06.COEFF.dtype.names[1:])


# Model behavior
################


def test_gmm_ids():
    assert set(Zea06.GMM_IDS) == {
        "ZHAO_06_INTERFACE",
        "ZHAO_06_INTERFACE_BASIN",
        "ZHAO_06_INTERFACE_BASIN_M9",
        "ZHAO_06_SLAB",
        "ZHAO_06_SLAB_BASIN",
    }
    s = Scenario(mag=8.0, dist_rup=100.0, depth_tor=50.0, v_s30=400.0, depth_2_5=7.0)
    for gmm in Zea06.GMM_IDS:
        event_type, options = split_gmm_id(gmm)
        m = Zea06(s.copy_with(event_type=event_type), **options)
        assert m.basin == ("BASIN" in gmm)
        assert m.m9 == ("M9" in gmm)


def test_m9_requires_basin():
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0, event_type="interface")
    with pytest.raises(ValueError, match="m9=True requires basin=True"):
        Zea06(s, m9=True)


def test_invalid_options_raise():
    s = Scenario(mag=9.0, dist_rup=100.0, v_s30=400.0, event_type="interface")
    with pytest.raises(TypeError):
        Zea06(s, site_fix=True)


def test_m9_only_affects_interface_events():
    s = Scenario(
        mag=8.0,
        dist_rup=100.0,
        depth_tor=50.0,
        v_s30=400.0,
        depth_2_5=7.0,
        event_type=np.array(["interface", "intraslab"]),
    )
    basin = Zea06(s, basin=True)
    m9 = Zea06(s, basin=True, m9=True)
    np.testing.assert_array_equal(m9.spec_accels[1], basin.spec_accels[1])
    long = basin.periods > 1.9
    np.testing.assert_array_equal(m9.spec_accels[0, ~long], basin.spec_accels[0, ~long])
    assert np.all(m9.spec_accels[0, long] > basin.spec_accels[0, long])


def test_site_classes():
    params = dict(mag=7.0, dist_rup=50.0, event_type="interface")
    s = Scenario(v_s30=np.array([150.0, 299.0, 300.0, 599.0, 600.0, 1500.0]), **params)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = Zea06(s)
    # The extrapolated periods use the V_s30 scaling of the reference models
    cols = m.periods < 7.5
    np.testing.assert_array_equal(m.spec_accels[0, cols], m.spec_accels[1, cols])
    np.testing.assert_array_equal(m.spec_accels[2, cols], m.spec_accels[3, cols])
    np.testing.assert_array_equal(m.spec_accels[4, cols], m.spec_accels[5, cols])
    assert np.all(m.spec_accels[4, ~cols] != m.spec_accels[5, ~cols])
    np.testing.assert_allclose(
        np.log(m.pga[[0, 2, 4]] / m.pga[4]),
        [
            Zea06.COEFF[0].C3 - Zea06.COEFF[0].C1,
            Zea06.COEFF[0].C2 - Zea06.COEFF[0].C1,
            0,
        ],
        atol=1e-14,
    )


def test_depths():
    params = dict(mag=7.0, dist_rup=100.0, v_s30=400.0)
    # Interface events use a depth of 20 km
    interface = Zea06(
        Scenario(
            event_type="interface", depth_tor=np.array([5.0, 20.0, 80.0]), **params
        )
    )
    np.testing.assert_array_equal(interface.spec_accels[0], interface.spec_accels[2])
    # Intraslab events use min(depth_tor, 125), with no depth term for depths
    # less than 15 km. The extrapolated periods use the BC Hydro (2012) depth
    # term, min(depth_tor, 120) - 60.
    slab = Zea06(
        Scenario(
            event_type="intraslab",
            depth_tor=np.array([5.0, 15.0, 125.0, 200.0]),
            **params,
        )
    )
    cols = slab.periods < 7.5
    np.testing.assert_array_equal(slab.spec_accels[0, cols], slab.spec_accels[1, cols])
    np.testing.assert_array_equal(slab.spec_accels[2], slab.spec_accels[3])
    assert np.all(slab.pga[2] > slab.pga[1])


def test_minimum_distance():
    s = Scenario(
        mag=7.0, dist_rup=np.array([0.0, 0.5, 1.0]), v_s30=400.0, event_type="interface"
    )
    m = Zea06(s, ims=["pga", "psa_1p000"])
    np.testing.assert_array_equal(m.pga[0], m.pga[2])
    np.testing.assert_array_equal(m.spec_accels[1], m.spec_accels[2])


def test_extrapolation():
    s = Scenario(mag=8.0, dist_rup=100.0, v_s30=400.0, event_type="interface")
    m = Zea06(s)
    # 7.5 and 10 s are smaller than 5 s
    sa = dict(zip(m.periods, m.spec_accels))
    assert sa[10.0] < sa[7.5] < sa[5.0]


def test_pgv_uses_extrapolated_periods():
    # T_PGV = exp(-4.09 + 0.66 * 9.4) = 8.3 s, between 7.5 and 10 s
    for event_type in EVENT_TYPES:
        s = Scenario(
            mag=9.4, dist_rup=150.0, depth_tor=60.0, v_s30=400.0, event_type=event_type
        )
        m = Zea06(s)
        expected = java_calc(
            PGV,
            event_type == "intraslab",
            False,
            False,
            9.4,
            150.0,
            60.0,
            400.0,
            np.nan,
        )
        np.testing.assert_allclose([np.log(m.pgv), m.ln_std_pgv], expected, rtol=1e-13)


# Vectorized scenarios
######################

KEYS = ["pga", "ln_pga", "ln_std_pga", "pgv", "ln_std_pgv", "spec_accels", "ln_stds"]


def scalar_results(row, **options):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = Zea06(Scenario(**dict(zip(NAMES, row))), **options)
    return {key: getattr(m, key) for key in KEYS}


@pytest.mark.parametrize("options", [OPTIONS[0], OPTIONS[2]])
def test_vectorized_matches_scalar(grid, options):
    rows, _ = grid
    rows = rows[::7]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = Zea06(grid_scenario(rows), **options)
    expected = [scalar_results(row, **options) for row in rows]
    for key in KEYS:
        desired = np.array([e[key] for e in expected])
        actual = getattr(m, key)
        assert actual.shape == desired.shape, key
        np.testing.assert_array_equal(actual, desired, err_msg=key)


def test_scalar_shapes():
    for event_type in EVENT_TYPES:
        m = Zea06(
            Scenario(
                mag=7.0,
                dist_rup=80.0,
                depth_tor=50.0,
                v_s30=400.0,
                event_type=event_type,
            )
        )
        for key in ["pga", "ln_pga", "ln_std_pga", "pgv", "ln_std_pgv"]:
            assert isinstance(getattr(m, key), float), key
        assert m.spec_accels.shape == (21,)
        assert m.ln_stds.shape == (21,)


def test_mixed_event_types():
    event_type = np.array(["interface", "intraslab", "intraslab", "interface"])
    mag = np.array([9.0, 7.0, 6.0, 8.0])
    depth_tor = np.array([np.nan, 50.0, 120.0, np.nan])
    depth_2_5 = np.array([7.0, np.nan, 4.0, np.nan])
    m = Zea06(
        Scenario(
            mag=mag,
            dist_rup=100.0,
            depth_tor=depth_tor,
            depth_2_5=depth_2_5,
            v_s30=360.0,
            event_type=event_type,
        ),
        basin=True,
        m9=True,
    )
    for i in range(4):
        row = (mag[i], 100.0, depth_tor[i], 360.0, depth_2_5[i], event_type[i])
        expected = scalar_results(row, basin=True, m9=True)
        for key in KEYS:
            np.testing.assert_array_equal(
                getattr(m, key)[i], expected[key], err_msg=key
            )
    assert np.all(np.isfinite(m.spec_accels))


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [6.0, 7.0, 8.0], [20.0, 50.0, 100.0, 300.0], indexing="ij"
    )
    m = Zea06(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            depth_tor=60.0,
            v_s30=400.0,
            event_type="intraslab",
        )
    )
    assert m.pga.shape == (3, 4)
    assert m.pgv.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    expected = scalar_results((8.0, 50.0, 60.0, 400.0, None, "intraslab"))
    assert m.pga[2, 1] == expected["pga"]
    np.testing.assert_array_equal(m.spec_accels[2, 1], expected["spec_accels"])


def test_interface_depth_can_be_none():
    params = dict(
        mag=np.array([7.0, 8.5]), dist_rup=np.array([30.0, 90.0]), v_s30=500.0
    )
    m = Zea06(Scenario(event_type="interface", **params))
    expected = Zea06(
        Scenario(event_type="interface", depth_tor=np.array([10.0, 500.0]), **params)
    )
    np.testing.assert_array_equal(m.spec_accels, expected.spec_accels)


@pytest.mark.parametrize(
    "event_type", ["intraslab", np.array(["interface", "intraslab"])]
)
def test_missing_depth_raises(event_type):
    with pytest.raises(ValueError, match="depth_tor is required"):
        Zea06(Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0, event_type=event_type))


@pytest.mark.filterwarnings("ignore:Using default value")
def test_invalid_event_types():
    with pytest.warns(UserWarning, match="1 of 3 values"):
        m = Zea06(
            Scenario(
                mag=7.0,
                dist_rup=50.0,
                depth_tor=60.0,
                v_s30=400.0,
                event_type=np.array(["interface", "crustal", "intraslab"]),
            )
        )
    assert np.all(np.isnan(m.spec_accels[1])) and np.isnan(m.pgv[1])
    assert np.all(np.isfinite(m.spec_accels[[0, 2]]))


def vector_scenario(n=500, seed=0):
    rng = np.random.default_rng(seed)
    return Scenario(
        mag=rng.uniform(5.0, 9.5, n),
        dist_rup=rng.uniform(0.0, 700.0, n),
        depth_tor=rng.uniform(10.0, 150.0, n),
        depth_2_5=np.where(rng.uniform(size=n) < 0.2, np.nan, rng.uniform(0, 9, n)),
        v_s30=rng.uniform(150.0, 1000.0, n),
        event_type=rng.choice(EVENT_TYPES, n),
    )


@pytest.mark.parametrize(
    "ims",
    [
        ["pga"],
        "pga",
        ["pgv"],
        ["psa_all"],
        ["psa_10p000"],
        ["pga", "psa_0p750", "psa_7p500", "pgv"],
    ],
)
@pytest.mark.parametrize("options", OPTIONS)
def test_ims_match_default(ims, options):
    s = vector_scenario()
    full = Zea06(s, **options)
    m = Zea06(s, ims=ims, **options)
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
    if computed in (["pga"], ["pgv"], ["psa_10p000"]):
        assert m._ln_resp.shape == (500, 1)
        assert m._ln_std.shape == (500, 1)


def test_ln_pga():
    m = Zea06(vector_scenario(), ims=["pga"])
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
