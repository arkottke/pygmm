"""Test the NGA-East seed model logic tree of USGS nshmp-lib."""

import itertools
import math
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from pygmm.ground_motion import nga_east
from pygmm.ground_motion.nga_east import TABLE_DISTS, TABLE_MAGS, NgaEast, NgaEastSeeds
from pygmm.model import Scenario

from .test_nga_east import REF_CASES, ref_cpa_ratio, ref_sigma, ref_site_amp, ref_table

PERIODS = [0.05, 0.3, 1.0, 2.5]

DATA = os.path.join(os.path.dirname(__file__), "data")
# Results of nshmp-lib for the 432 inputs of nga-east-inputs.csv and 8 intensity
# measures, given with 10 decimals
NSHMP = pd.read_csv(os.path.join(DATA, "nga_east_seeds-nshmp.csv.gz"), comment="#")

VARIANTS = list(NgaEastSeeds.GMM_IDS.items())


def compute(version="2026", adjusted=False, cpa=False, ims=None, **kwds):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return NgaEastSeeds(
            Scenario(**kwds), version=version, adjusted=adjusted, cpa=cpa, ims=ims
        )


def full_ln_resp(m):
    """ln response at all of the NgaEast periods (NaN if not computed)."""
    out = np.full(np.shape(m._ln_resp)[:-1] + (len(NgaEast.PERIODS),), np.nan)
    out[..., m._indices] = m._ln_resp
    return out


def full_ln_std(m):
    """ln standard deviation at all of the NgaEast periods."""
    out = np.full(np.shape(m._ln_std)[:-1] + (len(NgaEast.PERIODS),), np.nan)
    out[..., m._indices] = m._ln_std
    return out


def nshmp_columns(imts):
    period = [0.0 if s == "PGA" else float(s[2:].replace("P", ".")) for s in imts]
    cols = np.searchsorted(NgaEast.PERIODS, period)
    np.testing.assert_array_equal(NgaEast.PERIODS[cols], period)
    return cols


def nshmp_scenario(df):
    kwds = dict(
        mag=df["mag"].values,
        dist_jb=df["dist_jb"].values,
        dist_rup=df["dist_rup"].values,
        v_s30=df["v_s30"].values,
    )
    if "depth_sed" in df:
        kwds["depth_sed"] = df["depth_sed"].values
    return kwds


def test_nshmp_inputs():
    assert len(NSHMP) == 432 * 8
    assert NSHMP.loc[NSHMP.depth_sed.notna(), "index"].nunique() == 11
    # The reference inputs have dist_jb equal to dist_rup
    np.testing.assert_array_equal(NSHMP.dist_jb, NSHMP.dist_rup)


@pytest.mark.parametrize("gmm,options", VARIANTS)
def test_nshmp(gmm, options):
    m = compute(**options, **nshmp_scenario(NSHMP))
    rows = np.arange(len(NSHMP))
    cols = nshmp_columns(NSHMP.imt)
    median = np.exp(full_ln_resp(m)[rows, cols])
    sigma = full_ln_std(m)[rows, cols]
    # The reference values are rounded to 10 decimals
    np.testing.assert_allclose(median, NSHMP["median_" + gmm], rtol=0, atol=5.1e-11)
    np.testing.assert_allclose(sigma, NSHMP["sigma"], rtol=0, atol=5.1e-11)


def test_nshmp_variants_differ():
    on_cpa = NSHMP.depth_sed.notna()
    for a, b in [
        ("2018", "2023"),
        ("2023", "2026"),
        ("2023", "2023_ADJUSTED"),
        ("2026", "2026_ADJUSTED"),
    ]:
        diff = NSHMP[f"median_NGA_EAST_SEEDS_{a}"] - NSHMP[f"median_NGA_EAST_SEEDS_{b}"]
        assert np.max(np.abs(diff)) > 1e-3
    for v in ["2023", "2026"]:
        diff = (
            NSHMP[f"median_NGA_EAST_SEEDS_{v}_CPA"]
            - NSHMP[f"median_NGA_EAST_SEEDS_{v}"]
        )
        assert np.all(diff[~on_cpa] == 0)
        assert np.all(np.abs(diff[on_cpa & (NSHMP.v_s30 < 3000)]) > 0)


def test_seed_weights():
    weights = NgaEastSeeds.SEED_WEIGHTS
    assert len(weights) == 14
    assert math.isclose(sum(weights.values()), 1.0, rel_tol=1e-14)
    assert weights["SP16"] == 0.1089
    # The table-based seeds and SP16
    assert set(weights) == set(nga_east.SEED_TABLE_IDS) | {"SP16"}
    assert len(nga_east.SEED_TABLE_IDS) == 13
    assert set(nga_east.SEEDS_WITHOUT_PGV) <= set(nga_east.SEED_TABLE_IDS)


# Scalar transliteration of nshmp-lib NgaEast.UsgsSeeds_2018, UsgsSeeds_2023,
# NgaEast.Seed, and ShahjoueiPezeshk_2016, used to check the vectorized model at
# cases that are not in the nshmp-lib results


def ref_pezeshk_mean(c, mag, dist):
    r = math.sqrt(dist * dist + c["c11"] * c["c11"])
    mu = (
        c["c1"]
        + (c["c2"] * mag)
        + (c["c3"] * mag * mag)
        + (c["c4"] + c["c5"] * mag) * min(math.log10(r), math.log10(60.0))
        + (c["c6"] + c["c7"] * mag)
        * max(min(math.log10(r / 60.0), math.log10(2.0)), 0.0)
        + (c["c8"] + c["c9"] * mag) * max(math.log10(r / 120.0), 0.0)
        + (c["c10"] * r)
    )
    return mu * math.log(10.0)


def ref_row(coeff, j):
    return {name: coeff[name][j] for name in coeff.dtype.names}


def ref_seed_rock(seed, j, mag, dist_rup, dist_jb):
    """ln hard-rock median at period j and ln hard-rock PGA of a seed."""
    j_pga = NgaEast.INDEX_PGA
    if seed in nga_east.SEED_TABLE_IDS:
        n_r, n_m = len(TABLE_DISTS), len(TABLE_MAGS)
        tables = nga_east._load_seed_tables()[nga_east.SEED_TABLE_IDS.index(seed)]
        mu = ref_table(tables[:, j].reshape(n_r, n_m), mag, dist_rup)
        mu_pga = ref_table(tables[:, j_pga].reshape(n_r, n_m), mag, dist_rup)
        return mu, mu_pga
    assert seed == "SP16"
    coeff, dist = NgaEastSeeds.COEFF_SP16, dist_jb
    return (
        ref_pezeshk_mean(ref_row(coeff, j), mag, dist),
        ref_pezeshk_mean(ref_row(coeff, j_pga), mag, dist),
    )


def ref_apply(mu, f_s, s_s, z_score):
    """nshmp-lib NgaEast.SiteTerm.apply."""
    if f_s == 0.0:
        return mu
    mu_amp = mu + f_s
    return math.log(
        0.185 * math.exp(mu_amp + s_s * z_score)
        + 0.63 * math.exp(mu_amp)
        + 0.185 * math.exp(mu_amp - s_s * z_score)
    )


def ref_seeds(version, adjusted, cpa, mag, dist_rup, dist_jb, v_s30, depth_sed):
    """ln median and sigma at the periods of the nshmp-lib seed tree."""
    z_score = 1.0 if version == "2018" else 1.645
    if depth_sed is None or math.isnan(depth_sed):
        z_scale = 0.0
    else:
        z_scale = (1.0 - math.exp(-depth_sed / 0.2)) ** 4
    ln_resps, ln_stds = [], []
    coeff = NgaEast.COEFF
    for j, period in enumerate(NgaEast.PERIODS):
        if j == NgaEast.INDEX_PGV:
            continue
        c = ref_row(coeff, j)
        nga_adj = c["nga_adj"]
        if v_s30 > 1000.0:
            nga_adj += c["vs30_b"] * math.log(min(v_s30, 2000.0) / 1000.0)
        mu_adj = (1.0 - z_scale) * nga_adj if adjusted else 0.0
        apply_cpa = cpa and z_scale > 0.0
        f_cpa = (
            math.log(ref_cpa_ratio(j, depth_sed, mag, dist_jb)) if apply_cpa else 0.0
        )
        mus, wts = [], []
        for seed, weight in NgaEastSeeds.SEED_WEIGHTS.items():
            mu, mu_pga = ref_seed_rock(seed, j, mag, dist_rup, dist_jb)
            mu += mu_adj
            pga = math.exp(mu_pga + mu_adj)
            f_s, s_s = ref_site_amp(c, version, period, pga, v_s30)
            mu_site = ref_apply(mu, f_s, s_s, z_score)
            if apply_cpa:
                f_ref, _ = ref_site_amp(c, version, period, pga, 1000.0)
                mu_site += f_cpa - f_ref * z_scale
            mus.append(mu_site)
            wts.append(weight)
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


# Includes SP16 at short and zero Joyner-Boore distances, magnitudes beyond the
# SP16 magnitude break, and the cases of the NgaEast tests
SEED_CASES = REF_CASES + [
    (4.3, 3.0, 0.0, 760.0, np.nan),
    (5.0, 1.0, 0.5, 300.0, 0.0),
    (6.5, 12.0, 2.0, 1100.0, 0.25),
    (6.6, 140.0, 139.0, 400.0, np.nan),
]


@pytest.mark.parametrize("gmm,options", VARIANTS)
@pytest.mark.parametrize("mag,dist_rup,dist_jb,v_s30,depth_sed", SEED_CASES)
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
    ln_resp, ln_std = ref_seeds(
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


def test_sp16_uses_dist_jb():
    kwds = dict(mag=6.0, dist_rup=10.0, v_s30=760.0)
    # The tree depends on dist_jb through SP16
    a = compute(dist_jb=10.0, **kwds)
    b = compute(dist_jb=2.0, **kwds)
    assert np.all(b._ln_resp > a._ln_resp)
    # Only the SP16 branches
    branches_a = a.ln_branches("psa")
    branches_b = b.ln_branches("psa")
    seeds = [s for s in NgaEastSeeds.SEED_WEIGHTS for _ in range(2)]
    for seed, (_, mu_a, _), (_, mu_b, _) in zip(seeds, branches_a, branches_b):
        if seed == "SP16":
            assert np.all(mu_b > mu_a)
        else:
            np.testing.assert_array_equal(mu_a, mu_b)


def test_dist_jb_is_required():
    kwds = dict(mag=6.0, dist_rup=10.0, v_s30=760.0)
    with pytest.raises(ValueError, match="dist_jb is required"):
        compute(**kwds)
    with pytest.raises(ValueError, match="dist_jb is required"):
        compute(mag=np.array([5.0, 6.0]), dist_rup=10.0, v_s30=760.0)


def test_graizer_first_distance():
    # The Graizer16 and Graizer17 tables give the first distance as 0.01 km,
    # which nshmp-lib uses as 1e-5 km like the other tables
    kwds = dict(mag=6.0, dist_jb=0.0, v_s30=3000.0)
    a = compute(dist_rup=0.0, **kwds)
    b = compute(dist_rup=1e-5, **kwds)
    np.testing.assert_array_equal(a._ln_resp, b._ln_resp)
    for branch_a, branch_b in zip(a.ln_branches("psa"), b.ln_branches("psa")):
        np.testing.assert_array_equal(branch_a[1], branch_b[1])


def test_seed_tables():
    tables = nga_east._load_seed_tables()
    n = len(TABLE_DISTS) * len(TABLE_MAGS)
    assert tables.shape == (13, n, 23)
    for k, seed in enumerate(nga_east.SEED_TABLE_IDS):
        pgv = tables[k, :, NgaEast.INDEX_PGV]
        if seed in nga_east.SEEDS_WITHOUT_PGV:
            assert np.all(np.isnan(pgv))
        else:
            assert np.all(np.isfinite(pgv))
        assert np.all(np.isfinite(tables[k, :, 1:]))
    # Table node: B_ab95 PGA at 5 km and M 4.0 is 0.3168 g
    k = nga_east.SEED_TABLE_IDS.index("B_ab95")
    assert math.isclose(np.exp(tables[k, 2 * 11, NgaEast.INDEX_PGA]), 0.3168)


def test_sigma_is_that_of_nga_east():
    kwds = dict(mag=6.0, dist_rup=50.0, dist_jb=48.0, v_s30=1300.0, depth_sed=0.3)
    expected = NgaEast(Scenario(**kwds), ims=["pga", "psa_all"])._ln_std
    for _, options in VARIANTS:
        np.testing.assert_allclose(compute(**options, **kwds)._ln_std, expected, 1e-15)


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
        NgaEastSeeds(
            Scenario(mag=6.0, dist_rup=20.0, dist_jb=20.0, v_s30=760.0), **options
        )


def test_version_is_a_string():
    s = Scenario(mag=6.0, dist_rup=20.0, dist_jb=20.0, v_s30=760.0)
    assert NgaEastSeeds(s, version=2023).version == "2023"


def test_pgv():
    s = Scenario(mag=6.0, dist_rup=20.0, dist_jb=20.0, v_s30=760.0)
    with pytest.raises(ValueError, match="does not provide 'pgv'"):
        NgaEastSeeds(s, ims=["pgv"])
    with pytest.raises(NotImplementedError):
        NgaEastSeeds(s).pgv
    assert NgaEastSeeds.INDEX_PGV is None


def test_out_of_range_values_warn():
    with pytest.warns(UserWarning, match="v_s30"):
        NgaEastSeeds(
            Scenario(
                mag=6.0, dist_rup=20.0, dist_jb=20.0, v_s30=np.array([100.0, 760.0])
            )
        )
    with pytest.warns(UserWarning, match="mag"):
        NgaEastSeeds(Scenario(mag=8.5, dist_rup=20.0, dist_jb=20.0, v_s30=760.0))


# Vectorized models


MAGS = [4.0, 6.75, 8.5]
DISTS = [0.0, 75.0, 1500.0]
V_S30S = [140.0, 260.0, 1300.0, 3000.0]
DEPTHS = [np.nan, 0.3, 5.0]
GRID_ROWS = list(itertools.product(MAGS, DISTS, V_S30S, DEPTHS))
KEYS = [
    "pga",
    "ln_pga",
    "ln_std_pga",
    "pgv",
    "ln_std_pgv",
    "spec_accels",
    "ln_stds",
    "interp_spec_accels",
    "interp_ln_stds",
]
GRID_MODELS = [("tree", gmm) for gmm, _ in VARIANTS]


def make(kind, name, **kwds):
    assert kind == "tree"
    return compute(**NgaEastSeeds.GMM_IDS[name], **kwds)


def results(m):
    out = {}
    for key in KEYS:
        try:
            if key.startswith("interp"):
                out[key] = getattr(m, key)(PERIODS)
            else:
                out[key] = getattr(m, key)
        except NotImplementedError:
            out[key] = None
    return out


@pytest.fixture(scope="module", params=GRID_MODELS, ids=lambda p: p[1])
def grid(request):
    kind, name = request.param
    mag, dist_rup, v_s30, depth_sed = (np.array(c) for c in zip(*GRID_ROWS))
    m = make(
        kind,
        name,
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=0.9 * dist_rup,
        v_s30=v_s30,
        depth_sed=depth_sed,
    )
    expected = [
        results(
            make(kind, name, mag=mg, dist_rup=d, dist_jb=0.9 * d, v_s30=v, depth_sed=z)
        )
        for mg, d, v, z in GRID_ROWS
    ]
    return m, expected


@pytest.mark.parametrize("key", KEYS)
def test_vectorized_matches_scalar(grid, key):
    m, expected = grid
    actual = results(m)[key]
    if actual is None:
        assert all(e[key] is None for e in expected)
        return
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
    assert m.ln_std_pga.shape == (n,)
    assert m.spec_accels.shape == (n, 21)
    assert m.ln_stds.shape == (n, 21)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))


@pytest.mark.parametrize("kind,name", GRID_MODELS)
def test_scalar_shapes(kind, name):
    m = make(
        kind, name, mag=6.5, dist_rup=20.0, dist_jb=19.0, v_s30=760.0, depth_sed=2.0
    )
    assert isinstance(m.pga, float)
    assert isinstance(m.ln_pga, float)
    assert isinstance(m.ln_std_pga, float)
    assert m.spec_accels.shape == (21,)
    assert m.ln_stds.shape == (21,)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)


@pytest.mark.parametrize("kind,name", GRID_MODELS)
def test_scalars_broadcast_with_arrays(kind, name):
    mag = np.array([5.5, 6.5, 7.5])
    dist_rup = np.array([5.0, 20.0, 80.0])

    def scalar(**kwds):
        return make(kind, name, **kwds).spec_accels

    m = make(
        kind,
        name,
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=0.9 * dist_rup,
        v_s30=400.0,
        depth_sed=0.3,
    )
    expected = [
        scalar(mag=mg, dist_rup=d, dist_jb=0.9 * d, v_s30=400.0, depth_sed=0.3)
        for mg, d in zip(mag, dist_rup)
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)

    # Only dist_jb is an array
    dist_jb = np.array([0.0, 3.0, 10.0])
    m = make(kind, name, mag=6.0, dist_rup=10.0, dist_jb=dist_jb, v_s30=400.0)
    assert m.spec_accels.shape == (3, 21)
    assert m.ln_stds.shape == (3, 21)
    np.testing.assert_array_equal(
        m.spec_accels,
        [scalar(mag=6.0, dist_rup=10.0, dist_jb=d, v_s30=400.0) for d in dist_jb],
    )

    # Only V_S30 is an array
    v_s30 = np.array([200.0, 760.0, 2500.0])
    m = make(kind, name, mag=6.0, dist_rup=20.0, dist_jb=18.0, v_s30=v_s30)
    assert m.spec_accels.shape == (3, 21)
    np.testing.assert_array_equal(
        m.spec_accels,
        [scalar(mag=6.0, dist_rup=20.0, dist_jb=18.0, v_s30=v) for v in v_s30],
    )

    # Only the magnitude is an array, the standard deviation has the full shape
    m = make(kind, name, mag=mag, dist_rup=20.0, dist_jb=18.0, v_s30=760.0)
    assert m.ln_stds.shape == (3, 21)


@pytest.mark.parametrize("kind,name", GRID_MODELS)
def test_multidimensional_inputs(kind, name):
    mag, depth_sed = np.meshgrid(
        [5.5, 6.5, 7.5], [np.nan, 0.1, 0.5, 2.0], indexing="ij"
    )
    m = make(
        kind,
        name,
        mag=mag,
        dist_rup=20.0,
        dist_jb=18.0,
        v_s30=400.0,
        depth_sed=depth_sed,
    )
    assert m.pga.shape == (3, 4)
    assert m.spec_accels.shape == (3, 4, 21)
    assert m.ln_stds.shape == (3, 4, 21)
    expected = make(
        kind, name, mag=7.5, dist_rup=20.0, dist_jb=18.0, v_s30=400.0, depth_sed=0.1
    ).pga
    assert m.pga[2, 1] == expected


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


def make_vector(kind, name, ims=None):
    assert kind == "tree"
    return NgaEastSeeds(vector_scenario(), ims=ims, **NgaEastSeeds.GMM_IDS[name])


@pytest.mark.parametrize("kind,name", GRID_MODELS)
@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
def test_pga_only_matches_default(kind, name, ims):
    full = make_vector(kind, name)
    pga_only = make_vector(kind, name, ims=ims)
    np.testing.assert_array_equal(pga_only.pga, full.pga)
    np.testing.assert_array_equal(pga_only.ln_std_pga, full.ln_std_pga)
    # Only one period is computed
    assert pga_only._ln_resp.shape == (500, 1)
    assert pga_only._ln_std.shape == (500, 1)


@pytest.mark.parametrize("kind,name", GRID_MODELS)
@pytest.mark.parametrize("ims", [["psa_all"], ["psa_ngawest2_21"], ["pga", "psa_all"]])
def test_psa_matches_default(kind, name, ims):
    full = make_vector(kind, name)
    m = make_vector(kind, name, ims=ims)
    np.testing.assert_array_equal(m.periods, full.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels)
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds)
    np.testing.assert_allclose(
        m.interp_spec_accels(PERIODS), full.interp_spec_accels(PERIODS), rtol=1e-12
    )


@pytest.mark.parametrize("kind,name", GRID_MODELS)
def test_selected_periods_match_default(kind, name):
    full = make_vector(kind, name)
    m = make_vector(kind, name, ims=["psa_0p200", "psa_1p000"])
    assert m.psa_ims == ["psa_0p200", "psa_1p000"]
    cols = np.searchsorted(full.periods, [0.2, 1.0])
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])
    with pytest.raises(ValueError, match="pga was not computed"):
        m.pga


@pytest.mark.parametrize(
    "ims,match",
    [
        (["pgd"], "does not provide 'pgd'"),
        (["pgv"], "does not provide 'pgv'"),
        (["psa_0p040"], "does not provide 'psa_0p040'"),
        (["sa"], "not a valid intensity measure"),
    ],
)
def test_invalid_ims_raise(ims, match):
    with pytest.raises(ValueError, match=match):
        NgaEastSeeds(vector_scenario(), ims=ims)


@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(ims):
    s = vector_scenario()
    m = NgaEastSeeds(s, cpa=True, ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_array_equal(m.ln_pga, NgaEastSeeds(s, cpa=True).ln_pga)
