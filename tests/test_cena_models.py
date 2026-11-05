"""Test the vectorized CENA models against scalar scenarios.

The models are Atkinson and Boore (2006), Campbell (2003), Pezeshk et al.
(2011), and Tavakoli and Pezeshk (2005).
"""

import itertools
import warnings

import numpy as np
import pytest

import pygmm

# Magnitudes cover the standard deviation breaks at 7.0 (Pea11), 7.16 (C03),
# and 7.2 (TP05)
MAGS = [4.5, 5.5, 6.5, 7.0, 7.1, 7.16, 7.2, 7.5, 8.2]
# Distances cover the segments at 10, 70, 130, and 140 km
DISTS = [1.0, 5.0, 10.0, 50.0, 70.0, 100.0, 130.0, 135.0, 140.0, 150.0, 500.0]
# v_s30 covers the nonlinear site segments at 180, 300, and 760 m/s in AB06.
# Zero and 2000 m/s and greater use the hard-rock coefficients.
V_S30S = [0.0, 150.0, 180.0, 250.0, 300.0, 500.0, 760.0, 1999.0, 2000.0, 2500.0]
PERIODS = [0.05, 0.3, 1.0, 2.5]

# Model, whether it uses v_s30, and a single spectral period
MODELS = [
    (pygmm.AtkinsonBoore2006, True, "psa_1p000"),
    (pygmm.Campbell2003, False, "psa_1p000"),
    (pygmm.PezeshkZandiehTavakoli2011, False, "psa_1p000"),
    (pygmm.TavakoliPezeshk05, False, "psa_1p000"),
]
MODEL_IDS = [m[0].ABBREV for m in MODELS]

KEYS = [
    "pga",
    "ln_std_pga",
    "pgv",
    "pgd",
    "spec_accels",
    "ln_stds",
    "interp_spec_accels",
    "interp_ln_stds",
]


@pytest.fixture(autouse=True)
def ignore_limit_warnings():
    # Some magnitudes are outside of the recommended limits
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        yield


def make_scenario(uses_v_s30, mag, dist_rup, v_s30=760.0):
    if uses_v_s30:
        return pygmm.Scenario(mag=mag, dist_rup=dist_rup, v_s30=v_s30)
    return pygmm.Scenario(mag=mag, dist_rup=dist_rup)


def get_value(m, key):
    if key.startswith("interp"):
        return getattr(m, key)(PERIODS)
    try:
        return getattr(m, key)
    except NotImplementedError:
        return None


def grid_rows(uses_v_s30):
    v_s30s = V_S30S if uses_v_s30 else [760.0]
    dists = DISTS if uses_v_s30 else [0.0] + DISTS
    return list(itertools.product(MAGS, dists, v_s30s))


@pytest.mark.parametrize("key", KEYS)
@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_vectorized_matches_scalar(cls, uses_v_s30, im, key):
    rows = grid_rows(uses_v_s30)
    mag, dist_rup, v_s30 = (np.array(c) for c in zip(*rows))
    actual = get_value(cls(make_scenario(uses_v_s30, mag, dist_rup, v_s30)), key)
    expected = [get_value(cls(make_scenario(uses_v_s30, *row)), key) for row in rows]
    if expected[0] is None:
        # Not provided by the model
        assert actual is None
        return
    expected = np.array(expected)
    assert actual.shape == expected.shape
    if key.startswith("interp"):
        # interp1d evaluates 2-D arrays with slightly different floating point round-off
        np.testing.assert_allclose(actual, expected, rtol=1e-12)
    else:
        np.testing.assert_array_equal(actual, expected)


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_vectorized_shapes(cls, uses_v_s30, im):
    rows = grid_rows(uses_v_s30)
    mag, dist_rup, v_s30 = (np.array(c) for c in zip(*rows))
    m = cls(make_scenario(uses_v_s30, mag, dist_rup, v_s30))
    n = len(rows)
    n_periods = len(m.periods)
    assert m.spec_accels.shape == (n, n_periods)
    assert m.ln_stds.shape == (n, n_periods)
    assert m.interp_spec_accels(PERIODS).shape == (n, len(PERIODS))
    assert m.interp_ln_stds(PERIODS).shape == (n, len(PERIODS))
    if cls.INDEX_PGA is not None:
        assert m.pga.shape == (n,)
        assert m.ln_std_pga.shape == (n,)


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_scalar_shapes_are_unchanged(cls, uses_v_s30, im):
    m = cls(make_scenario(uses_v_s30, 6.5, 20.0))
    assert m.spec_accels.shape == (len(m.periods),)
    assert m.ln_stds.shape == (len(m.periods),)
    assert m.interp_spec_accels(PERIODS).shape == (len(PERIODS),)
    if cls.INDEX_PGA is not None:
        assert isinstance(m.pga, float)
        assert isinstance(m.ln_std_pga, float)


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_scalars_broadcast_with_arrays(cls, uses_v_s30, im):
    # A scalar distance (and v_s30) with an array of magnitudes
    mag = np.array([5.5, 6.5, 7.5])
    m = cls(make_scenario(uses_v_s30, mag, 80.0, 400.0))
    expected = [
        cls(make_scenario(uses_v_s30, mg, 80.0, 400.0)).spec_accels for mg in mag
    ]
    np.testing.assert_array_equal(m.spec_accels, expected)


def test_v_s30_array_broadcasts_with_scalars():
    v_s30 = np.array([0.0, 200.0, 400.0, 1000.0])
    m = pygmm.AtkinsonBoore2006(pygmm.Scenario(mag=6.5, dist_rup=30.0, v_s30=v_s30))
    expected = [
        pygmm.AtkinsonBoore2006(pygmm.Scenario(mag=6.5, dist_rup=30.0, v_s30=v)).pga
        for v in v_s30
    ]
    np.testing.assert_array_equal(m.pga, expected)


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_multidimensional_inputs(cls, uses_v_s30, im):
    # A grid of magnitudes and distances keeps its shape
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [5.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    m = cls(make_scenario(uses_v_s30, mag, dist_rup, 400.0))
    assert m.spec_accels.shape == (3, 4, len(m.periods))
    assert m.ln_stds.shape == (3, 4, len(m.periods))
    np.testing.assert_array_equal(
        m.spec_accels[2, 1],
        cls(make_scenario(uses_v_s30, 7.5, 20.0, 400.0)).spec_accels,
    )


def test_ab06_multidimensional_v_s30():
    # The site term interpolates over periods for each scenario
    mag, v_s30 = np.meshgrid([5.5, 7.5], [0.0, 200.0, 400.0], indexing="ij")
    m = pygmm.AtkinsonBoore2006(pygmm.Scenario(mag=mag, dist_rup=50.0, v_s30=v_s30))
    assert m.pga.shape == (2, 3)
    assert m.pgv.shape == (2, 3)
    expected = pygmm.AtkinsonBoore2006(
        pygmm.Scenario(mag=7.5, dist_rup=50.0, v_s30=200.0)
    )
    np.testing.assert_array_equal(m.spec_accels[1, 1], expected.spec_accels)
    assert m.pga[1, 1] == expected.pga


def vector_scenario(uses_v_s30, n=500, seed=0):
    rng = np.random.default_rng(seed)
    return make_scenario(
        uses_v_s30,
        rng.uniform(5.0, 8.0, n),
        rng.uniform(1.0, 300.0, n),
        rng.uniform(150.0, 1500.0, n),
    )


@pytest.mark.parametrize("ims", [["pga"], "pga", ("pga",)])
@pytest.mark.parametrize(
    "cls,uses_v_s30,im",
    [m for m in MODELS if m[0].INDEX_PGA is not None],
    ids=[m[0].ABBREV for m in MODELS if m[0].INDEX_PGA is not None],
)
def test_pga_only_matches_default(cls, uses_v_s30, im, ims):
    s = vector_scenario(uses_v_s30)
    full = cls(s)
    pga_only = cls(s, ims=ims)
    np.testing.assert_array_equal(pga_only.pga, full.pga)
    np.testing.assert_array_equal(pga_only.ln_pga, full.ln_pga)
    np.testing.assert_array_equal(pga_only.ln_std_pga, full.ln_std_pga)
    # Only one period is computed
    assert pga_only._ln_resp.shape == (500, 1)
    assert pga_only._ln_std.shape == (500, 1)


@pytest.mark.parametrize(
    "cls,uses_v_s30,im",
    [m for m in MODELS if m[0].INDEX_PGA is not None],
    ids=[m[0].ABBREV for m in MODELS if m[0].INDEX_PGA is not None],
)
def test_pga_only_scalar_matches_default(cls, uses_v_s30, im):
    s = make_scenario(uses_v_s30, 6.5, 20.0, 400.0)
    full = cls(s)
    pga_only = cls(s, ims=["pga"])
    assert pga_only.pga == full.pga
    assert pga_only.ln_std_pga == full.ln_std_pga
    assert isinstance(pga_only.pga, float)


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_single_period_matches_default(cls, uses_v_s30, im):
    s = vector_scenario(uses_v_s30)
    full = cls(s)
    m = cls(s, ims=[im])
    assert m.psa_ims == [im]
    index = full.psa_ims.index(im)
    assert m.spec_accels.shape == (500, 1)
    np.testing.assert_array_equal(m.spec_accels[:, 0], full.spec_accels[:, index])
    np.testing.assert_array_equal(m.ln_stds[:, 0], full.ln_stds[:, index])


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_psa_all_matches_default(cls, uses_v_s30, im):
    s = vector_scenario(uses_v_s30)
    full = cls(s)
    m = cls(s, ims=["psa_all"])
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels)
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds)
    np.testing.assert_allclose(
        m.interp_spec_accels(PERIODS), full.interp_spec_accels(PERIODS), rtol=1e-12
    )


@pytest.mark.parametrize("im", ["pgv", "pgd"])
def test_ab06_pgv_pgd_only_matches_default(im):
    s = vector_scenario(True)
    full = pygmm.AtkinsonBoore2006(s)
    m = pygmm.AtkinsonBoore2006(s, ims=[im])
    np.testing.assert_array_equal(getattr(m, im), getattr(full, im))
    with pytest.raises(ValueError, match="include 'pga' in ims"):
        m.pga


def test_c03_pga_not_provided():
    with pytest.raises(ValueError, match="does not provide 'pga'"):
        pygmm.Campbell2003(vector_scenario(False), ims=["pga"])


@pytest.mark.parametrize(
    "cls,uses_v_s30,im",
    [m for m in MODELS if m[0].INDEX_PGA is not None],
    ids=[m[0].ABBREV for m in MODELS if m[0].INDEX_PGA is not None],
)
@pytest.mark.parametrize("ims", [None, ["pga"]])
def test_ln_pga(cls, uses_v_s30, im, ims):
    s = vector_scenario(uses_v_s30)
    m = cls(s, ims=ims)
    assert m.ln_pga.shape == (500,)
    np.testing.assert_array_equal(np.exp(m.ln_pga), m.pga)
    np.testing.assert_array_equal(m.ln_pga, cls(s).ln_pga)
    np.testing.assert_allclose(m.ln_pga, np.log(m.pga), rtol=1e-14)


@pytest.mark.parametrize(
    "cls,uses_v_s30,im",
    [m for m in MODELS if m[0].INDEX_PGA is not None],
    ids=[m[0].ABBREV for m in MODELS if m[0].INDEX_PGA is not None],
)
def test_ln_pga_scalar(cls, uses_v_s30, im):
    m = cls(make_scenario(uses_v_s30, 6.5, 20.0, 400.0))
    assert isinstance(m.ln_pga, float)
    assert np.exp(m.ln_pga) == m.pga


@pytest.mark.parametrize("cls,uses_v_s30,im", MODELS, ids=MODEL_IDS)
def test_mismatched_lengths_raise(cls, uses_v_s30, im):
    with pytest.raises(ValueError):
        cls(
            make_scenario(
                uses_v_s30, np.array([6.0, 7.0]), np.array([10.0, 20.0, 30.0])
            )
        )


@pytest.mark.parametrize(
    "periods", [[-2.0, 0.0, 0.07, 1.1, 5.0, 6.0], [0.5], [0.025], [-1.0], [5.0]]
)
def test_interp_last_axis_matches_numpy(periods):
    from pygmm.ground_motion.atkinson_boore_2006 import interp_last_axis, interp_rows

    rng = np.random.default_rng(1)
    xp = np.array([-1.0, 0.0101, 0.025, 0.1, 1.0, 5.0])
    fp = rng.normal(size=(7, len(xp)))
    expected = [np.interp(periods, xp, row) for row in fp]
    np.testing.assert_array_equal(interp_last_axis(periods, xp, fp), expected)
    # Only the rows needed for the periods give the same values
    rows = interp_rows(periods, xp)
    np.testing.assert_array_equal(
        interp_last_axis(periods, xp[rows], fp[:, rows]), expected
    )
    np.testing.assert_array_equal(
        interp_last_axis(periods, xp[rows], fp[0, rows]), expected[0]
    )


@pytest.mark.parametrize("v_s30", [180.0, 300.0])
def test_ab06_site_term_is_continuous(v_s30):
    # The nonlinear site coefficient (Eq. 8 of Atkinson and Boore, 2006) is
    # continuous at the velocity breakpoints
    def calc(v):
        return pygmm.AtkinsonBoore2006(
            pygmm.Scenario(mag=7.0, dist_rup=10.0, v_s30=v)
        ).spec_accels

    np.testing.assert_allclose(calc(v_s30 - 1e-6), calc(v_s30 + 1e-6), rtol=1e-6)


def test_ab06_standard_deviation_in_natural_log_units():
    # Atkinson and Boore (2006) give a standard deviation of 0.30 in log10 units
    m = pygmm.AtkinsonBoore2006(pygmm.Scenario(mag=6.5, dist_rup=30.0, v_s30=760.0))
    np.testing.assert_allclose(m.ln_stds, np.log(10**0.30), rtol=1e-12)
    np.testing.assert_allclose(m.ln_std_pga, np.log(10**0.30), rtol=1e-12)


def test_ab06_pgv_and_pgd_units():
    # PGV (cm/sec) and PGD (cm) are not converted to gravity like PGA and PSA
    m = pygmm.AtkinsonBoore2006(pygmm.Scenario(mag=7.0, dist_rup=10.0, v_s30=760.0))
    log10_resp = m._ln_resp / np.log(10)
    np.testing.assert_allclose(m.pgv, 10 ** log10_resp[m.INDEX_PGV], rtol=1e-12)
    np.testing.assert_allclose(m.pgv, 44.34, rtol=1e-3)
    np.testing.assert_allclose(m.pgd, 19.05, rtol=1e-3)
    np.testing.assert_allclose(m.pga, 1.038, rtol=1e-3)


@pytest.mark.parametrize("v_s30", [2000.0, 2500.0])
def test_ab06_hard_rock_sites_use_rock_coefficients(v_s30):
    # Hard-rock sites use the hard-rock coefficients without site amplification,
    # the same as v_s30 of zero
    kwds = dict(mag=6.5, dist_rup=30.0)
    rock = pygmm.AtkinsonBoore2006(pygmm.Scenario(v_s30=0.0, **kwds))
    m = pygmm.AtkinsonBoore2006(pygmm.Scenario(v_s30=v_s30, **kwds))
    np.testing.assert_array_equal(m._ln_resp, rock._ln_resp)

    # Just below 2000 m/s, the B/C coefficients are used with site amplification
    soft = pygmm.AtkinsonBoore2006(pygmm.Scenario(v_s30=1999.0, **kwds))
    assert np.all(soft._ln_resp != m._ln_resp)


@pytest.mark.parametrize(
    "mag,period,c_12,c_13,c_14,sigma_reg",
    [
        # Coefficients for PGA and 1.0 s from OpenQuake (pezeshk_2011.py)
        (6.0, 0.0, -0.021050254, 0.37776584, 0.27905505, 0.020605025),
        (7.5, 0.0, -0.021050254, 0.37776584, 0.27905505, 0.020605025),
        (6.0, 1.0, -0.011795872, 0.35880759, 0.32487746, 0.022183642),
    ],
)
def test_pzt11_standard_deviation_in_natural_log_units(
    mag, period, c_12, c_13, c_14, sigma_reg
):
    # The standard deviations of Pezeshk et al. (2011) are in log10 units
    m = pygmm.PezeshkZandiehTavakoli2011(pygmm.Scenario(mag=mag, dist_rup=30.0))
    log10_std_mean = c_12 * mag + c_13 if mag <= 7.0 else -6.95e-3 * mag + c_14
    expected = np.log(10 ** np.sqrt(log10_std_mean**2 + sigma_reg**2))
    actual = m.ln_std_pga if period == 0 else m.interp_ln_stds([period])[0]
    # pygmm's coefficients are rounded to 4 significant digits
    np.testing.assert_allclose(actual, expected, rtol=1e-3)
