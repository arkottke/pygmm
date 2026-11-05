"""Test selecting intensity measures with the ``ims`` argument."""

import inspect
import warnings

import numpy as np
import pytest

import pygmm
from pygmm.model import GroundMotionModel

NGAWEST2_21 = [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4]
NGAWEST2_21 += [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]


def idriss(ims=None, n=300):
    rng = np.random.default_rng(0)
    s = pygmm.Scenario(
        mag=rng.uniform(5.0, 8.0, n),
        dist_rup=rng.uniform(0.0, 150.0, n),
        v_s30=rng.uniform(450.0, 1200.0, n),
        mechanism=rng.choice(["SS", "RS"], n),
    )
    return pygmm.Idriss2014(s, ims=ims)


def bssa14(ims=None, n=300):
    rng = np.random.default_rng(0)
    s = pygmm.Scenario(
        mag=rng.uniform(4.0, 8.0, n),
        dist_jb=rng.uniform(0.0, 200.0, n),
        v_s30=rng.uniform(180.0, 1300.0, n),
        mechanism=rng.choice(["U", "SS", "NS", "RS"], n),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pygmm.BooreStewartSeyhanAtkinson2014(s, ims=ims)


def nga_west2_scenario(n):
    rng = np.random.default_rng(0)
    dist_jb = rng.uniform(0.0, 200.0, n)
    on_hanging_wall = rng.choice([True, False], n)
    return pygmm.Scenario(
        mag=rng.uniform(4.0, 8.0, n),
        dist_jb=dist_jb,
        dist_rup=dist_jb + rng.uniform(0.0, 10.0, n),
        dist_x=np.where(on_hanging_wall, dist_jb, -dist_jb),
        on_hanging_wall=on_hanging_wall,
        dip=rng.choice([45.0, 90.0], n),
        v_s30=rng.uniform(180.0, 1300.0, n),
        mechanism=rng.choice(["SS", "NS", "RS"], n),
        region=rng.choice(["global", "california", "japan", "china"], n),
    )


def nga_west2_factory(model):
    def factory(ims=None, n=300):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return model(nga_west2_scenario(n), ims=ims)

    return factory


def nga_sub_final(model, **options):
    def factory(ims=None, n=300):
        rng = np.random.default_rng(0)
        s = pygmm.Scenario(
            mag=rng.uniform(5.0, 9.0, n),
            dist_rup=rng.uniform(10.0, 500.0, n),
            v_s30=rng.uniform(150.0, 1000.0, n),
            depth_tor=rng.uniform(5.0, 100.0, n),
            depth_2_5=rng.choice([np.nan, 1.0, 4.0], n),
            event_type=rng.choice(["interface", "intraslab"], n),
            region=rng.choice(["global", "alaska", "cascadia", "prvi"], n),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return model(s, ims=ims, **options)

    return factory


def nshmp_subduction(model, event_types, **options):
    def factory(ims=None, n=300):
        rng = np.random.default_rng(0)
        s = pygmm.Scenario(
            mag=rng.uniform(5.0, 9.5, n),
            dist_rup=rng.uniform(0.0, 500.0, n),
            v_s30=rng.uniform(150.0, 1000.0, n),
            depth_tor=rng.uniform(5.0, 150.0, n),
            depth_2_5=rng.choice([np.nan, 1.0, 4.0, 7.0], n),
            event_type=rng.choice(event_types, n),
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return model(s, ims=ims, **options)

    return factory


def nga_west2_nshmp(model, gmm):
    """Factory of the nshmp-lib variant (Gmm id) of an NGA-West2 model."""
    params = set(inspect.signature(model.__init__).parameters)
    options = {k: v for k, v in model.GMM_IDS[gmm].items() if k in params}
    values = {k: v for k, v in model.GMM_IDS[gmm].items() if k not in params}

    def factory(ims=None, n=300):
        rng = np.random.default_rng(0)
        dist_jb = rng.uniform(0.0, 200.0, n)
        depth_tor = rng.uniform(0.0, 10.0, n)
        s = pygmm.Scenario(
            mag=rng.uniform(5.0, 8.0, n),
            dist_jb=dist_jb,
            dist_rup=np.hypot(dist_jb, depth_tor),
            dist_x=rng.choice([-1.0, 1.0], n) * dist_jb,
            dip=rng.choice([45.0, 90.0], n),
            depth_tor=depth_tor,
            depth_bor=depth_tor + rng.uniform(5.0, 15.0, n),
            mechanism=rng.choice(["U", "SS", "NS", "RS"], n),
            v_s30=rng.uniform(180.0, 1300.0, n),
            depth_1_0=rng.choice([np.nan, 0.2, 0.6], n),
            depth_2_5=rng.choice([np.nan, 1.0, 4.0], n),
            **values,
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return model(s, ims=ims, **options)

    return factory


def nga_east_2026(ims=None, n=300):
    rng = np.random.default_rng(0)
    s = pygmm.Scenario(
        mag=rng.uniform(4.0, 8.2, n),
        dist_rup=rng.uniform(0.0, 1000.0, n),
        v_s30=rng.uniform(200.0, 3000.0, n),
    )
    return pygmm.NgaEast(s, version="2026", adjusted=True, ims=ims)


def nga_east_seeds_2026(ims=None, n=300):
    rng = np.random.default_rng(0)
    s = pygmm.Scenario(
        mag=rng.uniform(4.0, 8.2, n),
        dist_rup=rng.uniform(0.0, 1000.0, n),
        dist_jb=rng.uniform(0.0, 1000.0, n),
        v_s30=rng.uniform(200.0, 3000.0, n),
    )
    return pygmm.NgaEastSeeds(s, version="2026", adjusted=True, ims=ims)


MODELS = {
    "idriss": (idriss, 22),
    "bssa14": (bssa14, 105),
    "ask14": (nga_west2_factory(pygmm.AbrahamsonSilvaKamai2014), 22),
    "cb14": (nga_west2_factory(pygmm.CampbellBozorgnia2014), 21),
    "cy14": (nga_west2_factory(pygmm.ChiouYoungs2014), 24),
    "nga_east_2026": (nga_east_2026, 21),
    "nga_east_seeds_2026": (nga_east_seeds_2026, 21),
    "ag20": (nga_sub_final(pygmm.AbrahamsonGulerce2020, basin=True), 21),
    "kbcg20": (nga_sub_final(pygmm.KuehnEtAl2020, basin=True), 21),
    "psbah20": (nga_sub_final(pygmm.ParkerEtAl2020, basin=True), 21),
    "am09": (
        nshmp_subduction(pygmm.AtkinsonMacias2009, ["interface"], basin=True, m9=True),
        21,
    ),
    "zea06": (
        nshmp_subduction(
            pygmm.ZhaoEtAl2006, ["interface", "intraslab"], basin=True, m9=True
        ),
        21,
    ),
    "ask14_nshmp": (
        nga_west2_nshmp(pygmm.AbrahamsonSilvaKamai2014, "ASK_14_BASIN"),
        22,
    ),
    "bssa14_nshmp": (
        nga_west2_nshmp(pygmm.BooreStewartSeyhanAtkinson2014, "BSSA_14_BASIN"),
        105,
    ),
    "cb14_nshmp": (nga_west2_nshmp(pygmm.CampbellBozorgnia2014, "CB_14_BASIN"), 21),
    "cy14_nshmp": (nga_west2_nshmp(pygmm.ChiouYoungs2014, "CY_14_BASIN"), 24),
}


@pytest.fixture(params=sorted(MODELS))
def make(request):
    return MODELS[request.param]


@pytest.mark.parametrize(
    "period,name", [(0.01, "psa_0p010"), (1.0, "psa_1p000"), (10.0, "psa_10p000")]
)
def test_psa_name_and_period(period, name):
    assert GroundMotionModel.psa_name(period) == name
    assert GroundMotionModel.psa_period(name) == period


@pytest.mark.parametrize("name", ["psa_1p0", "psa_1p000", "psa_1", "psa_1p00000"])
def test_psa_period_accepts_any_number_of_decimals(name):
    assert GroundMotionModel.psa_period(name) == 1.0


@pytest.mark.parametrize("name", ["psa_", "psa_abc", "psa_1p0p0", "psa1p000"])
def test_invalid_psa_names_raise(name):
    with pytest.raises(ValueError):
        GroundMotionModel.psa_period(name)


def test_default_computes_all_periods(make):
    factory, n_periods = make
    m = factory()
    assert len(m.periods) == n_periods
    assert len(m.psa_ims) == n_periods
    assert m.spec_accels.shape == (300, n_periods)


def test_ngawest2_21(make):
    factory, _ = make
    full = factory()
    m = factory(ims=["psa_ngawest2_21"])
    np.testing.assert_array_equal(m.periods, NGAWEST2_21)
    assert m.psa_ims == [GroundMotionModel.psa_name(p) for p in NGAWEST2_21]
    cols = np.searchsorted(full.periods, NGAWEST2_21)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds[:, cols])
    assert m._ln_resp.shape == (300, 21)


def test_psa_all_matches_default(make):
    factory, _ = make
    full = factory()
    m = factory(ims=["psa_all"])
    np.testing.assert_array_equal(m.periods, full.periods)
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels)
    np.testing.assert_array_equal(m.ln_stds, full.ln_stds)


def test_pga_and_ngawest2_21(make):
    factory, _ = make
    full = factory()
    m = factory(ims=["pga", "psa_ngawest2_21"])
    np.testing.assert_array_equal(m.pga, full.pga)
    np.testing.assert_array_equal(m.ln_std_pga, full.ln_std_pga)
    assert m.spec_accels.shape == (300, 21)


def test_single_periods_are_sorted(make):
    factory, _ = make
    full = factory()
    m = factory(ims=["psa_1p000", "psa_0p2", "pga"])
    np.testing.assert_array_equal(m.periods, [0.2, 1.0])
    assert m.psa_ims == ["psa_0p200", "psa_1p000"]
    cols = np.searchsorted(full.periods, [0.2, 1.0])
    np.testing.assert_array_equal(m.spec_accels, full.spec_accels[:, cols])
    np.testing.assert_array_equal(m.pga, full.pga)


def test_single_period_is_one_column(make):
    factory, _ = make
    m = factory(ims=["psa_1p000"])
    assert m._ln_resp.shape == (300, 1)
    assert m.spec_accels.shape == (300, 1)


def test_interpolation_uses_computed_periods(make):
    factory, _ = make
    m = factory(ims=["psa_0p200", "psa_1p000"])
    np.testing.assert_allclose(
        m.interp_spec_accels([0.2, 1.0]), m.spec_accels, rtol=1e-12
    )
    # Outside of the computed periods, interpolation is not possible
    assert np.all(np.isnan(m.interp_spec_accels([0.1])))


def test_duplicate_and_overlapping_ims(make):
    factory, _ = make
    m = factory(ims=["psa_1p000", "psa_ngawest2_21", "psa_1p0"])
    np.testing.assert_array_equal(m.periods, NGAWEST2_21)


def test_missing_period_raises(make):
    factory, _ = make
    with pytest.raises(ValueError, match="Missing periods \\(s\\): 0.123"):
        factory(ims=["psa_0p123"])


@pytest.mark.parametrize("ims", [["psa"], ["sa_1p000"], ["PGA"]])
def test_invalid_ims_raise(make, ims):
    factory, _ = make
    with pytest.raises(ValueError, match="not a valid intensity measure"):
        factory(ims=ims)


@pytest.mark.parametrize("attr", ["periods", "psa_ims", "spec_accels", "ln_stds"])
def test_psa_not_computed_raises(make, attr):
    factory, _ = make
    m = factory(ims=["pga"])
    with pytest.raises(ValueError, match="Spectral accelerations were not computed"):
        getattr(m, attr)


@pytest.mark.parametrize(
    "model",
    [i.cls for i in pygmm.find_models() if issubclass(i.cls, GroundMotionModel)],
    ids=lambda m: m.__name__,
)
def test_all_ground_motion_models_accept_ims(model):
    assert "ims" in inspect.signature(model.__init__).parameters
