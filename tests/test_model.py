#!/usr/bin/env python
"""Test model interface using C03 model."""

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from pygmm import Campbell2003 as C03
from pygmm.model import (
    Coefficients,
    Scenario,
    as_column,
    equals,
    symmetric_branches,
    take_periods,
)


@pytest.fixture
def model():
    return C03(Scenario(mag=6.5, dist_rup=20))


def test_ln_std(model):
    assert_array_equal(model._ln_std, model.ln_stds)


def test_scenario():
    s = Scenario(mag=6, dist_jb=20)
    assert_allclose(s.mag, s["mag"])
    assert_allclose(s.dist_jb, s["dist_jb"])


@pytest.mark.parametrize(
    "attr", ["pga", "ln_std_pga", "pgv", "ln_std_pgv", "pgd", "ln_std_pgd"]
)
def test_pga(model, attr):
    with pytest.raises(NotImplementedError):
        getattr(model, attr)


@pytest.mark.parametrize(
    "values",
    [
        np.array(["SS", "RS", "NS", "U"]),
        np.array([["SS", "RS"], ["U", "RS"]]),
        np.array(["SS", "RS", "NS", "U"])[::2],
        np.array(["SS", "RS"], dtype=object),
        ["SS", "RS"],
        "RS",
        np.array([1, 2, 3]),
    ],
)
@pytest.mark.parametrize("option", ["RS", "U", "LONGER", 2])
def test_equals_matches_comparison(values, option):
    assert_array_equal(equals(values, option), np.asarray(values) == option)


def test_take_periods_and_as_column():
    one_d = np.arange(5.0)
    two_d = np.arange(10.0).reshape(2, 5)
    assert take_periods(one_d, 2) == 2.0
    assert_array_equal(take_periods(two_d, 2), [2.0, 7.0])
    assert as_column(3.0).shape == (1,)
    assert as_column(np.zeros(4)).shape == (4, 1)


def test_models_without_ims_are_unchanged(model):
    # Models that do not support ims compute all intensity measures by default
    assert model._indices is None
    assert_array_equal(model.spec_accels, np.exp(model._ln_resp))


def test_ln_pga_not_implemented(model):
    # Campbell (2003) does not provide PGA
    with pytest.raises(NotImplementedError):
        model.ln_pga


def test_ln_branches_single(model):
    [(w, ln, std)] = model.ln_branches("psa")
    assert w == 1.0
    assert_allclose(ln, np.log(model.spec_accels))
    assert_allclose(std, model.ln_stds)
    if model.INDEX_PGA is None:
        with pytest.raises(NotImplementedError):
            model.ln_branches("pga")
    with pytest.raises(ValueError):
        model.ln_branches("sa")


def test_symmetric_branches():
    w = (0.185, 0.63, 0.185)
    mu = np.array([np.log(0.1), np.log(0.5)])
    delta = np.array([0.3, 0.4])
    collapsed = np.log(
        w[0] * np.exp(mu - delta) + w[1] * np.exp(mu) + w[2] * np.exp(mu + delta)
    )
    branches = symmetric_branches(collapsed, 0.6, delta, w)
    assert [b[0] for b in branches] == list(w)
    for (_, ln, std), sign in zip(branches, [-1, 0, 1]):
        assert_allclose(ln, mu + sign * delta, rtol=1e-12)
        assert_allclose(std, [0.6, 0.6])


SUBDUCTION_SCENARIO = dict(
    mag=np.array([7.0, 8.0, 9.0]),
    dist_rup=np.array([30.0, 60.0, 100.0]),
    depth_tor=np.array([5.0, 10.0, 10.0]),
    depth_hyp=np.array([20.0, 20.0, 20.0]),
    v_s30=760.0,
    event_type="interface",
    region="cascadia",
)


@pytest.mark.filterwarnings("ignore::UserWarning")
@pytest.mark.parametrize(
    "name", ["AbrahamsonGulerce2020", "KuehnEtAl2020", "ParkerEtAl2020"]
)
def test_subduction_epistemic_branches(name):
    import pygmm

    cls = getattr(pygmm, name)
    s = Scenario(**SUBDUCTION_SCENARIO)
    epi = cls(s, ims=["pga"])
    center = cls(s, ims=["pga"], epistemic=False)
    branches = epi.ln_branches("pga")
    assert [b[0] for b in branches] == [0.185, 0.63, 0.185]
    lo, mid, hi = (b[1] for b in branches)
    assert_allclose(mid, center.ln_pga, rtol=1e-12)
    assert np.all(hi - mid > 0)
    assert_allclose(hi - mid, mid - lo, rtol=1e-10)
    w = np.array([b[0] for b in branches])
    assert_allclose(np.log(w @ np.exp(np.array([lo, mid, hi]))), epi.ln_pga, rtol=1e-12)
    for b in branches:
        assert_allclose(b[2], center.ln_std_pga)
    [(w_c, _, _)] = center.ln_branches("pga")
    assert w_c == 1.0


def test_coefficients_missing_attribute():
    """A missing name raises AttributeError, so hasattr and getattr work."""
    c = Coefficients(a=1.0)
    assert c.a == 1.0
    assert c["a"] == 1.0
    assert not hasattr(c, "__wrapped__")
    assert getattr(c, "b", None) is None
    with pytest.raises(KeyError):
        c["b"]
