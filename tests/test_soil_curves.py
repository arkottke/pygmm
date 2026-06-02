"""Tests for pygmm.soil_curves models and backward-compat top-level imports."""

import numpy as np
import pytest
from numpy.testing import assert_array_less

from pygmm.contracts import NonlinearSoilCurves
from pygmm.soil_curves import (
    AlemuEtAlSoilType,
    DarendeliSoilType,
    KishidaSoilType,
    MenqSoilType,
    RollinsEtAlSoilType,
    WangSoilType,
)


def _check_curves(c: NonlinearSoilCurves):
    """Shared sanity checks for any NonlinearSoilCurves instance."""
    assert isinstance(c, NonlinearSoilCurves)
    assert c.strains.shape == c.mod_reduc.shape == c.damping.shape
    assert np.all(np.isfinite(c.strains))
    assert np.all(np.isfinite(c.mod_reduc))
    assert np.all(np.isfinite(c.damping))
    # G/Gmax ∈ (0, 1]
    assert_array_less(0, c.mod_reduc)
    assert c.mod_reduc[0] <= 1.0 + 1e-10
    # Damping ≥ 0
    assert np.all(c.damping >= 0)
    # damping_min is positive
    assert c.damping_min > 0


def test_darendeli_defaults():
    m = DarendeliSoilType()
    _check_curves(m.curves())


def test_darendeli_params():
    m = DarendeliSoilType(unit_wt=18.0, plas_index=15, ocr=2, stress_mean=50.0)
    c = m.curves()
    _check_curves(c)
    assert c.unit_wt == 18.0


def test_menq_defaults():
    m = MenqSoilType()
    _check_curves(m.curves())


def test_menq_params():
    m = MenqSoilType(coef_unif=5, diam_mean=2, stress_mean=200)
    _check_curves(m.curves())


def test_rollins_without_coef_unif():
    m = RollinsEtAlSoilType(stress_mean=100.0)
    _check_curves(m.curves())


def test_rollins_with_coef_unif():
    m = RollinsEtAlSoilType(stress_mean=100.0, coef_unif=5.0)
    _check_curves(m.curves())


def test_alemu_defaults():
    m = AlemuEtAlSoilType()
    _check_curves(m.curves())


def test_alemu_params():
    m = AlemuEtAlSoilType(plas_index=15, ocr=2, stress_mean=60.0)
    _check_curves(m.curves())


def test_wang_clean_sand():
    m = WangSoilType("clean_sand_and_gravel", stress_mean=100.0)
    _check_curves(m.curves())


def test_wang_nonplastic_silt():
    m = WangSoilType("nonplastic_silty_sand", stress_mean=80.0, void_ratio=0.7)
    _check_curves(m.curves())


def test_wang_clayey():
    m = WangSoilType("clayey_soil", stress_mean=60.0, void_ratio=0.8)
    _check_curves(m.curves())


def test_wang_bad_group():
    with pytest.raises(ValueError):
        WangSoilType("unknown_group", stress_mean=100.0)


def test_kishida_defaults():
    m = KishidaSoilType()
    c = m.curves()
    _check_curves(c)
    assert c.unit_wt > 0


def test_kishida_with_unit_wt():
    m = KishidaSoilType(unit_wt=14.0, stress_vert=50.0, organic_content=30)
    c = m.curves()
    _check_curves(c)
    assert c.unit_wt == 14.0


def test_backward_compat_imports():
    import pygmm

    assert pygmm.DarendeliSoilType is DarendeliSoilType
    assert pygmm.MenqSoilType is MenqSoilType
    assert pygmm.RollinsEtAlSoilType is RollinsEtAlSoilType
    assert pygmm.AlemuEtAlSoilType is AlemuEtAlSoilType
    assert pygmm.WangSoilType is WangSoilType
    assert pygmm.KishidaSoilType is KishidaSoilType
    assert pygmm.kea16_profile is not None
