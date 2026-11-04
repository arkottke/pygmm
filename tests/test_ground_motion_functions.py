"""Tests for the function-style relationships in pygmm.ground_motion."""

import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm.ground_motion import (
    calc_damping_scaling_rea15,
    calc_period_rea05,
    calc_pulse_proportion_hea12,
)

from . import FPATH_DATA


@pytest.mark.parametrize(
    "kind,mag,dist_rup,site_class,directivity,expected",
    [
        # Example from Rathje et. al (2014)
        ("period_mean", 7, 5, "b", False, 0.45),
        # Example from Rathje and Antonakos (2011)
        ("period_mean", 8, 2, "b", False, 0.46),
    ],
)
def test_calc_period_rea05(kind, mag, dist_rup, site_class, directivity, expected):
    actual = calc_period_rea05(kind, mag, dist_rup, site_class, directivity)
    assert_allclose(actual, expected, rtol=0.01)


with open(FPATH_DATA / "rezaeian_et_al-2014.json") as fp:
    REA15_CASES = json.load(fp)


@pytest.mark.parametrize("case", REA15_CASES)
def test_calc_damping_scaling_rea15(case):
    ret = calc_damping_scaling_rea15(
        case["damping"], case["mag"], case["dist_rup"], comp=case["comp"]
    )

    assert_allclose(ret.dsf, case["dsf"], rtol=0.01)
    assert_allclose(ret.ln_std, case["ln_std"], rtol=0.01)


def test_calc_damping_scaling_rea15_periods():
    full = calc_damping_scaling_rea15(2, 6, 10)
    periods = full.period[[3, 10]]
    ret = calc_damping_scaling_rea15(2, 6, 10, periods=periods)
    assert_allclose(ret.dsf, full.dsf[[3, 10]])


def test_calc_damping_scaling_rea15_bounds():
    with pytest.raises(ValueError, match="damping"):
        calc_damping_scaling_rea15(50, 6, 10)


def test_calc_pulse_proportion_hea12():
    prop = calc_pulse_proportion_hea12(np.array([0.0, 10.0, 30.0]), 1.0)
    assert np.all((0 < prop) & (prop < 1))
    # Proportion decreases with distance
    assert np.all(np.diff(prop) < 0)
