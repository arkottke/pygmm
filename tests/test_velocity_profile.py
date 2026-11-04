"""Tests for pygmm.velocity_profile models."""

from __future__ import annotations

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm.contracts import VelocityProfile
from pygmm.velocity_profile import (
    calc_density_bea16,
    calc_vel_comp_bea16,
    calc_vel_shear_spt_wds12,
    kea16_profile,
    sa18_profile,
)
from pygmm.velocity_profile.boore_2016 import (
    bj97gr760_density,
    bj97gr760_profile,
    bj97gr_profile,
    bj97gvhr_profile,
)
from pygmm.velocity_profile.boore_thompson_cadet_2011 import btc11_profile


def test_kea16_profile_smoke():
    depth = np.linspace(0, 100, 11)
    prof = kea16_profile(depth, vs30=400.0, region="california")

    assert isinstance(prof, VelocityProfile)
    assert prof.depth.shape == prof.vs_median.shape == prof.std_vs_ln.shape
    assert np.all(np.isfinite(prof.vs_median))
    assert np.all(prof.vs_median > 0)
    assert np.all(prof.std_vs_ln > 0)


@pytest.mark.parametrize("region", ["california", "japan"])
def test_kea16_profile_known_values(region):
    # Regression values computed directly from the Kamai et al. (2016)
    # equations for Vs30 = 400 m/s.
    depth = np.array([0.0, 5, 10, 30, 50])

    expected_vs_median = {
        "california": np.array(
            [204.0, 297.564475, 365.78945191, 534.40152007, 634.85789197]
        ),
        "japan": np.array(
            [120.0, 285.46213319, 398.39292794, 661.53223732, 811.92318939]
        ),
    }
    expected_std_vs_ln = {"california": 0.275, "japan": 0.4}

    prof = kea16_profile(depth, vs30=400.0, region=region)

    assert_allclose(prof.vs_median, expected_vs_median[region], rtol=1e-6)
    assert_allclose(
        prof.std_vs_ln, np.full_like(depth, expected_std_vs_ln[region]), rtol=1e-6
    )


def test_kea16_profile_region_case_insensitive():
    depth = np.array([0.0, 10, 30])
    prof_lower = kea16_profile(depth, vs30=400.0, region="california")
    prof_upper = kea16_profile(depth, vs30=400.0, region="CALIFORNIA")

    assert_allclose(prof_lower.vs_median, prof_upper.vs_median)
    assert_allclose(prof_lower.std_vs_ln, prof_upper.std_vs_ln)


def test_kea16_profile_invalid_region_raises():
    with pytest.raises(ValueError, match="region must be"):
        kea16_profile([0, 10], vs30=400.0, region="texas")


@pytest.mark.parametrize("vs30", [100.0, 900.0])
def test_kea16_profile_warns_outside_recommended_range(vs30):
    with pytest.warns(UserWarning, match="recommended"):
        kea16_profile([0, 10], vs30, region="california")


def test_sa18_profile_shape_and_monotonicity():
    depth = np.array([0.0, 1.25, 2.5, 5, 10, 30, 50])
    prof = sa18_profile(depth, vs30=300.0)

    assert isinstance(prof, VelocityProfile)
    assert prof.depth.shape == prof.vs_median.shape == prof.std_vs_ln.shape
    assert np.all(np.isfinite(prof.vs_median))
    assert np.all(np.isfinite(prof.std_vs_ln))
    assert np.all(prof.std_vs_ln > 0)
    # Vs is constant over the top 2.5 m and monotonically increasing below.
    assert prof.vs_median[0] == prof.vs_median[1] == prof.vs_median[2]
    assert np.all(np.diff(prof.vs_median[2:]) > 0)
    assert prof.region == "california"


def test_sa18_profile_known_values():
    # Regression values computed directly from equations (3), (4), and (9)
    # of Shi and Asimaki (2018) for Vs30 = 300 m/s.
    depth = np.array([0.0, 1.25, 2.5, 5, 10, 30, 50])
    prof = sa18_profile(depth, vs30=300.0)

    expected_vs_median = np.array(
        [
            205.3928,
            205.3928,
            205.3928,
            239.54615298,
            288.92766567,
            406.23029423,
            480.82679063,
        ]
    )
    expected_std_vs_ln = np.array(
        [
            0.31533521,
            0.32459952,
            0.33382248,
            0.30435504,
            0.28108959,
            0.27930994,
            0.30199437,
        ]
    )

    assert_allclose(prof.vs_median, expected_vs_median, rtol=1e-5)
    assert_allclose(prof.std_vs_ln, expected_std_vs_ln, rtol=1e-5)


@pytest.mark.parametrize("vs30", [100.0, 1100.0])
def test_sa18_profile_warns_outside_recommended_range(vs30):
    with pytest.warns(UserWarning, match="recommended"):
        sa18_profile([0, 10], vs30)


def test_sa18_profile_raises_above_max_vs30():
    with pytest.raises(ValueError, match="must not exceed"):
        sa18_profile([0, 10], 1600.0)


def test_bj97gr760_profile_regression():
    # Regression values interpolated from Table 1 of Boore (2016).
    depth = np.array([0.0, 1.0, 30.0, 300.0, 1000.0, 7850.0])
    prof = bj97gr760_profile(depth)

    assert isinstance(prof, VelocityProfile)
    expected_vs_median = np.array([314.0, 427.0, 1020.0, 2147.0, 2663.0, 3519.0])
    assert_allclose(prof.vs_median, expected_vs_median, rtol=1e-6)
    assert_allclose(prof.std_vs_ln, np.zeros_like(expected_vs_median))


def test_bj97gr760_density_regression():
    depth = np.array([0.0, 1.0, 30.0, 300.0, 1000.0, 7850.0])
    density = bj97gr760_density(depth)

    expected_density = np.array([1.934, 1.989, 2.184, 2.426, 2.535, 2.723])
    assert_allclose(density, expected_density, rtol=1e-6)


def test_bj97gvhr_profile_regression():
    depth = np.array([0.0, 30.0, 300.0, 750.0])
    prof = bj97gvhr_profile(depth)

    assert isinstance(prof, VelocityProfile)
    expected_vs_median = np.array([2768.0, 2792.0, 2993.0, 3260.0])
    assert_allclose(prof.vs_median, expected_vs_median, rtol=1e-6)


def test_bj97gvhr_profile_raises_beyond_max_depth():
    with pytest.raises(ValueError, match="must not exceed"):
        bj97gvhr_profile([800.0])


def test_bj97gr_profile_regression():
    # BJ97gr reconstructed by inverting the Boore (2016) slowness
    # interpolation of BJ97gr760 and BJ97gvhr.
    depth = np.array([0.0, 30.0, 300.0, 750.0])
    prof = bj97gr_profile(depth)

    assert isinstance(prof, VelocityProfile)
    expected_vs_median = np.array(
        [245.65051193, 850.57742584, 1972.05860579, 2363.4912982]
    )
    assert_allclose(prof.vs_median, expected_vs_median, rtol=1e-6)


def test_bj97gr_profile_matches_vs30_target():
    # Integrating the reconstructed BJ97gr slowness to 30 m should recover
    # the reference V-bar_S(30 m) = 618 m/s to within the discretization
    # error of the tabulated profiles.
    depth = np.linspace(0, 30, 3000)
    prof = bj97gr_profile(depth)
    vs30_bar = 30.0 / np.trapezoid(1 / prof.vs_median, depth)

    assert vs30_bar == pytest.approx(618.0, rel=1e-2)


def test_bj97gr_profile_raises_beyond_max_depth():
    with pytest.raises(ValueError, match="must not exceed"):
        bj97gr_profile([800.0])


def test_btc11_profile_regression():
    depth = np.array([5.0, 10.0, 20.0, 29.0])
    prof = btc11_profile(depth, vs30=400.0)

    assert isinstance(prof, VelocityProfile)
    expected_vs_median = np.array(
        [50.08434349, 43.22469183, 162.17888979, 369.25337284]
    )
    expected_std_vs_ln = np.array([0.119, 0.084, 0.035, 0.003])
    assert_allclose(prof.vs_median, expected_vs_median, rtol=1e-6)
    assert_allclose(prof.std_vs_ln, expected_std_vs_ln, rtol=1e-6)


@pytest.mark.parametrize("depth", [3.0, 30.0])
def test_btc11_profile_raises_outside_depth_range(depth):
    with pytest.raises(ValueError, match="must be between"):
        btc11_profile([depth], vs30=400.0)


def test_bj97gr_profile_matches_analytic():
    # Boore & Joyner (1997) define BJ97gr by piecewise power laws (valid to
    # 8 km); the reconstruction should agree over its 750 m range.
    depth = np.array([0.5, 5, 20, 30, 50, 100, 190, 300, 500, 750.0])
    z = depth / 1000
    expected = 1000 * np.piecewise(
        z,
        [z <= 0.001, (0.001 < z) & (z <= 0.030), (0.030 < z) & (z <= 0.190), 0.190 < z],
        [
            0.245,
            lambda z: 2.206 * z**0.272,
            lambda z: 3.542 * z**0.407,
            lambda z: 2.505 * z**0.199,
        ],
    )
    assert_allclose(bj97gr_profile(depth).vs_median, expected, rtol=0.005)


def test_calc_vel_comp_bea16_regression():
    assert_allclose(calc_vel_comp_bea16([200.0, 760.0]), [1329.1, 2168.3], atol=0.1)


def test_calc_density_bea16_regression():
    # One value from each of the three velocity ranges.
    assert_allclose(
        calc_density_bea16([200.0, 760.0, 4000.0]), [1.8245, 2.1114, 2.9496], atol=1e-4
    )


@pytest.mark.parametrize(
    "soil_type,age,expected",
    [
        ("all", None, 203.4),
        ("sand", "holocene", 155.6),
        ("gravel", "pleistocene", 466.6),
    ],
)
def test_calc_vel_shear_spt_wds12(soil_type, age, expected):
    vs = calc_vel_shear_spt_wds12(20, 1.0, soil_type=soil_type, age=age)
    assert_allclose(vs, expected, atol=0.1)


def test_calc_vel_shear_spt_wds12_invalid_soil_type():
    with pytest.raises(ValueError, match="soil_type"):
        calc_vel_shear_spt_wds12(20, 1.0, soil_type="rock")
