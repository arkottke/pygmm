"""Tests for pygmm.slope_displacement."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm import slope_displacement as sd
from pygmm.slope_displacement.watson_lamprey_abrahamson_2006 import (
    _calc_wla06_ln_a_rms,
    _calc_wla06_ln_dur_key,
)
from pygmm.time_series import TimeSeries

from . import FPATH_DATA


@pytest.mark.parametrize(
    "yield_coef,pga,mag,pgv,expected",
    [
        # Examples from Rathje et. al (2014)
        (0.1, 0.54, 6.75, None, 43),
        (0.1, 0.54, None, 42, 29),
        (0.1, 0.88, 6.75, None, 113),
        (0.1, 0.88, None, 71, 81),
        # Examples from Rathje & Antonakos (2011)
        # These do not work...
        # (0.1, 0.38, 8.0, None, 63.1),
        # (0.1, 0.38, None, 80, 36.9),
    ],
)
def test_calc_disp_rs08(yield_coef, pga, mag, pgv, expected):
    disp = sd.calc_disp_rs08(yield_coef, pga, mag=mag, pgv=pgv)
    # Compare to 2 cm of displacement
    assert_allclose(disp, expected, atol=2)


@pytest.mark.parametrize(
    "yield_coef,pga,period_slide,period_mean,mag,pgv,expected",
    [
        # Examples from Rathje et. al (2014)
        # Values read from Figure 8c
        (0.05, 0.35, 0.2, 0.45, 7.0, None, 87.5),
        (0.10, 0.35, 0.2, 0.45, 7.0, None, 22.0),
        (0.05, 0.35, 0.2, 0.45, None, 30, 43.5),
        (0.10, 0.35, 0.2, 0.45, None, 30, 10.5),
        # Examples from Rathje & Antonakos (2011)
        # These do not work...
        # (0.1, 0.48, 0.2, 0.46, 8.0, None, 126),
        # (0.1, 0.48, 0.2, 0.46, None, 74, 49)
    ],
)
def test_calc_disp_ra11(yield_coef, pga, period_slide, period_mean, mag, pgv, expected):
    disp = sd.calc_disp_ra11(
        yield_coef, pga, period_slide, period_mean, mag=mag, pgv=pgv
    )
    assert_allclose(disp, expected, atol=1)


def test_calc_wla06_ln_a_rms():
    a_rms = np.exp(_calc_wla06_ln_a_rms(0.35))
    assert_allclose(a_rms, 0.1067, atol=0.0001)


def test_calc_wla06_ln_dur_key():
    dur_ky = np.exp(_calc_wla06_ln_dur_key(0.11, 0.35, 0.18, 6.5))
    assert_allclose(dur_ky, 0.5896, atol=0.0001)


def test_calc_disp_wla06():
    median = sd.calc_disp_wla06(0.11, 0.35, 0.18, 6.5)
    assert_allclose(median, 3.001, atol=0.001)


# Test values from Table 1 in BT07
@pytest.mark.parametrize(
    "yield_coef,period_slide,psa_dts,expected",
    [
        (0.35, 0.45, 0.43, 0.85),
        (0.08, 0.00, 0.24, 0.1),
        (0.14, 0.33, 0.94, 0.0),
    ],
)
def test_calc_prob_disp_bt07(yield_coef, period_slide, psa_dts, expected):
    # Convert from probability of zero displacement to probability of non-zero
    # displacement
    expected = 1 - expected

    calc = sd.calc_prob_disp_bt07(yield_coef, period_slide, psa_dts)
    assert_allclose(calc, expected, atol=0.02, rtol=0.02)


# Test values from Table 1 in BT07. Authors provided a range of values, which
# are interpreted to be log-normally distributed. The median is computed and
# tested against
@pytest.mark.parametrize(
    "yield_coef,period_slide,psa_dts,expected",
    [
        (0.08, 0.00, 0.24, np.sqrt(4 * 15)),
        (0.14, 0.33, 0.94, np.sqrt(20 * 70)),
    ],
)
def test_calc_disp_bt07(yield_coef, period_slide, psa_dts, expected):
    median = sd.calc_disp_bt07(yield_coef, period_slide, psa_dts)
    assert_allclose(median, expected, rtol=0.10)


@pytest.mark.parametrize(
    "invert,yield_coef,expected",
    [
        (False, 0.05, 69.4),
        (True, 0.05, 66.1),
        (False, 0.10, 28.5),
        (True, 0.10, 30.5),
    ],
)
def test_calc_block_displacement(invert, yield_coef, expected):
    time_step = 0.01
    # Acceleration data in g
    accels = np.loadtxt(str(FPATH_DATA / "kobe-nis-000.dat"))
    yield_coef = yield_coef

    disps = sd.calc_rigid_disp(time_step, accels, yield_coef, invert)["disps"]
    disp_max = disps[-1]
    # Large rtol because of slightly different implementations
    assert_allclose(disp_max, expected, rtol=0.04)


# Test values from Table 2 of Bray et al. (2017). Coastline slope and Nishigo dam
@pytest.mark.parametrize(
    "yield_coef,period_slide,psa_dts,mag,expected",
    [
        (0.10, 0.6, 0.25, 8.0, np.sqrt(3 * 12)),
        (0.26, 0.15, 1.51, 9.0, np.sqrt(14 * 58)),
    ],
)
def test_calc_disp_bea17(yield_coef, period_slide, psa_dts, mag, expected):
    ret = sd.calc_disp_bea17(yield_coef, period_slide, psa_dts, mag, stats=False)
    # Limited number of digits in table make testing difficult
    assert_allclose(ret, expected, rtol=0.50)


# Test values from Table 2 of Bray et al. (2017). Esperanza and Tutuven dams.
@pytest.mark.parametrize(
    "yield_coef,period_slide,psa_dts,expected",
    [
        (0.24, 0.40, 0.43, 0.50),
        (0.39, 0.15, 0.75, 0.60),
    ],
)
def test_calc_prob_disp_bea17(yield_coef, period_slide, psa_dts, expected):
    ret = sd.calc_prob_disp_bea17(yield_coef, period_slide, psa_dts)
    # Reported values are probability of no displacement
    # Limited number of digits in table make testing difficult
    assert_allclose(ret, 1 - expected, rtol=0.25)


# Test values from Figure 13 of Cho and Rathje (20222)
@pytest.mark.parametrize(
    "yield_coef,period_slide,height_ratio,pgv,expected",
    [
        (0.20, 0.30, 0.7, 20, 5.7),
        (0.20, 0.30, 0.5, 20, 3.1),
        (0.20, 0.45, 0.5, 20, 6.8),
        (0.20, 0.60, 0.5, 20, 8.3),
        (0.20, 0.30, 0.7, 100, 50),
        (0.20, 0.30, 0.5, 100, 22),
        (0.20, 0.45, 0.5, 100, 53),
        (0.20, 0.60, 0.5, 100, 58),
    ],
)
def test_calc_disp_cr22_pgv(yield_coef, period_slide, height_ratio, pgv, expected):
    ret = sd.calc_disp_cr22(yield_coef, period_slide, height_ratio, pgv=pgv)
    # Reported values are probability of no displacement
    # Limited number of digits in table make testing difficult
    assert_allclose(ret, expected, rtol=0.05)


# Test values from Figure 14 of Cho and Rathje (2022)
@pytest.mark.parametrize(
    "yield_coef,period_slide,height_ratio,mag,pga,expected",
    [
        (0.10, 0.30, 0.7, 7.5, 0.7, 71),
        (0.20, 0.30, 0.7, 7.5, 0.7, 36),
        (0.40, 0.30, 0.7, 7.5, 0.7, 19),
    ],
)
def test_calc_disp_cr22_mag_pga(
    yield_coef, period_slide, height_ratio, mag, pga, expected
):
    ret = sd.calc_disp_cr22(yield_coef, period_slide, height_ratio, mag=mag, pga=pga)
    # Reported values are probability of no displacement
    # Limited number of digits in table make testing difficult
    assert_allclose(ret, expected, rtol=0.05)


# Tests from Figure 11 in the Cho and Rathje (2021)
@pytest.mark.parametrize(
    "yield_coef,period_slide,height_ratio,pgv,expected",
    [
        (0.1, 0.4, 0.2, 20, 6.0),
        (0.2, 0.4, 0.8, 20, 5.6),
        (0.4, 0.4, 1.6, 20, 3.1),
        (0.1, 0.3, 0.2, 20, 3.9),
        (0.2, 0.4, 0.2, 20, 5.3),
        (0.4, 0.5, 0.2, 20, 6.3),
    ],
)
def test_calc_disp_cr21(yield_coef, period_slide, height_ratio, pgv, expected):
    ret = sd.calc_disp_cr21(pgv, yield_coef, period_slide, height_ratio)
    # Reported values are probability of no displacement
    # Limited number of digits in table make testing difficult
    assert_allclose(ret, expected, rtol=0.05)


def test_ha19_trans_func():
    pytest.importorskip("pyrotd")
    ha19 = sd.HaleAbrahamson19(100, freq_nat=5.0)
    np.testing.assert_allclose(ha19.shear_vel, 1500)

    result = ha19.transfer_function([0.1], 0.01)
    np.testing.assert_allclose(result.tf, [1])


# FIXME: Add test value
def test_ha19_displacement():
    pytest.importorskip("pyrotd")
    # Acceleration data in g
    ts = TimeSeries.read_at2(FPATH_DATA / "RSN4863_CHUETSU_65036EW.AT2")

    ha19 = sd.HaleAbrahamson19(50, shear_vel=800)
    _ = ha19(ts.time_step, ts.accels, 0.1)


@pytest.mark.parametrize(
    "pgv,pga,mag,expected",
    # Test values selected from Figure 3.13
    [
        (10, None, None, 5.9),
        (None, 0.1, 7, 6.1),
    ],
)
def test_calc_disp_cr20(pgv, pga, mag, expected):
    disp = sd.calc_disp_cr20(0.1, 0.6, pga, mag, pgv)
    np.testing.assert_allclose(disp, expected, rtol=0.01)


@pytest.mark.parametrize(
    "func,args,kwds",
    [
        (sd.calc_disp_rs08, (0.1, 0.54), {"pgv": 42}),
        (sd.calc_disp_ra11, (0.1, 0.35, 0.2, 0.45), {"mag": 7.0}),
        (sd.calc_disp_wla06, (0.11, 0.35, 0.18, 6.5), {}),
        (sd.calc_disp_bt07, (0.14, 0.33, 0.94), {}),
        (sd.calc_disp_bea17, (0.10, 0.6, 0.25, 8.0), {}),
        (sd.calc_disp_cr20, (0.1, 0.6), {"pgv": 10}),
        (sd.calc_disp_cr21, (20, 0.1, 0.4, 0.2), {}),
        (sd.calc_disp_cr22, (0.2, 0.3, 0.5), {"pgv": 20}),
        (sd.calc_dam_period_pk19, (50,), {}),
    ],
)
def test_stats(func, args, kwds):
    ln_mean, ln_std = func(*args, stats=True, **kwds)
    assert_allclose(np.exp(ln_mean), func(*args, **kwds))
    assert np.isfinite(ln_std) and ln_std > 0


def test_calc_disp_wla06_below_yield():
    assert np.isnan(sd.calc_disp_wla06(0.5, 0.35, 0.18, 6.5))
    disps = sd.calc_disp_wla06(0.2, np.array([0.1, 0.35]), 0.18, 6.5)
    assert np.isnan(disps[0]) and np.isfinite(disps[1])


def test_calc_dam_period_pk19_aliases():
    assert_allclose(
        sd.calc_dam_period_pk19(50, "trans"), sd.calc_dam_period_pk19(50, "transverse")
    )


def test_ha19_worked_example():
    """Example calculation in Section 4.5.2 of Hale (2019)."""
    # 100 ft high, 1160 ft/s average shear-wave velocity, SA(T_NL) = 1.392 g
    ha19 = sd.HaleAbrahamson19(100, shear_vel=1160)
    ln_sa = np.log(1.392)

    np.testing.assert_allclose(1 / ha19.freq_nat, 0.259, atol=5e-4)  # T_o
    np.testing.assert_allclose(ha19.period_nl, 0.508, atol=5e-4)

    damping, period_ratio = ha19._calc_params(ln_sa)
    np.testing.assert_allclose(damping, 0.363, atol=5e-4)
    np.testing.assert_allclose(period_ratio, 2.601, atol=5e-4)

    period_1 = period_ratio / ha19.freq_nat
    np.testing.assert_allclose(period_1, 0.673, atol=5e-4)
    omega_1 = 2 * np.pi / period_1
    np.testing.assert_allclose(omega_1, 9.341, atol=5e-3)
    np.testing.assert_allclose(omega_1 / np.sqrt(1 - damping**2), 10.023, atol=5e-3)


def test_ha19_trans_func_shape():
    ha19 = sd.HaleAbrahamson19(100, shear_vel=1160)
    ln_sa = np.log(1.392)
    damping, period_ratio = ha19._calc_params(ln_sa)
    freq_first = ha19.freq_nat / period_ratio
    freq_nat = freq_first / np.sqrt(1 - damping**2)

    freqs = np.array([0.1, 0.5 * freq_first, freq_first, 2 * freq_first, 20.0])
    trans_func = ha19.transfer_function(freqs, np.exp(ln_sa)).tf

    # Equation 4.1
    w, wn = 2 * np.pi * freqs, 2 * np.pi * freq_nat
    expected = 0.85 * (1 + (1j * w) ** 2 / (w**2 - wn**2 - 2j * damping * wn * w))
    np.testing.assert_allclose(trans_func[2:], expected[2:])
    # Amplitude is at least one below the first-mode frequency
    assert np.all(np.abs(trans_func[:2]) >= 1.0)
    # Resonance near the first-mode frequency and decay to alpha-scaled zero
    assert np.abs(trans_func[2]) > np.abs(trans_func[1])
    assert np.abs(trans_func[-1]) < 0.85


def test_ha19_transfer_function_result():
    """The result carries the parameters of the worked example, Section 4.5.2."""
    ha19 = sd.HaleAbrahamson19(100, shear_vel=1160)
    freqs = np.logspace(-1, 1.5, 50)
    result = ha19.transfer_function(freqs, 1.392)

    np.testing.assert_allclose(result.freqs, freqs)
    assert result.tf.shape == freqs.shape and np.iscomplexobj(result.tf)
    np.testing.assert_allclose(result.damping, 0.363, atol=5e-4)
    np.testing.assert_allclose(result.period_ratio, 2.601, atol=5e-4)
    np.testing.assert_allclose(1 / result.freq_first, 0.673, atol=5e-4)
    np.testing.assert_allclose(2 * np.pi * result.freq_nat_eff, 10.023, atol=5e-3)
    # The first-mode frequency is below the natural frequency of the undamped system
    assert result.freq_first < result.freq_nat_eff < ha19.freq_nat


def test_ha19_transfer_function_strengthens_with_shaking():
    """Stronger shaking lowers the frequency and raises the damping."""
    ha19 = sd.HaleAbrahamson19(100, shear_vel=1160)
    weak, strong = (
        ha19.transfer_function([1.0], 0.05),
        ha19.transfer_function([1.0], 2.0),
    )
    assert strong.freq_first < weak.freq_first
    assert strong.damping > weak.damping


@pytest.mark.parametrize("spec_accel", [0.0, -1.0, np.nan])
def test_ha19_transfer_function_rejects_bad_spec_accel(spec_accel):
    with pytest.raises(ValueError, match="spec_accel"):
        sd.HaleAbrahamson19(100, shear_vel=1160).transfer_function([1.0], spec_accel)


def test_ha19_from_metric():
    metric = sd.HaleAbrahamson19.from_metric(30.48, shear_vel=353.568)
    feet = sd.HaleAbrahamson19(100.0, shear_vel=1160.0)
    np.testing.assert_allclose(metric.height, feet.height)
    np.testing.assert_allclose(metric.shear_vel, feet.shear_vel)
    np.testing.assert_allclose(metric.freq_nat, feet.freq_nat)
    np.testing.assert_allclose(
        metric.transfer_function([1.0], 1.392).tf,
        feet.transfer_function([1.0], 1.392).tf,
    )
    # The natural frequency does not depend on the units
    by_freq = sd.HaleAbrahamson19.from_metric(30.48, freq_nat=5.0)
    assert by_freq.freq_nat == 5.0
    np.testing.assert_allclose(by_freq.height, 100.0)
