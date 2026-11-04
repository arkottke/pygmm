"""Tests for pygmm.kappa."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm.kappa import calc_site_atten_vhea11, calc_site_atten_xr21


@pytest.mark.parametrize(
    "func,v_s30,expected",
    [
        (calc_site_atten_vhea11, 760.0, np.exp(3.49 - 1.062 * np.log(760.0))),
        (calc_site_atten_xr21, 760.0, 0.0386),
    ],
)
def test_median(func, v_s30, expected):
    assert_allclose(func(v_s30), expected, rtol=1e-3)


@pytest.mark.parametrize(
    "func,ln_std", [(calc_site_atten_vhea11, 0.55), (calc_site_atten_xr21, 0.3)]
)
def test_stats(func, ln_std):
    ln_mean, actual = func(760.0, stats=True)
    assert_allclose(np.exp(ln_mean), func(760.0))
    assert actual == ln_std


def test_xr21_clips_v_s30():
    assert_allclose(calc_site_atten_xr21(3000.0), calc_site_atten_xr21(2000.0))
    assert_allclose(calc_site_atten_xr21(100.0), calc_site_atten_xr21(155.0))


def test_vectorized():
    v_s30 = np.array([300.0, 760.0, 1500.0])
    kappa = calc_site_atten_xr21(v_s30)
    assert kappa.shape == v_s30.shape
    # Kappa decreases with increasing velocity
    assert np.all(np.diff(kappa) < 0)
