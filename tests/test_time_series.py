"""Tests for pygmm.time_series."""

import pytest
from numpy.testing import assert_allclose

from pygmm.time_series import TimeSeries

from . import FPATH_DATA


@pytest.fixture()
def timeseries():
    return TimeSeries.read_at2(FPATH_DATA / "RSN4863_CHUETSU_65036EW.AT2")


def test_at2_time_step(timeseries):
    assert_allclose(timeseries.time_step, 0.01)


def test_at2_accels(timeseries):
    assert_allclose(timeseries.accels.size, 6000)
    assert_allclose(timeseries.accels[[0, -1]], [-0.2674464e-03, -0.4684600e-04])
