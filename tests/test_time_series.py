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


@pytest.mark.parametrize(
    "header,time_step",
    [
        ("NPTS=  6000, DT= .0100 SEC", 0.01),
        ("NPTS=  16396, DT= 0.005 SEC", 0.005),
        ("NPTS=6000, DT=0.0100 SEC", 0.01),
        ("   4096   0.0100   NPTS, DT", 0.01),
    ],
)
def test_at2_header_formats(tmp_path, header, time_step):
    fpath = tmp_path / "test.AT2"
    fpath.write_text(
        "PEER NGA\nEVENT\nACCELERATION TIME SERIES IN UNITS OF G\n"
        f"{header}\n 0.1 0.2 0.3\n"
    )
    ts = TimeSeries.read_at2(fpath)
    assert_allclose(ts.time_step, time_step)
    assert_allclose(ts.accels, [0.1, 0.2, 0.3])
