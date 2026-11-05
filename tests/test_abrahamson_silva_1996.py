"""Reference case from AbrahamsonSilva1996 based on spreadsheet by N.
Gregor."""

import itertools
import json
import os
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

import pygmm

fpath = os.path.join(os.path.dirname(__file__), "data", "as96.json")
with open(fpath) as fp:
    data = json.load(fp)


@pytest.mark.parametrize("case", data)
def test_calc_duration(case):
    s = pygmm.Scenario(**case["scenario"])
    m = pygmm.AbrahamsonSilva1996(s)

    nias = case["nias"]
    stds = case["stds"]

    actual = m.interp(nias, stds)

    # The spreadsheet values are rounded
    assert_allclose(actual, case["durs"], atol=0.0001, rtol=1e-4)


@pytest.mark.parametrize("site_cond", ["soil", "rock"])
def test_interp_matches_duration(site_cond):
    # The 5 to 75% duration from interp() is the duration of the model
    m = pygmm.AbrahamsonSilva1996(
        pygmm.Scenario(mag=6.5, dist_rup=20.0, site_cond=site_cond)
    )
    assert_allclose(m.interp([0.75]), [m.duration], rtol=1e-12)


# Tests of the vectorized model against scalar scenarios

MAGS = [3.5, 4.0, 5.5, 6.0, 7.5, 8.0]
DISTS = [0.0, 5.0, 10.0, 10.5, 50.0, 260.0]
# "bad" is not an option, so it is treated like rock
SITE_CONDS = ["soil", "rock", "bad"]
NIAS = [0.05, 0.1, 0.3, 0.5, 0.75, 0.9, 0.95, 0.99]
STDS = [-1.0, 0.0, 1.5]


def scalar_model(mag, dist_rup, site_cond):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return pygmm.AbrahamsonSilva1996(
            pygmm.Scenario(mag=mag, dist_rup=dist_rup, site_cond=site_cond)
        )


@pytest.fixture(scope="module")
def grid():
    rows = list(itertools.product(MAGS, DISTS, SITE_CONDS))
    mag, dist_rup, site_cond = (np.array(c) for c in zip(*rows))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = pygmm.AbrahamsonSilva1996(
            pygmm.Scenario(mag=mag, dist_rup=dist_rup, site_cond=site_cond)
        )
    return rows, m, [scalar_model(*row) for row in rows]


def test_vectorized_matches_scalar(grid):
    rows, m, scalars = grid
    n = len(rows)
    assert m.duration.shape == (n,)
    assert m.std_err.shape == (n,)
    assert_array_equal(m.duration, [s.duration for s in scalars])
    assert_array_equal(m.std_err, [s.std_err for s in scalars])

    actual = m.interp(NIAS)
    assert actual.shape == (n, len(NIAS))
    assert_array_equal(actual, [s.interp(NIAS) for s in scalars])

    actual = m.interp(NIAS, STDS)
    assert actual.shape == (n, len(STDS), len(NIAS))
    assert_array_equal(actual, [s.interp(NIAS, STDS) for s in scalars])


def test_interp_scalar_nia(grid):
    rows, m, scalars = grid
    actual = m.interp(0.5)
    assert actual.shape == (len(rows),)
    assert_array_equal(actual, [s.interp(0.5) for s in scalars])
    actual = m.interp(0.5, STDS)
    assert actual.shape == (len(rows), len(STDS))
    assert_array_equal(actual, [s.interp(0.5, STDS)[:, 0] for s in scalars])


def test_scalar_shapes_are_unchanged():
    m = scalar_model(6.0, 20.0, "soil")
    assert isinstance(m.duration, float)
    assert m.std_err == 0.55
    assert isinstance(m.std_err, float)
    assert m.interp(NIAS).shape == (len(NIAS),)
    assert m.interp(NIAS, STDS).shape == (len(STDS), len(NIAS))
    assert np.ndim(m.interp(0.5)) == 0
    assert m.interp(0.5, STDS).shape == (len(STDS), 1)


def test_scalars_broadcast_with_arrays():
    mag = np.array([5.0, 6.0, 7.0])
    m = pygmm.AbrahamsonSilva1996(
        pygmm.Scenario(mag=mag, dist_rup=30.0, site_cond="rock")
    )
    assert_array_equal(
        m.duration, [scalar_model(mg, 30.0, "rock").duration for mg in mag]
    )


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [5.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    m = pygmm.AbrahamsonSilva1996(
        pygmm.Scenario(mag=mag, dist_rup=dist_rup, site_cond="soil")
    )
    assert m.duration.shape == (3, 4)
    assert m.std_err.shape == (3, 4)
    assert m.interp(NIAS).shape == (3, 4, len(NIAS))
    assert m.interp(NIAS, STDS).shape == (3, 4, len(STDS), len(NIAS))
    assert m.duration[2, 1] == scalar_model(7.5, 20.0, "soil").duration
    assert_array_equal(
        m.interp(NIAS, STDS)[2, 1], scalar_model(7.5, 20.0, "soil").interp(NIAS, STDS)
    )


def test_invalid_site_cond_entries_use_default():
    site_cond = np.array(["soil", "rock", "bad"])
    with pytest.warns(UserWarning) as record:
        m = pygmm.AbrahamsonSilva1996(
            pygmm.Scenario(mag=6.0, dist_rup=20.0, site_cond=site_cond)
        )
    messages = [str(r.message) for r in record]
    assert sum("site_cond has 1 of 3 values" in msg for msg in messages) == 1
    assert_array_equal(
        m.duration, [scalar_model(6.0, 20.0, sc).duration for sc in site_cond]
    )


@pytest.mark.parametrize("vectorized", [False, True])
def test_nias_are_not_modified(vectorized):
    mag = np.array([6.0, 7.0]) if vectorized else 6.0
    m = pygmm.AbrahamsonSilva1996(
        pygmm.Scenario(mag=mag, dist_rup=20.0, site_cond="soil")
    )
    nias = np.array([0.05, 0.5, 0.99])
    m.interp(nias)
    m.interp(nias, STDS)
    pygmm.AbrahamsonSilva1996.calc_ln_dur_incr(nias)
    assert_array_equal(nias, [0.05, 0.5, 0.99])


def test_integer_and_list_nias():
    m = scalar_model(6.0, 20.0, "soil")
    # Values outside of the range are NaN, including integers
    assert np.isnan(m.interp(0)) and np.isnan(m.interp(1))
    actual = m.interp([0, 0.5, 1])
    assert np.isnan(actual[0]) and np.isnan(actual[2])
    assert actual[1] == m.interp(0.5)
    assert np.isnan(pygmm.AbrahamsonSilva1996.calc_ln_dur_incr(1))
