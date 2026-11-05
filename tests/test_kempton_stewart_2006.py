import itertools
import warnings

import numpy as np
import pytest

from pygmm import Scenario
from pygmm.ground_motion import KemptonStewart2006


def test_run():
    """Simple test to make sure we can get outputs."""

    s = Scenario(
        mag=6,
        dist_rup=50,
        v_s30=300,
    )

    m = KemptonStewart2006(s)

    m.duration

    m.std_err


# Tests of the vectorized model against scalar scenarios

MAGS = [4.5, 5.0, 6.0, 6.8, 7.6, 8.0]
DISTS = [0.0, 10.0, 50.0, 250.0]
V_S30S = [150.0, 200.0, 500.0, 1000.0, 1200.0]
FIELDS = ("D_5t75a", "D_5t95a", "D_5t75v", "D_5t95v")


def scalar_model(mag, dist_rup, v_s30):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return KemptonStewart2006(Scenario(mag=mag, dist_rup=dist_rup, v_s30=v_s30))


def as_array(recarray):
    """Stack the fields of a duration recarray along the last axis."""
    return np.stack([recarray[n] for n in recarray.dtype.names], axis=-1)


def test_vectorized_matches_scalar():
    rows = list(itertools.product(MAGS, DISTS, V_S30S))
    mag, dist_rup, v_s30 = (np.array(c) for c in zip(*rows))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = KemptonStewart2006(Scenario(mag=mag, dist_rup=dist_rup, v_s30=v_s30))
    scalars = [scalar_model(*row) for row in rows]
    for attr in ["duration", "std_err"]:
        expected = np.array([as_array(getattr(s, attr)) for s in scalars])
        actual = as_array(getattr(m, attr))
        assert actual.shape == (len(rows), 4)
        np.testing.assert_array_equal(actual, expected)


def test_vectorized_shapes():
    n = 3
    m = KemptonStewart2006(
        Scenario(
            mag=np.array([5.5, 6.5, 7.5]),
            dist_rup=np.array([5.0, 20.0, 80.0]),
            v_s30=300.0,
        )
    )
    assert m._ln_dur.shape == (n, 4)
    for rec in [m.duration, m.std_err]:
        assert rec.dtype.names == FIELDS
        assert rec.shape == (n,)
        for name in FIELDS:
            assert rec[name].shape == (n,)


def test_scalar_shapes_are_unchanged():
    m = scalar_model(6.0, 50.0, 300.0)
    assert m._ln_dur.shape == (4,)
    for rec in [m.duration, m.std_err]:
        assert isinstance(rec, np.recarray)
        assert rec.shape == ()
        # Fields are 0-d arrays
        assert rec.D_5t95a.shape == ()


def test_scalars_broadcast_with_arrays():
    v_s30 = np.array([200.0, 400.0, 800.0])
    m = KemptonStewart2006(Scenario(mag=6.5, dist_rup=30.0, v_s30=v_s30))
    expected = np.array([as_array(scalar_model(6.5, 30.0, v).duration) for v in v_s30])
    np.testing.assert_array_equal(as_array(m.duration), expected)


def test_multidimensional_inputs():
    mag, dist_rup = np.meshgrid(
        [5.5, 6.5, 7.5], [5.0, 20.0, 80.0, 150.0], indexing="ij"
    )
    m = KemptonStewart2006(Scenario(mag=mag, dist_rup=dist_rup, v_s30=300.0))
    assert m._ln_dur.shape == (3, 4, 4)
    assert m.duration.shape == (3, 4)
    assert m.std_err.D_5t95v.shape == (3, 4)
    assert m.duration.D_5t95v[2, 1] == scalar_model(7.5, 20.0, 300.0).duration.D_5t95v


def test_out_of_range_values_warn_once():
    with pytest.warns(UserWarning, match="1 of 3 values") as record:
        KemptonStewart2006(
            Scenario(mag=np.array([4.0, 6.0, 8.0]), dist_rup=10.0, v_s30=300.0)
        )
    assert len(record) == 2  # One warning each for the values below and above
