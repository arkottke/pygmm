#!/usr/bin/python
# -*- coding: utf-8 -*-
import gzip
import json
import os

import numpy as np
import pytest

from pygmm import DerrasBardCotton2014 as DBC14
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 1e-2

# Load the tests
fname = os.path.join(os.path.dirname(__file__), "data", "dbc14_tests.json.gz")
with gzip.open(fname, "rt") as fp:
    tests = json.load(fp)

testdata = [(t["params"], t["results"]) for t in tests]


@pytest.mark.parametrize("params,expected", testdata)
def test_spec_accels(params, expected):
    m = DBC14(Scenario(**params))
    np.testing.assert_allclose(
        m.interp_spec_accels(expected["periods"]),
        # Need to convert from m/sec to g
        np.array(expected["spec_accels"]) / m.GRAVITY,
        rtol=RTOL,
        err_msg="Spectral accelerations",
    )


@pytest.mark.parametrize("params,expected", testdata)
@pytest.mark.parametrize("key", ["pga", "pgv"])
def test_im_values(params, expected, key):
    m = DBC14(Scenario(**params))
    # PGA needs to be converted from m/sec² to g, and PGV need to be
    # converted from m/sec into cm/sec
    scale = m.GRAVITY if key == "pga" else 0.01

    np.testing.assert_allclose(
        getattr(m, key),
        expected[key] / scale,
        rtol=RTOL,
    )


@pytest.fixture
def model():
    """Instance of the DBC13 model."""
    return DBC14(Scenario(dist_jb=10, mag=6, v_s30=600, depth_hyp=10, mechanism="SS"))


def test_ln_std(model):
    # Log10 total standard deviations from the paper
    expected = [
        0.298,
        0.309,
        0.31,
        0.313,
        0.319,
        0.323,
        0.325,
        0.328,
        0.335,
        0.338,
        0.338,
        0.337,
        0.338,
        0.337,
        0.335,
        0.334,
        0.333,
        0.33,
        0.328,
        0.328,
        0.326,
        0.325,
        0.322,
        0.322,
        0.326,
        0.328,
        0.329,
        0.33,
        0.33,
        0.329,
        0.327,
        0.328,
        0.328,
        0.327,
        0.33,
        0.331,
        0.332,
        0.332,
        0.332,
        0.331,
        0.33,
        0.331,
        0.333,
        0.335,
        0.339,
        0.343,
        0.346,
        0.353,
        0.356,
        0.359,
        0.362,
        0.365,
        0.368,
        0.368,
        0.37,
        0.37,
        0.373,
        0.375,
        0.377,
        0.378,
        0.378,
        0.376,
        0.375,
        0.375,
    ]

    np.testing.assert_allclose(np.log10(np.exp(model._ln_std)), expected, rtol=RTOL)
