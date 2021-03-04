#!/usr/bin/python
# -*- coding: utf-8 -*-
import gzip
import json
import os

import numpy as np
import pytest

from pygmm import AkkarSandikkayaBommer2014 as ASB14
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 1e-2

# Load the tests
fname = os.path.join(os.path.dirname(__file__), "data", "asb14_tests.json.gz")
with gzip.open(fname, "rt") as fp:
    tests = json.load(fp)

testdata = [
    (dist, t["params"], results) for t in tests for dist, results in t["results"]
]


def create_model(params, dist):
    params = dict(params)
    params[dist] = params.pop("dist")
    s = Scenario(**params)
    m = ASB14(s)
    return m


@pytest.mark.parametrize("dist,params,expected", testdata)
def test_spec_accels(dist, params, expected):
    m = create_model(params, dist)
    np.testing.assert_allclose(
        m.interp_spec_accels(expected["periods"]),
        expected["spec_accels"],
        rtol=RTOL,
        err_msg="Spectral accelerations",
    )


@pytest.mark.parametrize("dist,params,expected", testdata)
@pytest.mark.parametrize("key", ["pga", "pgv"])
def test_im_values(dist, params, expected, key):
    m = create_model(params, dist)
    np.testing.assert_allclose(
        getattr(m, key),
        expected[key],
        rtol=RTOL,
    )
