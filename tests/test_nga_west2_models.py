#!/usr/bin/python
# -*- coding: utf-8 -*-
import gzip
import itertools
import json
import os

import numpy as np
import pytest

import pygmm

models = [
    pygmm.AbrahamsonSilvaKamai2014,
    pygmm.BooreStewartSeyhanAtkinson2014,
    pygmm.CampbellBozorgnia2014,
    pygmm.ChiouYoungs2014,
    pygmm.Idriss2014,
]

# Relative tolerance for all tests
RTOL = 1e-2

# Load the tests
fname = os.path.join(os.path.dirname(__file__), "data", "ngaw2_tests.json.gz")
with gzip.open(fname, "rt") as fp:
    tests = json.load(fp)

testdata = [
    (m, t["params"], t["results"][m.ABBREV])
    for m, t in itertools.product(models, tests)
]


@pytest.mark.parametrize("model,params,expected", testdata)
def test_ln_stds(model, params, expected):
    m = model(pygmm.model.Scenario(**params))
    np.testing.assert_allclose(
        m.interp_ln_stds(expected["periods"]),
        expected["ln_stds"],
        rtol=RTOL,
        err_msg="Logarithmic standard deviations",
    )


@pytest.mark.parametrize("model,params,expected", testdata)
def test_spec_accels(model, params, expected):
    m = model(pygmm.model.Scenario(**params))
    np.testing.assert_allclose(
        m.interp_spec_accels(expected["periods"]),
        expected["spec_accels"],
        rtol=RTOL,
        err_msg="Spectral accelerations",
    )


@pytest.mark.parametrize("model,params,expected", testdata)
@pytest.mark.parametrize("key", ["pga", "ln_std_pga", "pgv", "ln_std_pgv"])
def test_im_values(model, params, expected, key):
    if expected.get(key, None) is None:
        return

    m = model(pygmm.model.Scenario(**params))

    try:
        value = getattr(m, key)
    except NotImplementedError:
        return

    if value is None:
        return

    np.testing.assert_allclose(
        getattr(m, key),
        expected[key],
        rtol=RTOL,
    )
