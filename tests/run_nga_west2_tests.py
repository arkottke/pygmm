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

# Number of decimal places to test against.
DECIMAL = 4

# Load the tests
fname = os.path.join(
    os.path.dirname(__file__), 'data', 'ngaw2_tests.json.gz')
with gzip.open(fname, 'rt') as fp:
    tests = json.load(fp)

testdata = [(m, t['params'], t['results'][m.ABBREV])
            for m, t in itertools.product(models, tests)]


@pytest.mark.parametrize('model,params,expected', testdata)
def test_model(model, params, expected):
    m = model(**params)

    np.testing.assert_array_almost_equal(
        m.interp_spec_accels(expected['periods']),
        expected['spec_accels'],
        decimal=DECIMAL,
        err_msg='Spectral accelerations'
    )

    np.testing.assert_array_almost_equal(
        m.interp_ln_stds(expected['periods']),
        expected['ln_stds'],
        decimal=DECIMAL,
        err_msg='Logarithmic standard deviations'
    )

    for key in ['pga', 'pga_ln_std', 'pgv', 'pgv_ln_std']:
        if expected[key] is None:
            continue

        if not hasattr(m, key) or getattr(m, key) is None:
            continue

        np.testing.assert_almost_equal(
            getattr(m, key),
            expected[key],
            decimal=DECIMAL,
        )
