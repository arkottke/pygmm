#!/usr/bin/python
# -*- coding: utf-8 -*-

import gzip
import json
import os

import numpy as np
import pytest

from pygmm import DerrasBardCotton2013 as DBC13

# Relative tolerance for all tests
RTOL = 1e-2

# Load the tests
fname = os.path.join(
    os.path.dirname(__file__), 'data', 'dbc13_tests.json.gz')
with gzip.open(fname, 'rt') as fp:
    tests = json.load(fp)

testdata = [(t['params'], t['results']) for t in tests]


@pytest.mark.parametrize('params,expected', testdata)
def test_spec_accels(params, expected):
    m = DBC13(**params)
    np.testing.assert_allclose(
        m.interp_spec_accels(expected['periods']),
        # Need to convert from m/sec to g
        np.array(expected['spec_accels']) / m.GRAVITY,
        rtol=RTOL,
        err_msg='Spectral accelerations'
    )


@pytest.mark.parametrize('params,expected', testdata)
@pytest.mark.parametrize('key', ['pga', 'pgv'])
def test_im_values(params, expected, key):
    m = DBC13(**params)

    # PGA needs to be converted from m/sec² to g, and PGV need to be
    # converted from m/sec into cm/sec
    scale = m.GRAVITY if key == 'pga' else 0.01

    np.testing.assert_allclose(
        getattr(m, key),
        expected[key] / scale,
        rtol=RTOL,
    )
