#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import os

import pytest

from numpy.testing import assert_allclose

from pygmm import GulerceAbrahamson2011 as GA11


fpath = os.path.join(
    os.path.dirname(__file__), 'data', 'gulerce_abrahamson_2011.json')

with open(fpath) as fp:
    tests = json.load(fp)


def idfn(test):
    if isinstance(test, dict):
        p = test['params']
        return f'M={p["mag"]}, Rjb={["dist_rup"]} km'
    else:
        return test


@pytest.mark.parametrize('test', tests, ids=idfn)
@pytest.mark.parametrize('key', ['ratio', 'ln_std'])
def test_calc_cond_spectrum(test, key):
    m = GA11(test['params'])
    actual = getattr(m, key)
    expected = test['output'][key]
    assert_allclose(actual, expected, atol=0.001, rtol=0.005)
