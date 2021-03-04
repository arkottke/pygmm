"""Reference case from AbrahamsonSilva1996 based on spreadsheet by N.
Gregor."""
import json
import os

import pytest
from numpy.testing import assert_allclose

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

    assert_allclose(actual, case["durs"], atol=0.0001, rtol=0.01)
