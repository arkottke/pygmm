#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test calculation of CY14 static methods."""
import csv

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import g as GRAVITY

from . import FPATH_DATA
from pygmm import BaylessAbrahamson2019 as BA19

# Convert from cm/s to g-s
TO_G = 1 / (GRAVITY * 100)


def read_scenarios(fpath):
    with fpath.open() as fp:
        reader = csv.reader(fp)
        rows = list(reader)

    # Read the parameters
    stop = 6
    parser = {"mechanism": str}
    scenarios = [
        {r[0]: parser.get(r[0], float)(r[i + 1]) for r in rows[:stop]}
        for i in range(len(rows[0]) - 1)
    ]

    # Read the spectrum
    start = 8
    freqs = np.array([float(r[0]) for r in rows[start:]])
    for i, s in enumerate(scenarios):
        s["freqs"] = freqs
        s["eas_median"] = np.array([float(r[i + 1]) for r in rows[start:]])

    return scenarios


def iter_scenarios():
    for fpath in FPATH_DATA.glob("BA19-*.csv"):
        yield from read_scenarios(fpath)


@pytest.mark.parametrize(
    "scenario", iter_scenarios(), ids=lambda s: f'M {s["mag"]}, Vs30 {s["v_s30"]} m/s'
)
def test_ba19(scenario):
    m = BA19(scenario)
    assert_allclose(m.freqs, scenario["freqs"])
    assert_allclose(m.eas, scenario["eas_median"] * TO_G, rtol=0.005)
