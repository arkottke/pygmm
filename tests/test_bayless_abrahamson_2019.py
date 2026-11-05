#!/usr/bin/env python
"""Test calculation of CY14 static methods."""

import csv

import numpy as np
import pytest
from numpy.testing import assert_allclose
from scipy.constants import g as GRAVITY

from pygmm import BaylessAbrahamson2019 as BA19

from . import FPATH_DATA

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
    "scenario", iter_scenarios(), ids=lambda s: f"M {s['mag']}, Vs30 {s['v_s30']} m/s"
)
def test_ba19(scenario):
    m = BA19(scenario)
    assert_allclose(m.freqs, scenario["freqs"])
    assert_allclose(m.eas, scenario["eas_median"] * TO_G, rtol=0.005)


def _piecewise(mag, left, right):
    """Eqs. 19a-c written out as separate branches."""
    if mag < 4:
        return left
    if mag <= 6:
        return left + (right - left) / 2 * (mag - 4)
    return right


@pytest.mark.parametrize("mag", [3.0, 4.0, 5.0, 6.0, 7.0])
def test_ba19_ln_std(mag):
    # The first coefficient is larger than the second at some frequencies and
    # smaller at others, so both orderings are covered.
    base = next(iter_scenarios())
    scenario = {k: v for k, v in base.items() if k not in ("freqs", "eas_median")}
    m = BA19({**scenario, "mag": mag})
    c = BA19.COEFF
    tau = _piecewise(mag, c.s1, c.s2)
    phi_s2s = _piecewise(mag, c.s3, c.s4)
    phi_ss = _piecewise(mag, c.s5, c.s6)
    expected = np.sqrt(tau**2 + phi_s2s**2 + phi_ss**2 + c.c1a**2)
    assert np.any(c.s1 > c.s2) and np.any(c.s1 < c.s2)
    assert_allclose(m.ln_std, expected)
