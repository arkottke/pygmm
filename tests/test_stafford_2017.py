#!/usr/bin/env python
"""Test calculation of CY14 static methods."""

import os

import numpy as np
import pandas as pd

from pygmm import Stafford2017


def test_spec_accels():
    fname = os.path.join(
        os.path.dirname(__file__), "data", "PJS2017_cor_M4pt5_fmin0pt06.csv"
    )
    expected = pd.read_csv(fname, sep=",", header=0)

    cor = Stafford2017.cor(
        expected["Freq"], sigma_E=None, sigma_S=None, sigma_A=None, magnitude=4.5
    )
    np.testing.assert_allclose(
        cor,
        # Need to convert from m/sec to g
        expected["Cor"],
        rtol=0.05,
        err_msg="Correlations",
    )
