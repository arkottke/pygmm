#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import os

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm.baker_jayaram_2008 import calc_cond_mean_spectrum
from pygmm.baker_jayaram_2008 import calc_correls

fpath = os.path.join(os.path.dirname(__file__), "data", "baker_jayaram_2008.json")

with open(fpath) as fp:
    data = json.load(fp)


def idfn(case):
    if isinstance(case, dict):
        return "T_2={period_cond}s".format(**case)
    else:
        return case


@pytest.mark.parametrize("case", data["cases_correl"], ids=idfn)
def test_calc_correl(case):
    expected = case["correls"]
    actual = calc_correls(case["periods"], case["period_cond"]).tolist()
    # Points were digitized by hand so need to use larger tolerances
    assert_allclose(actual, expected, atol=0.005, rtol=0.005)


@pytest.mark.parametrize("case", data["cases_cms"], ids=idfn)
@pytest.mark.parametrize("param", ["psas_cms", "ln_stds_cms"])
def test_calc_cond_spectrum(case, param):
    model = data["model"]
    results = calc_cond_mean_spectrum(
        model["periods"],
        np.log(model["psas"]),
        model["ln_stds"],
        case["period_cond"],
        np.log(case["psa_cond"]),
    )
    # Need to go from ln_psa to psa
    actual = np.exp(results[0]) if param == "psas_cms" else results[1]
    expected = case[param]
    assert_allclose(actual, expected, atol=0.0001, rtol=0.0005)
