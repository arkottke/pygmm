import json
import os

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm.kishida_2017 import calc_cond_mean_spectrum_vector

fpath = os.path.join(os.path.dirname(__file__), "data", "kishida_2017.json")
with open(fpath) as fp:
    data = json.load(fp)


def idfn(case):
    if isinstance(case, dict):
        return str(case["periods_cond"])


@pytest.mark.parametrize("case", data["cases"], ids=idfn)
@pytest.mark.parametrize("param", ["psas_cmsv", "ln_stds_cmsv"])
def test_calc_cond_mean_spectrum_vector(case, param):
    model = data["model"]
    ln_psas_cond = np.ma.array(
        np.log(model["psas_target"]),
        mask=(~np.in1d(model["periods"], case["periods_cond"])),
    )
    results = calc_cond_mean_spectrum_vector(
        model["periods"], np.log(model["psas"]), model["ln_stds"], ln_psas_cond
    )
    # Need to go from ln_psa to psa
    actual = np.exp(results[0]) if param == "psas_cmsv" else results[1]
    expected = case[param]
    assert_allclose(actual, expected, atol=0.0001, rtol=0.0005)
