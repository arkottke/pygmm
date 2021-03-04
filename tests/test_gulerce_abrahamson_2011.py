import numpy as np
import pytest

from . import load_tests
from pygmm import GulerceAbrahamson2011 as GA11


TESTS = load_tests("gulerce_abrahamson_2011.json.gz")


def idfn(test):
    if isinstance(test, dict):
        p = test["params"]
        return f'M={p["mag"]}, Rjb={["dist_rup"]} km'
    else:
        return test


@pytest.mark.parametrize("test", TESTS, ids=idfn)
@pytest.mark.parametrize("key", ["ratio", "ln_std"])
def test_calc_cond_spectrum(test, key):
    m = GA11(test["params"])
    actual = getattr(m, key)
    expected = test["output"][key]
    np.testing.assert_allclose(actual, expected, atol=0.001, rtol=0.005)
