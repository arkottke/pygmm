import numpy as np
import pytest

from . import load_tests
from pygmm import AbrahamsonGregorAddo2016 as AGA16
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 2e-2

TESTS = load_tests("abrahamson_gregor_addo_2016.json.gz")


def create_model(params):
    s = Scenario(**params)
    m = AGA16(s)
    return m


@pytest.mark.parametrize("test", TESTS)
@pytest.mark.parametrize("key", ["spec_accels", "ln_stds", "pga", "ln_std_pga"])
def test_spec_accels(test, key):
    m = create_model(test["params"])
    np.testing.assert_allclose(
        getattr(m, key), test["results"][key], rtol=RTOL, err_msg=key
    )
