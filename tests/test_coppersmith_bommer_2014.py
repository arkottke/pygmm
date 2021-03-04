import numpy as np
import pytest

from . import load_tests
from pygmm import CoppersmithBommer2014 as Model
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 2e-2

TESTS = load_tests("coppersmith_bommer_2014.json.gz")


def create_model(params):
    s = Scenario(**params)
    m = Model(s)
    return m


@pytest.mark.parametrize("test", TESTS)
@pytest.mark.parametrize("key", ["spec_accels", "ln_stds"])
def test_spec_accels(test, key):
    m = create_model(test["params"])
    np.testing.assert_allclose(
        getattr(m, key), test["results"][key], rtol=RTOL, err_msg=key
    )
