
import gzip
import json
import os

import numpy as np
import pytest

from pygmm.model import Scenario
from pygmm import AbrahamsonGregorAddo2016 as AGA16

# Relative tolerance for all tests
RTOL = 2e-2

# Load the tests
fname = os.path.join(os.path.dirname(__file__), 'data', 'aga16_tests.json.gz')
with gzip.open(fname, 'rt') as fp:
    tests = json.load(fp)


def create_model(params):
    s = Scenario(**params)
    m = AGA16(s)
    return m


@pytest.mark.parametrize('test', tests)
@pytest.mark.parametrize(
    'key', ['spec_accels', 'ln_stds', 'pga', 'ln_std_pga'])
def test_spec_accels(test, key):
    m = create_model(test['params'])
    np.testing.assert_allclose(
        getattr(m, key),
        test['results'][key],
        rtol=RTOL,
        err_msg=key
    )
