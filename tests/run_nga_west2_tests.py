#!/usr/bin/python
# -*- coding: utf-8 -*-

import gzip

import json
import os

import numpy as np

from ..abrahamson_silva_kamai_2014 import AbrahamsonSilvaKamai2014
from ..boore_stewart_seyhan_atkinson_2014 import BooreStewartSeyhanAtkinson2014
from ..campbell_bozorgnia_2014 import CampbellBozorgnia2014
from ..chiou_youngs_2014 import ChiouYoungs2014
from ..idriss_2014 import Idriss2014

models = [
    AbrahamsonSilvaKamai2014,
    BooreStewartSeyhanAtkinson2014,
    CampbellBozorgnia2014,
    ChiouYoungs2014,
    Idriss2014,
]

# Number of decimal places to test against.
DECIMAL = 4

# Load the tests
fname = os.path.join(
    os.path.dirname(__file__), 'data', 'ngaw2_tests.json.gz')
with gzip.open(fname, 'rt') as fp:
    tests = json.load(fp)


def test_generator():
    for t in tests:
        for model in models:
            yield check_model, model, t['params'], t['results'][model.ABBREV]


def check_model(model, params, values):
    m = model(**params)

    np.testing.assert_array_almost_equal(
        m.interp_spec_accels(values['periods']),
        values['spec_accels'],
        decimal=DECIMAL,
        err_msg='Spectral accelerations'
    )

    np.testing.assert_array_almost_equal(
        m.interp_ln_stds(values['periods']),
        values['ln_stds'],
        decimal=DECIMAL,
        err_msg='Logarithmic standard deviations'
    )

    for key in ['pga', 'pga_ln_std', 'pgv', 'pgv_ln_std']:
        if values[key] is None:
            continue

        if not hasattr(m, key) or getattr(m, key) is None:
            continue

        np.testing.assert_almost_equal(
            getattr(m, key),
            values[key],
            decimal=DECIMAL,
        )
