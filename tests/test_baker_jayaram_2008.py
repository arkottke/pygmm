#!/usr/bin/env python
# -*- coding: utf-8 -*-

import json
import pathlib

import pytest
import numpy as np

from numpy.testing import assert_almost_equal

from pygmm.baker_jayaram_2008 import calc_correl, calc_cond_mean_spectrum


# def load_cases_correl():
#     fpath = pathlib.Path(__file__).parent.joinpath(
#         'data', 'baker_jayaram_2008-figure_4.csv')
#
#     cases = []
#     with fpath.open() as fp:
#         lines = list(fp)
#         case = None
#         for line in lines:
#             if line[0] == '#':
#                 period_cond = float(re.search('T_2=(\S+)s', line).group(1))
#                 case = dict(period_cond=period_cond, periods=[], correls=[])
#             elif line.isspace():
#                 cases.append(case)
#                 continue
#             else:
#                 parts = line.split(', ')
#                 for i, key in enumerate(['periods', 'correls']):
#                     case[key].append(float(parts[i]))
#     return cases

fpath = pathlib.Path(__file__).parent.joinpath(
    'data', 'baker_jayaram_2008.json')

with fpath.open() as fp:
    data = json.load(fp)


def idfn(case):
    return 'T_2={period_cond}s'.format(**case)


@pytest.mark.parameterize('case', data['cases_correl'], ids=idfn)
def test_calc_correl(case):
    expected = case['correls']
    actual = calc_correl(case['periods'], case['period_cond']).tolist()
    # Points were digitized by hand so need to use larger tolerances
    assert_almost_equal(actual, expected, atol=0.005, rtol=0.005)


@pytest.mark.parameterize('case', data['cases_cms'])
@pytest.mark.parameterize('key,attr', [('psas_cms', 0), ('ln_stds_cms', 1)])
def test_calc_cond_spectrum(case, key, attr):
    results = calc_cond_mean_spectrum(
        data['periods'],
        np.log(data['psas']),
        data['ln_stds'],
        case['period_cond'],
        np.log(case['psa_cond']),
    )
    # Need to go from ln_psa to psa
    results = list(results)
    results[0] = np.exp(results[0])

    assert_almost_equal(
        results[attr],
        case[key],
        atol=0.0001,
        rtol=0.0005,
    )
