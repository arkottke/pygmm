#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import pytest

from pygmm.hermkes_kuehn_riggelsen_2014 import HermkesKuehnRiggelsen2014 as HKR13
from pygmm.model import Scenario

# Relative tolerance for all tests
RTOL = 1e-2

names = ("mag", "depth_hyp", "mechanism", "dist_jb", "v_s30")

# These test cases are extracted from the HKC13 data using the
# create_hkc14_data.py script
events = [
    [4.0, 5.0, "RS", 0.0, 100.0],  # Line 2 - first
    [8.0, 5.0, "SS", 29.0, 100.0],  # Line 9472
    [7.3, 15.0, "SS", 16.0, 900.0],  # Line 1245082
    [8.0, 30.0, "NS", 200.0, 760.0],  # Line 1780057 - last
]

# Predictions. The predictions are provided in pairs of predicted variable
# and variance for PGV, PGA, PSA(T=0.1s), PSA(T=0.5s), PSA(T=1.0s), and PSA(
# T=4.0s).
predictions = [
    [
        -4.07083934,
        0.67692656,
        -1.07231844,
        0.69332252,
        -0.33641352,
        0.7994422,
        -0.97409228,
        0.74814538,
        -1.82468142,
        0.92448446,
        -5.12087126,
        1.02068773,
    ],
    [
        -0.52099869,
        0.72063588,
        1.33454867,
        0.72938761,
        1.74372217,
        0.83268984,
        2.00780052,
        0.79246556,
        1.86559946,
        0.96531185,
        0.91670043,
        1.03366518,
    ],
    [
        -1.80883591,
        0.55683512,
        0.3117025,
        0.57986321,
        0.96330332,
        0.67941943,
        0.67846516,
        0.61702095,
        0.13766521,
        0.71312718,
        -0.60224743,
        0.90132807,
    ],
    [
        -3.50927882,
        0.96678294,
        -1.96656956,
        0.91876949,
        -1.75490823,
        0.97815117,
        -1.34961607,
        1.08108927,
        -0.96706897,
        1.30808437,
        -1.82250014,
        1.10876257,
    ],
]


@pytest.mark.slow
@pytest.mark.parametrize("event,prediction", zip(events, predictions))
@pytest.mark.parametrize(
    "indices,attr",
    [
        (0, "pgv"),
        (1, "ln_std_pgv"),
        (2, "pga"),
        (3, "ln_std_pga"),
        ([2, 4, 6, 8, 10], "spec_accels"),
        ([3, 5, 7, 9, 11], "ln_stds"),
    ],
)
def test_model(event, prediction, indices, attr):
    m = HKR13(Scenario(**dict(zip(names, event))))

    prediction = np.asarray(prediction)
    if "ln_" in attr:
        # Convert from variance to standard deviation
        prediction = np.sqrt(prediction[indices])
    else:
        # Convert from log to natural space.
        prediction = np.exp(prediction[indices])

    np.testing.assert_allclose(
        getattr(m, attr),
        prediction,
        rtol=RTOL,
    )
