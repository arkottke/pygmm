#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Move the predictor into a arrays for usage."""
import os

import numpy as np

periods = [
    (-1, "pgv"),
    (0, "pga"),
    (0.1, "T0.1"),
    (0.5, "T0.5"),
    (1.0, "T1"),
    (4.0, "T4"),
]

path = "./HermkesAl"

events = np.genfromtxt(
    os.path.join(path, "predictors.csv"),
    delimiter=",",
    skip_header=1,
    # Do not create records as it complicates using the data in the
    # interpolation
    # names=['mag', 'depth_hyp', 'flag_rs', 'flag_ss', 'flag_ns', 'dist_jb',
    # 'v_s30'],
)

fname = "gpselin_corrnoise_pred.csv"

predictions = []
for period, period_name in periods:
    p = np.genfromtxt(
        os.path.join(path, period_name, fname),
        delimiter=",",
    )
    predictions.append(p)

# Convert to a recarray
# names = [pn + suffix for p, pn in periods for suffix in ['', '_var']]

predictions = np.array([p[:, i] for p in predictions for i in range(2)]).T

np.savez_compressed(
    "../pygmm/data/hermkes_kuehn_riggelsen_2014.npz",
    events=events,
    predictions=predictions,
)

# Print test cases
for i in [0, 9470, 1245080, -1]:
    print("Event", i, "\n", events[i])
    print("Prediction", "\n", predictions[i])
