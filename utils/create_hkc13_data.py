#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Move the predictor into a arrays for usage."""

import os

import numpy as np

periods = [
    (-1,  'pgv'),
    (0,   'pga'),
    (0.1, 'T0.1'),
    (1.0, 'T1'),
    (4.0, 'T4'),
]

path = './HermkesAl'

events = np.genfromtxt(
    os.path.join(path, 'predictors.csv'),
    delimiter=',',
    skip_header=1,
    names=['mag', 'depth_hyp', 'flag_rs', 'flag_ss', 'flag_ns', 'dist_jb', 'v_s30'],
)

print(events[:3])
