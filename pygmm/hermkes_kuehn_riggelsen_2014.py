#!/usr/bin/env python3
# encoding: utf-8

from __future__ import division

import os

import numpy as np
from scipy.interpolate import NearestNDInterpolator

from . import model

__author__ = 'Albert Kottke'

fname_data = os.path.join(
    os.path.dirname(__file__),
    'data',
    'hermkes_kuehn_riggelsen_2014.npz'
)

if not os.path.exists(fname_data):
    # Download the model data if not found.
    from six.moves.urllib.request import urlretrieve
    from six.moves.urllib.error import HTTPError

    url = ('https://dl.dropboxusercontent.com/u/1401593/'
           'hermkes_kuehn_riggelsen_2014.npz')
    try:
        urlretrieve(url, fname_data)
    except HTTPError:
        print(
            'Hermkes, Kuehn, and Riggelsen (2013) model data required, '
            'which cannot be downloaded. Download the file from {url}'
            'to this location: {loc}'.format(
                url=url, loc=os.path.abspath(fname_data)
            )
        )
        raise RuntimeError

with np.load(fname_data) as data:
    # Need to transform the record arrays into flat numpy arrays
    INTERPOLATOR = NearestNDInterpolator(
        data['events'], data['predictions'])


class HermkesKuehnRiggelsen2014(model.Model):
    """Hermkes, Kuehn, Riggelsen (2014, :cite:`hermkes14`) model.

    Only the *GPSELinCorr* model is implemented. This model must be imported
    directly by::

        from pygmm.hermkes_kuehn_riggelsen_2014 import
        HermkesKuehnRiggelsen2014

    This is to due to the large file size of the model data, which takes
    time to load.
    """
    NAME = 'Hermkes, Kuehn, Riggelsen (2014)'
    ABBREV = 'HKR14'

    # Reference velocity (m/sec)
    V_REF = None

    PERIODS = np.array([-1, 0.01, 0.1, 0.5, 1.0, 4.0])
    INDICES_PSA = np.arange(1, 6)
    INDEX_PGA = 1
    INDEX_PGV = 0
    PARAMS = [
        model.NumericParameter('depth_hyp', False, 0, 40, 15),
        model.NumericParameter('dist_jb', False, 0, 200),
        model.NumericParameter('mag', True, 4, 8),
        model.NumericParameter('v_s30', True, 100, 1200),
        model.CategoricalParameter('mechanism', True, ['SS', 'NS', 'RS']),
    ]

    def __init__(self, scenario):
        """Initialize the model.

        Note that this model was developed using a Bayesian non-parametric
        method, which means it is should only be used over the data range
        used to develop the model. See the paper for more details.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
        """
        super(HermkesKuehnRiggelsen2014, self).__init__(scenario)

        s = self._scenario

        flag_rs = flag_ss = flag_ns = 0
        if s.mechanism == 'SS':
            flag_ss = 1
        elif s.mechanism == 'NS':
            flag_ns = 1
        elif s.mechanism == 'RS':
            flag_rs = 1

        event = (s.mag, s.depth_hyp, flag_rs, flag_ss, flag_ns, s.dist_jb, s.v_s30)
        prediction = INTERPOLATOR(event)

        self._ln_resp = prediction[0::2]
        self._ln_std = np.sqrt(prediction[1::2])
