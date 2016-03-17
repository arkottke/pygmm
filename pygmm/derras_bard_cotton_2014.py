#!/usr/bin/env python3
# encoding: utf-8

from __future__ import division

import json
import os

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class DerrasBardCotton2014(model.Model):
    """Derras, Bard and Cotton (2014, :cite:`derras14`) model.
    """
    NAME = 'Derras, Bard & Cotton (2014)'
    ABBREV = 'DBC13'

    # Load the coefficients for the model
    COEFF = json.load(
        open(
            os.path.join(
                os.path.dirname(__file__),
                'data',
                'derras_bard_cotton_2014.json')
        )
    )
    GRAVITY = 9.80665
    PERIODS = np.array(COEFF['period'])

    INDICES_PSA = np.arange(2, 64)
    INDEX_PGA = 1
    INDEX_PGV = 0
    PARAMS = [
        model.NumericParameter('dist_jb', True, 5, 200),
        model.NumericParameter('mag', True, 4, 7),
        model.NumericParameter('v_s30', True, 200, 800),
        model.NumericParameter('depth_hyp', False, 0, 25),
        model.CategoricalParameter('mechanism', True, ['SS', 'NS', 'RS']),
    ]

    def __init__(self, **kwds):
        """Initialize the model.

        Keyword Args:
            depth_hyp (float): depth of the hypocenter (km).

            dist_jb (float): Joyner-Boore distance to the rupture plane
                (:math:`R_\\text{JB}`, km)

            mag (float): moment magnitude of the event (:math:`M_w`)

            mechanism (str): fault mechanism. Valid options: "SS", "NS",
                "RS".

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

        """
        super(DerrasBardCotton2014, self).__init__(**kwds)
        c = self.COEFF
        p = dict(self.params)

        for k in ['v_s30', 'dist_jb']:
            p['log10_' + k] = np.log10(p[k])
        # Translate to mechanism integer
        p['mechanism'] = dict(NS=1, RS=3, SS=4)[p['mechanism']]

        # Create the normalized parameter matrix
        keys = ['log10_dist_jb', 'mag', 'log10_v_s30', 'depth_hyp',
                'mechanism']
        values = np.array([p[k] for k in keys])
        limits = np.rec.array([c['min_max'][k] for k in keys], names='min,max')
        p_n = np.matrix(
            2 * (values - limits['min']) /
            (limits['max'] - limits['min']) - 1).T

        # Compute the normalized response
        b_1 = np.matrix(c['b_1']).T
        b_2 = np.matrix(c['b_2']).T
        w_1 = np.matrix(c['w_1'])
        w_2 = np.matrix(c['w_2'])

        # Here .A1 returns the flatten matrix
        log10_resp_n = (b_2 + w_2 * np.tanh(b_1 + w_1 * p_n)).A1

        # Convert from normalized values
        log10_resp_limits = np.rec.array(c['min_max']['log10_resp'],
                                         names='min,max')
        log10_resp = (
            0.5 * (log10_resp_n + 1) * (log10_resp_limits['max'] -
                                        log10_resp_limits['min']) +
            log10_resp_limits['min']
        )
        # Convert from m/sec and m/sec² into cm/sec and g
        scale = np.log10(
            np.r_[0.01, self.GRAVITY * np.ones(self.PERIODS.size - 1)])

        self._ln_resp = np.log(10 ** (log10_resp - scale))
        self._ln_std = np.log(
            10 ** np.array(c['log10_std']['total']))
