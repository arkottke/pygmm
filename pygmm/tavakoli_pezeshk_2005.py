#!/usr/bin/env python
# encoding: utf-8

"""Tavakoli and Pezeshk (2005) ground motion model."""

from __future__ import division

import logging

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class TavakoliPezeshk05(model.Model):
    """Tavakoli and Pezeshk (2005) :cite:`tavakoli05` ground motion prediction
    model.

    Developed for the Eastern North America with a reference velocity of 2880
    m/s.
    """

    NAME = 'Tavakoli and Pezeshk (2005)'
    ABBREV = 'TP05'

    # Reference velocity (m/sec)
    V_REF = 2880.

    # Load the coefficients for the model
    COEFF = model.load_data_file('tavakoli_pezeshk_2005.csv', 1)
    PERIODS = COEFF['period']

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 14)

    PARAMS = [
        model.NumericParameter('dist_rup', True, None, 1000),
        model.NumericParameter('mag', True, 5.0, 8.2),
    ]

    def __init__(self, **kwds):
        """Initialize the model.

        Keyword Args:
            mag (float): moment magnitude of the event (:math:`M_w`)
            dist_rup (float): Closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)
        """
        super(TavakoliPezeshk05, self).__init__(**kwds)
        p = self.params

        def calc_ln_resp(period, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,
                         c12, c13, c14, c15, c16):
            del period, c14, c15, c16

            # Magnitude scaling
            f1 = c1 + c2 * p['mag'] + c3 * (8.5 - p['mag']) ** 2.5

            # Distance scaling
            f2 = c9 * np.log(p['dist_rup'] + 4.5)

            if p['dist_rup'] > 70:
                f2 += c10 * np.log(p['dist_rup'] / 70.)

            if p['dist_rup'] > 130:
                f2 += c11 * np.log(p['dist_rup'] / 130.)

            # Calculate scaled, magnitude dependent distance R for use when
            # calculating f3
            R = np.sqrt(p['dist_rup'] ** 2 +
                        (c5 * np.exp(c6 * p['mag'] +
                                     c7 * (8.5 - p['mag']) ** 2.5)) ** 2)

            f3 = ((c4 + c13 * p['mag']) * np.log(R) + (c8 + c12 * p['mag']) * R)

            # Compute the ground motion
            ln_resp = f1 + f2 + f3

            return ln_resp

        def calc_ln_std(period, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,
                        c12, c13, c14, c15, c16):
            del period, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13

            # Compute the standard deviation
            if p['mag'] < 7.2:
                ln_std = c14 + c15 * p['mag']
            else:
                ln_std = c16

            return ln_std

        self._ln_resp = np.array(
                [calc_ln_resp(*c) for c in self.COEFF])
        self._ln_std = np.array([calc_ln_std(*c) for c in self.COEFF])
