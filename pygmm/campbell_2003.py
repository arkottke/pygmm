#!/usr/bin/env python3
# encoding: utf-8

"""Model for the Campbell (2003) ground motion model."""

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class Campbell2003(model.Model):
    """Campbell (2003) :cite:`campbell03` ground motion model for Eastern US.
    """

    NAME = 'Campbell (2003)'
    ABBREV = 'C03'

    # Reference velocity (m/sec)
    V_REF = 2800.

    COEFF = model.load_data_file('campbell_2003.csv', 1)
    PERIODS = COEFF['period']

    INDICES_PSA = np.arange(16)

    PARAMS = [
        model.NumericParameter('mag', True, 5.0, 8.2),
        model.NumericParameter('dist_rup', True, None, 1000.),
    ]

    def __init__(self, **kwds):
        """Initialize the ground motion model.

        Keyword Args:
            mag (float): moment magnitude of the event (:math:`M_w`)

            dist_rup (float): closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)

        """
        super(Campbell2003, self).__init__(**kwds)
        p = self.params

        def calc_ln_resp(period, c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9,
                         c_10, c_11, c_12, c_13):
            # Magnitude scaling

            f_1 = c_2 * p['mag'] + c_3 * (8.5 - p['mag']) ** 2

            # Distance scaling
            f_2 = (c_4 * np.log(np.sqrt(
                p['dist_rup'] ** 2 + (c_7 * np.exp(c_8 * p['mag'])) ** 2))
                + (c_5 + c_6 * p['mag']) * p['dist_rup'])

            # Geometric attenuation
            r_1 = 70.0
            r_2 = 130.0
            if p['dist_rup'] <= r_1:
                f_3 = 0.
            else:
                f_3 = c_9 * (np.log(p['dist_rup']) - np.log(r_1))

                if r_2 < p['dist_rup']:
                    f_3 += c_10 * (np.log(p['dist_rup']) - np.log(r_2))

            # Compute the ground motion
            ln_resp = c_1 + f_1 + f_2 + f_3

            return ln_resp

        def calc_ln_std(period, c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9,
                        c_10, c_11, c_12, c_13):
            # Compute the standard deviation

            del period, c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10
            if p['mag'] < 7.16:
                ln_std = c_11 + c_12 * p['mag']
            else:
                ln_std = c_13

            return ln_std

        self._ln_resp = np.array([calc_ln_resp(*c) for c in self.COEFF])
        self._ln_std = np.array([calc_ln_std(*c) for c in self.COEFF])
