#!/usr/bin/env python3
# encoding: utf-8

"""Model for the Campbell (2003) ground motion model."""

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class Campbell2003(model.Model):
    """Campbell (2003, :cite:`campbell03`) model.

    This model was developed for the Eastern US.
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
        """Initialize the model.

        Keyword Args:
            mag (float): moment magnitude of the event (:math:`M_w`)

            dist_rup (float): closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)

        """
        super(Campbell2003, self).__init__(**kwds)
        self._ln_resp = self._calc_ln_resp()
        self._ln_std = self._calc_ln_std()

    def _calc_ln_resp(self):
        """Calculate the natural logarithm of the response.

        Returns:
            :class:`np.array`: Natural logarithm of the response.
        """
        p = self.params
        c = self.COEFF

        f_1 = c.c_2 * p['mag'] + c.c_3 * (8.5 - p['mag']) ** 2

        # Distance scaling
        f_2 = (c.c_4 * np.log(
            np.sqrt(p['dist_rup'] ** 2 +
                    (c.c_7 * np.exp(c.c_8 * p['mag'])) ** 2)) +
               (c.c_5 + c.c_6 * p['mag']) * p['dist_rup'])

        # Geometric attenuation
        r_1 = 70.0
        r_2 = 130.0
        if p['dist_rup'] <= r_1:
            f_3 = 0.
        else:
            f_3 = c.c_9 * (np.log(p['dist_rup']) - np.log(r_1))

            if r_2 < p['dist_rup']:
                f_3 += c.c_10 * (np.log(p['dist_rup']) - np.log(r_2))

        # Compute the ground motion
        ln_resp = c.c_1 + f_1 + f_2 + f_3

        return ln_resp

    def _calc_ln_std(self):
        """Calculate the logarithmic standard deviation.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        p = self.params
        c = self.COEFF

        if p['mag'] < 7.16:
            ln_std = c.c_11 + c.c_12 * p['mag']
        else:
            ln_std = c.c_13

        return ln_std
