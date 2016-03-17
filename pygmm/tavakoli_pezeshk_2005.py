#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class TavakoliPezeshk05(model.Model):
    """Tavakoli and Pezeshk (2005, :cite:`tavakoli05`) model.

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
        self._ln_resp = self._calc_ln_resp()
        self._ln_std = self._calc_ln_std()

    def _calc_ln_resp(self):
        """Calculate the natural logarithm of the response.

        Returns:
            :class:`np.array`: Natural log of the response.
        """
        p = self.params
        c = self.COEFF

        # Magnitude scaling
        f1 = c.c_1 + c.c_2 * p['mag'] + c.c_3 * (8.5 - p['mag']) ** 2.5

        # Distance scaling
        f2 = c.c_9 * np.log(p['dist_rup'] + 4.5)

        if p['dist_rup'] > 70:
            f2 += c.c_10 * np.log(p['dist_rup'] / 70.)

        if p['dist_rup'] > 130:
            f2 += c.c_11 * np.log(p['dist_rup'] / 130.)

        # Calculate scaled, magnitude dependent distance R for use when
        # calculating f3
        dist = np.sqrt(
            p['dist_rup'] ** 2 +
            (c.c_5 * np.exp(c.c_6 * p['mag'] +
                            c.c_7 * (8.5 - p['mag']) ** 2.5)) ** 2)

        f3 = ((c.c_4 + c.c_13 * p['mag']) * np.log(dist) +
              (c.c_8 + c.c_12 * p['mag']) * dist)

        # Compute the ground motion
        ln_resp = f1 + f2 + f3

        return ln_resp

    def _calc_ln_std(self):
        """Calculate the logarithmic standard deviation.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        p = self.params
        c = self.COEFF

        if p['mag'] < 7.2:
            ln_std = c.c_14 + c.c_15 * p['mag']
        else:
            ln_std = c.c_16

        return ln_std
