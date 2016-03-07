#!/usr/bin/env python
# encoding: utf-8

"""Pezeshk, Zandieh, and Tavakoli (2011) ground motion model."""

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class PezeshkZandiehTavakoli2011(model.Model):
    """Pezeshk, Zandieh, and Tavakoli (2011) :cite:`pezeshk11` ground motion
    prediction model.

    Developed for the Eastern North America with a reference velocity of 2000
    m/s.
    """

    NAME = 'Pezeshk et al. (2011)'
    ABBREV = 'Pea11'

    # Reference shear-wave velocity (m/sec)
    V_REF = 2000.

    # Load the coefficients for the model
    COEFF = model.load_data_file('pezeshk_zandieh_tavakoli_2011.csv', 1)
    PERIODS = COEFF['period']

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 23)

    PARAMS = [
        model.NumericParameter('mag', True, 5, 8),
        model.NumericParameter('dist_rup', True, None, 1000),
    ]

    def __init__(self, **kwds):
        """Initialize the model.

        Keyword Args:
            mag (float): moment magnitude of the event (:math:`M_w`)
            dist_rup (float): Closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)
        """
        super(PezeshkZandiehTavakoli2011, self).__init__(**kwds)
        p = self.params

        def calc_log10_resp(period, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
                            c11, c12, c13, c14, sigma_reg):
            # Distance scaling
            del period, c12, c13, c14, sigma_reg

            R = np.sqrt(p['dist_rup'] ** 2 + c11 ** 2)

            return (
                c1 + c2 * p['mag'] + c3 * p['mag'] ** 2 +
                (c4 + c5 * p['mag']) * min(np.log10(R), np.log10(70.)) +
                (c6 + c7 * p['mag']) * max(min(np.log10(R / 70.),
                                               np.log10(140. / 70.)), 0.) +
                (c8 + c9 * p['mag']) * max(np.log10(R / 140.), 0) + c10 * R
            )

        def calc_ln_std(period, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11,
                        c12, c13, c14, sigma_reg):
            # Compute the standard deviation
            del period, c1, c2, c4, c5, c6, c7, c8, c9, c10, c11

            if p['mag'] <= 7.:
                ln_std_mean = c12 * p['mag'] + c13
            else:
                ln_std_mean = -6.95e-3 * p['mag'] + c14

            ln_std = np.sqrt(ln_std_mean ** 2 + sigma_reg ** 2)

            return ln_std

        log10_resp = np.array([calc_log10_resp(*c) for c in self.COEFF])
        self._ln_resp = np.log(np.power(10, log10_resp))

        self._ln_std = np.array([calc_ln_std(*c) for c in self.COEFF])
