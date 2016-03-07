#!/usr/bin/env python
# encoding: utf-8

"""Model for the Idriss (2014) ground motion model."""

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class Idriss2014(model.Model):
    """Idriss (2014) :cite:`idriss14` ground motion model.
    """

    NAME = 'Idriss (2014)'
    ABBREV = 'I14'

    # Reference velocity (m/s)
    V_REF = 1200.

    # Load the coefficients for the model
    COEFF = dict(
        small=model.load_data_file('idriss_2014-small.csv', 2),
        large=model.load_data_file('idriss_2014-large.csv', 2), )
    PERIODS = COEFF['small']['period']

    INDEX_PGA = 0
    INDICES_PSA = np.arange(22)

    PARAMS = [
        model.NumericParameter(
            'dist_rup', True, None, 150),
        model.NumericParameter(
            'mag', True, 5, None),
        model.NumericParameter(
            'v_s30', True, 450, 1200),
        model.CategoricalParameter(
            'mechanism', True, ['SS', 'RS'], 'SS'),
    ]

    def __init__(self, **kwds):
        """Initialize the ground motion model.

        Keyword Args:
            dist_rup (float): Closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)

            mag (float): moment magnitude of the event (:math:`M_w`)

            mechanism (str): fault mechanism. *ss* and *ns* mechanism are
                treated the same with :math:`F=0` in the model.

                +------+--------------+-----------------------------+
                | name | meaning      | definition                  |
                +======+==============+=============================+
                | ss   | strike slip  |                             |
                +------+--------------+-----------------------------+
                | ns   | normal slip  | -120° :math:`\le`           |
                |      |              | rake angle :math:`\le` -60° |
                |      |              |                             |
                |      |              | (excludes normal-oblique)   |
                +------+--------------+-----------------------------+
                | rs   | reverse slip | 30° :math:`\le`             |
                |      |              | rake angle :math:`\le` 150° |
                |      |              |                             |
                |      |              | (combines reverse           |
                |      |              | and reverse-oblique)        |
                +------+--------------+-----------------------------+

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).
        """
        super(Idriss2014, self).__init__(**kwds)
        p = self.params

        if p['mechanism'] == 'RS':
            flag_mech = 1
        else:
            # SS/RS/U
            flag_mech = 0

        def calc_ln_resp(period, alpha_1, alpha_2, alpha_3, beta_1, beta_2,
                         epsilon, gamma, phi):
            # Equation 3 on page 1166
            del period
            f_mag = (
                alpha_1 + alpha_2 * p['mag'] + alpha_3 * (8.5 - p['mag']) ** 2)
            f_dst = (-(beta_1 + beta_2 * p['mag']) * np.log(
                p['dist_rup'] + 10) + gamma * p['dist_rup'])
            f_ste = epsilon * np.log(p['v_s30'])
            f_mec = phi * flag_mech

            return f_mag + f_dst + f_ste + f_mec

        coeff = self.COEFF['small'] if p['mag'] <= 6.75 else self.COEFF['large']
        self._ln_resp = np.array([calc_ln_resp(*c) for c in coeff])

        # Equation 4 on page 1168
        self._ln_std = (
            1.18 + 0.035 *
            np.log(np.clip(self.PERIODS, 0.05, 3.0)) -
            0.06 * np.clip(p['mag'], 5.0, 7.5)
        )
