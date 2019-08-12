#!/usr/bin/env python
# encoding: utf-8
"""Gülerce and Abrahamson (2011, :cite:`gulerce11`) V/H model."""

from __future__ import division

import numpy as np

from . import model
from .abrahamson_silva_kamai_2014 import AbrahamsonSilvaKamai2014 as ASK14

__author__ = 'Albert Kottke'


class GulerceAbrahamson2011(model.Model):
    """Gülerce and Abrahamson (2011, :cite:`gulerce11`) V/H model.

    This model was developed for active tectonic regions.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = 'Gülerce & Abrahamson (2011)'
    ABBREV = 'GA11'

    # Reference velocity (m/s)
    V_REF = 1100.

    # Load the coefficients for the model
    COEFF = model.load_data_file(
        'gulerce_abrahamson_2011-coefficients.csv', 2)
    PERIODS = COEFF['period']

    INDEX_PGA = 0
    INDEX_PGV = 1
    INDICES_PSA = 2 + np.arange(27)

    PARAMS = [
        model.NumericParameter('dist_rup', True, None, 200),
        model.NumericParameter('mag', True, 5, 8.5),
        model.NumericParameter('v_s30', True, 450, 1200),
        model.NumericParameter('pga_ref', False, 0.001, 10),
        model.CategoricalParameter('mechanism', True, ['SS', 'NS', 'RS']),
    ]

    def __init__(self, *args, **kwargs):
        """Initialize the model."""
        super(GulerceAbrahamson2011, self).__init__(*args, **kwargs)

        if 'pga_ref' not in self._scenario:
            s = model.Scenario(self._scenario)
            s['v_s30'] = self.V_REF
            m = ASK14(s)
            self._scenario['pga_ref'] = m.pga

        self._ln_ratio = self._calc_ln_ratio()
        self._ln_std = self._calc_ln_std()


    def _calc_ln_ratio(self):
        s = self._scenario
        c = self.COEFF

        dist = np.sqrt(s.dist_rup ** 2 + c.c4 ** 2)

        # Magnitude and distance scaling
        f1 = (
            c.a1 + c.a8 (8.5 - s.mag) ** 2 +
            (c.a2 + c.a3 * (s.mag - c.c1)) * np.log(dist)
        )
        mask = (s.mag < c.c1)
        f1[mask] += c.a4 * (s.mag - c.c1)
        f1[~mask] += c.a5 * (s.mag - c.c1)

        # Site model
        v_1 = np.select(
            [
                self.PERIODS < 0,
                self.PERIODS <= 0.50,
                self.PERIODS <= 1.00,
                self.PERIODS <= 2.00,
                self.PERIODS <= 10.0,
            ],
            [
                862,
                1500.,
                np.exp(8 - 0.795 * np.log(self.PERIODS / 0.21)),
                np.exp(6.76 - 0.297 * np.log(self.PERIODS)),
                700,
            ]
        )
        v_s30_lim = np.minimum(s.v_s30, v_1)

        f5 = c.a10 * np.log(v_s30_lim / c.v_lin)
        mask = v_s30_lim < c.v_lin
        f5[mask] -= (
            -c.b * np.log(s.pga_ref + c.c) +
            c.b * np.log(s.pga_ref + c * (v_s30_lim / c.v_lin) ** c.n)
        )
        f5[~mask] -= (
            c.b * c.n * np.log(v_s30_lim / c.v_lin)
        )

        ln_ratio = f1 + f5
        if s.mechanism == 'RS':
            ln_ratio += c.a6
        elif s.mechanism == 'NS':
            ln_ratio += c.a7

        return ln_ratio

    def _calc_ln_std_within(self):
        c = self.COEFF
        mag = self._scenario.mag
        return np.clip(
            c.s1 + 0.5 * (c.s2 - c.s1) * (mag - 5),
            c.s1, c.s2
        )

    def _calc_ln_std_between(self):
        c = self.COEFF
        mag = self._scenario.mag
        return np.clip(
            c.s3 + 0.5 * (c.s4 - c.s3) * (mag - 5),
            c.s3, c.s4
        )

    def _calc_ln_std(self):
        return np.sqrt(
            self._calc_ln_std_between() ** 2 +
            self._calc_ln_std_within() ** 2
        )
