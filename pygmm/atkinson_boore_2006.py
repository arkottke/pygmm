#!/usr/bin/env python3
# encoding: utf-8

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class AtkinsonBoore2006(model.Model):
    """Atkinson and Boore (2006, :cite:`atkinson06`) model.

    Developed for the Eastern North America with a reference velocity of 760
    or 2000 m/s.
    """

    NAME = 'Atkinson and Boore (2006)'
    ABBREV = 'AB06'

    # Load the coefficients for the model
    COEFF = dict(
        bc=model.load_data_file('atkinson_boore_2006-bc.csv', 2),
        rock=model.load_data_file('atkinson_boore_2006-rock.csv', 2),
    )

    PERIODS = COEFF['bc']['period']

    COEFF_SITE = model.load_data_file('atkinson_boore_2006-site.csv', 2)
    COEFF_SF = model.load_data_file('atkinson_boore_2006-sf.csv', 2)

    INDEX_PGD = 0
    INDEX_PGV = 1
    INDEX_PGA = 2
    INDICES_PSA = np.arange(3, 27)

    PARAMS = [
        model.NumericParameter('mag', True),
        model.NumericParameter('dist_rup', True),
        model.NumericParameter('v_s30', True)
    ]

    def __init__(self, **kwds):
        """Initialize the model.

        Args:
            mag (float): moment magnitude of the event (:math:`M_w`)
            dist_rup (float): Closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)
            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).
        """
        super(AtkinsonBoore2006, self).__init__(**kwds)
        self._ln_resp = self._calc_ln_resp()
        self._ln_std = self._calc_ln_std()

    def _calc_ln_resp(self):
        """Calculate the natural logarithm of the response.

        Returns:
            :class:`np.array`: Natural log of the response.
        """
        p = self.params
        c = self.COEFF['bc'] if p['v_s30'] else self.COEFF['rock']

        # Compute the response at the reference condition
        r0 = 10.0
        r1 = 70.0
        r2 = 140.0

        f0 = np.maximum(np.log10(r0 / p['dist_rup']), 0)
        f1 = np.minimum(np.log10(p['dist_rup']), np.log10(r1))
        f2 = np.maximum(np.log10(p['dist_rup'] / r2), 0)

        # Compute the log10 PSA in units of cm/sec/sec
        log10_resp = (
            c.c_1 +
            c.c_2 * p['mag'] +
            c.c_3 * p['mag'] ** 2 +
            (c.c_4 + c.c_5 * p['mag']) * f1 +
            (c.c_6 + c.c_7 * p['mag']) * f2 +
            (c.c_8 + c.c_9 * p['mag']) * f0 +
            c.c_10 * p['dist_rup']
        )

        # Apply stress drop correction
        log10_resp += self._calc_stress_factor()

        if p['v_s30']:
            # Compute the site amplification
            pga_bc = (10 ** log10_resp[self.INDEX_PGA])

            log10_site = self._calc_log10_site(pga_bc)

            log10_resp += log10_site

        # Convert from cm/sec/sec to gravity
        log10_resp -= np.log10(980.665)

        ln_resp = np.log(10 ** log10_resp)
        return ln_resp

    def _calc_ln_std(self):
        """Calculate the logarithmic standard deviation.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        ln_std = np.ones_like(self.PERIODS) * 0.30
        return ln_std

    def _calc_stress_factor(self):
        """Calculate the stress correction factor proposed by Atkinson and
        Boore (2011) :cite:`atkinson11`.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        p = self.params
        c = self.COEFF_SF

        stress_drop = 10. ** (3.45 - 0.2 * p['mag'])
        v1 = c.delta + 0.05
        v2 = (0.05 + c.delta * np.maximum(p['mag'] - c.m_1, 0) /
              (c.m_h - c.m_1))

        log10_stress_factor = (np.minimum(2., stress_drop / 140.) *
                               np.minimum(v1, v2))

        return np.interp(self.PERIODS, c.period, log10_stress_factor)

    def _calc_log10_site(self, pga_bc):
        """Calculate the log10 of the site amplification.

        Returns:
            :class:`np.array`: Log10 site amplification.
        """
        p = self.params
        c = self.COEFF_SITE
        VS_1 = 180.
        VS_2 = 300.
        VS_REF = 760.

        if p['v_s30'] <= VS_1:
            b_nl = c.b_1
        elif VS_1 < p['v_s30'] <= VS_2:
            b_nl = ((c.b_1 - c.b_2) * np.log(p['v_s30'] / VS_2) /
                    np.log(VS_1 / VS_2))
        elif VS_2 < p['v_s30'] <= VS_REF:
            b_nl = (c.b_2 * np.log(p['v_s30'] / VS_REF) /
                    np.log(VS_2 / VS_REF))
        else:
            # Vs30 > VS_REF
            b_nl = 0

        pga_bc = max(pga_bc, 60.)

        log10_site = np.log10(
            np.exp(c.b_lin * np.log(p['v_s30'] / VS_REF) + b_nl *
                   np.log(pga_bc / 100.)))

        return np.interp(self.PERIODS, c.period, log10_site)
