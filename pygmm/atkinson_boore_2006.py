#!/usr/bin/env python3
# encoding: utf-8

"""Model for the Atkinson and Boore (2006) ground motion model."""

from __future__ import division

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class AtkinsonBoore2006(model.Model):
    """Atkinson and Boore (2006) ground motion prediction model.

    Developed for the Eastern North America with a reference velocity of 760
    or 2000 m/s.

    Citation:
        Atkinson, G. M., & Boore, D. M. (2006). Earthquake ground-motion
        prediction equations for eastern North America. Bulletin of the
        seismological society of America, 96(6), 2181-2205.
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

        Inputs:
            mag: float
                moment magnitude of the event
            dist_rup: float
                closest distance to the rupture (km)
            v_s30: float
                time-averaged shear-wave velocity over the top 30 m of the site (m/s)
        """
        super(AtkinsonBoore2006, self).__init__(**kwds)
        p = self.params

        def compute_log10_resp(period, c_1, c_2, c_3, c_4, c_5, c_6, c_7,
                               c_8, c_9, c_10):
            del period
            R0 = 10.0
            R1 = 70.0
            R2 = 140.0

            f0 = max(np.log10(R0 / p['dist_rup']), 0)
            f1 = min(np.log10(p['dist_rup']), np.log10(R1))
            f2 = max(np.log10(p['dist_rup'] / R2), 0)

            log10_resp = (c_1 +
                          c_2 * p['mag'] +
                          c_3 * p['mag'] ** 2 +
                          (c_4 + c_5 * p['mag']) * f1 +
                          (c_6 + c_7 * p['mag']) * f2 +
                          (c_8 + c_9 * p['mag']) * f0 +
                          c_10 * p['dist_rup']
                          )

            return log10_resp

        def compute_log10_site(pga_bc, period, b_lin, b_1, b_2):
            del period
            VS_1 = 180.
            VS_2 = 300.
            VS_REF = 760.

            if p['v_s30'] <= VS_1:
                b_nl = b_1
            elif VS_1 < p['v_s30'] <= VS_2:
                b_nl = ((b_1 - b_2) * np.log(p['v_s30'] / VS_2) /
                        np.log(VS_1 / VS_2))
            elif VS_2 < p['v_s30'] <= VS_REF:
                b_nl = (b_2 * np.log(p['v_s30'] / VS_REF) /
                        np.log(VS_2 / VS_REF))
            else:
                # Vs30 > VS_REF
                b_nl = 0

            pga_bc = max(pga_bc, 60.)

            log10_site = np.log10(
                np.exp(b_lin * np.log(p['v_s30'] / VS_REF) + b_nl *
                       np.log(pga_bc / 100.)))

            return log10_site

        def compute_log10_stress_factor(stress_drop, period, delta, m1, mh):
            del period
            foo = delta + 0.05
            bar = 0.05 + delta * max(p['mag'] - m1, 0) / (mh - m1)
            log10_stress_factor = min(2., stress_drop / 140.) * min(foo, bar)

            return log10_stress_factor

        COEFF = self.COEFF['bc'] if p['v_s30'] else self.COEFF['rock']

        # Compute the log10 PSA in units of cm/sec/sec
        log10_resp = np.array([compute_log10_resp(*c) for c in COEFF])

        # Apply the stress correction factor as recommended by Atkinson and
        # Boore (2011)
        if p['mag'] >= 5:
            stress_drop = 10. ** (3.45 - 0.2 * p['mag'])
            log10_stress_factor = [compute_log10_stress_factor(stress_drop, *c)
                                   for c in self.COEFF_SF]
            log10_resp += np.interp(COEFF.period,
                                    self.COEFF_SF.period, log10_stress_factor)

        if p['v_s30']:
            # Compute the site amplification
            pga_bc = (10 ** log10_resp[self.INDEX_PGA])

            log10_site = [compute_log10_site(pga_bc, *c)
                          for c in self.COEFF_SITE]

            # Amplification is specified at periods the differ from the ground
            # motion model so we have to interpolate to a common period
            # spacing before adding the influence of the site
            log10_site = np.interp(COEFF.period,
                                   self.COEFF_SITE.period, log10_site)
            log10_resp += log10_site

        # Convert from cm/sec/sec to gravity
        log10_resp -= np.log10(980.665)

        # Convert to log-base 10
        self._ln_resp = np.log(10 ** log10_resp)
        self._ln_std = 0.30
