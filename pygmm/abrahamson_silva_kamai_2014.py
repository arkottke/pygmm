#!/usr/bin/env python
# encoding: utf-8

"""Abrahamson, Silva, and Kamai (2014) ground motion model."""

from __future__ import division

import numpy as np

from . import model
from .chiou_youngs_2014 import ChiouYoungs2014 as CY14

__author__ = 'Albert Kottke'


class AbrahamsonSilvaKamai2014(model.Model):
    """Abrahamson, Silva, and Kamai (2014) ground motion model.

    Citation:
        Abrahamson, N. A., Silva, W. J., & Kamai, R. (2014). Summary of the
        ASK14 ground motion relation for active crustal regions. Earthquake
        Spectra, 30(3), 1025-1055.
    """

    NAME = 'Abrahamson, Silva, & Kamai (2014)'
    ABBREV = 'ASK14'

    # Reference velocity (m/sec)
    V_REF = 1180.

    # Load the coefficients for the model
    COEFF = model.load_data_file('abrahamson_silva_kamai_2014.csv', 2)

    PERIODS = COEFF['period']

    # Period independent model coefficients
    COEFF_A7 = 0
    COEFF_M2 = 5
    COEFF_N = 1.5

    INDICES_PSA = np.arange(22)
    INDEX_PGA = -2
    INDEX_PGV = -1

    PARAMS = [
        model.NumericParameter('dist_rup', True, None, 300),
        model.NumericParameter('dist_jb', True),
        model.NumericParameter('mag', True, 3, 8.5),
        model.NumericParameter('v_s30', True, 180, 1000),

        model.NumericParameter('depth_1_0', False),
        model.NumericParameter('depth_tor', False),
        model.NumericParameter('dip', True),
        model.NumericParameter('dist_crjb', False, default=15),
        model.NumericParameter('dist_x', False),
        model.NumericParameter('dist_y0', False),
        model.NumericParameter('width', False),

        model.CategoricalParameter('mechanism', True, ['SS', 'NS', 'RS']),
        model.CategoricalParameter(
            'region', False,
            ['global', 'california', 'china', 'italy', 'japan', 'taiwan'],
            'global'
        ),
        model.CategoricalParameter(
            'vs_source', False, ['measured', 'inferred'], 'measured'),
        model.CategoricalParameter(
            'is_aftershock', False, [True, False], False),
        model.CategoricalParameter('on_hanging_wall', False,
                                   [True, False], False),
    ]

    def _check_inputs(self, **kwds):
        super(AbrahamsonSilvaKamai2014, self)._check_inputs(**kwds)

        p = self.params

        if p['width'] is None:
            p['width'] = self.calc_width(p['mag'], p['dip'])

        if p['depth_tor'] is None:
            p['depth_tor'] = self.calc_width(p['mag'])

    def __init__(self, **kwds):
        """
        Inputs:
            depth_1_0: float, default: None
                (optional) depth to the 1.0 km∕s shear wave velocity horizon
                beneath the site. Used to compute `depth_2_5` if it is not
                explicitly provided.
            depth_2_5: float, default: None
                (optional) depth to the 2.5 km∕s shear wave velocity horizon
                beneath the
                site (a.k.a. sediment depth). If `None`, then it is computed
                from `depth_1_0`, or `v_s30`. If `depth_2_5` is to be
                estimated from the `v_s30` and the `region` parameter.
            depth_tor: float, default: None
                (optional) Depth to the top of the rupture (km). If not provided,
                average model is used.
            depth_bor: float, default: None
                (optional) Depth to bottom of the rupture (km).
            depth_hyp: float, default: None
                Depth of the hypocenter (km).
            dip: float
                Fault dip angle (deg)
            dist_jb: float
                Joyner-Boore distance to the rupture plane (km)
            dist_rup: float
                closest distance to the rupture (km)
            dist_x: float
                site coordinate (km) measured perpendicular to the fault strike
                from the fault line with the down-dip direction being positive
                (see Figure 3.12 in Chiou and Youngs (2008).
            dist_y0: float
                (optional) the horizontal distance off the end of the rupture
                measured parallel to strike (km).
            mag: float
                moment magnitude of the event
            mechanism: str
                fault mechanism.

                    SS
                        strike slip

                    NS
                        normal slip
                        -120° <= rake angle <= -60°
                        (excludes normal-oblique)

                    RS
                        reverse slip
                        30° <= rake angle <= 150°
                        (combined reverse and reverse-oblique)
            on_hanging_wall: bool, default: False
                (optional) If the site is located on the hanging wall of the
                fault.
                    True
                        :math:`R_x >= 0`
                    False
                        :math:`R_x < 0`
            region: str, default: global
                region. Potential regions:
                    global
                    california
                    china
                    italy
                    japan
                    taiwan
            v_s30: float
                time-averaged shear-wave velocity over the top 30 m of the site
                (m/s)
            vs_source: str, default: measured
                source of the :math:`V_{s30}`. Potential sources:
                    measured
                    inferred
            width: float, default: None
                (optional) Down-dip width of the fault. If not provided,
                it is estimated.
        """
        super(AbrahamsonSilvaKamai2014, self).__init__(**kwds)
        p = self.params
        flag_nm = 1 if p['mechanism'] == 'NS' else 0
        flag_rv = 1 if p['mechanism'] == 'RS' else 0

        def calc_f1(mag, dist_rup, m1, a1, a2, a3, a4, a5, a6, a8, a17, c4):
            # Magnitude dependent taper
            m2 = self.COEFF_M2
            a7 = self.COEFF_A7
            c4_mag = c4 - (c4 - 1) * np.clip(5 - mag, 0, 1)
            dist = np.sqrt(dist_rup ** 2 + c4_mag ** 2)

            # Magnitude scaling
            base = a1
            if mag <= m2:
                base += (
                    a4 * (m2 - m1) + a8 * (8.5 - m2) ** 2 +
                    a6 * (mag - m2) + a7 * (mag - m2) +
                    (a2 + a3 * (m2 - m1)) * np.log(dist) + a17 * dist_rup
                )
            else:
                if mag <= m1:
                    base += a4 * (mag - m1)
                else:
                    base += a5 * (mag - m1)

                base += (a8 * (8.5 - mag) ** 2 + (a2 + a3 * (mag - m1)) *
                         np.log(dist) + a17 * dist_rup)

            return base

        def calc_f4(dist_jb, dist_x, dist_y0, width, dip, depth_tor, mag, a13):

            t1 = min(90 - dip, 60) / 45

            # Constant from page 1041
            a2hw = 0.2
            if mag <= 5.5:
                t2 = 0
            elif mag < 6.5:
                t2 = 1 + a2hw * (mag - 6.5) - (1 - a2hw) * (mag - 6.5) ** 2
            else:
                t2 = 1 + a2hw * (mag - 6.5)

            # Constants defined on page 1040
            r1 = width * np.cos(np.radians(dip))
            r2 = 3 * r1
            h1 = 0.25
            h2 = 1.5
            h3 = -0.75
            if dist_x < r1:
                t3 = h1 + h2 * (dist_x / r1) + h3 * (dist_x / r1) ** 2
            elif dist_x < r2:
                t3 = 1 - ((dist_x - r1) / (r2 - r1))
            else:
                t3 = 0

            t4 = np.clip(1 - depth_tor ** 2 / 100, 0, 1)

            if dist_y0 is None:
                t5 = np.clip(1 - dist_jb / 30, 0, 1)
            else:
                dist_y1 = dist_x * np.tan(np.radians(20))
                t5 = np.clip(1 - (dist_y0 - dist_y1) / 5, 0, 1)

            return a13 * t1 * t2 * t3 * t4 * t5

        def calc_region(region, v_s30, dist_rup, vs_ratio, a25, a28,
                        a29, a31, a36, a37, a38, a39, a40, a41, a42):
            if region == 'taiwan':
                return a31 * np.log(vs_ratio) + a25 * dist_rup
            elif region == 'china':
                return a28 * dist_rup
            elif region == 'japan':
                f13 = np.interp(
                    v_s30,
                    [150, 250, 350, 450, 600, 850, 1150],
                    [a36, a37, a38, a39, a40, a41, a42],
                    a36, a42
                )
                return f13 + a29 * dist_rup
            else:
                return 0

        def calc_ln_resp(v_s30, psa_ref, period, m1, v_lin, b, c, c4, a1, a2,
                         a3, a4, a5, a6, a8, a10, a11, a12, a13, a14, a15,
                         a17, a43, a44, a45, a46, a25, a28, a29, a31, a36,
                         a37, a38, a39, a40, a41, a42, s1e, s2e, s3, s4, s1m,
                         s2m, s5, s6):
            del s1e, s2e, s3, s4, s1m, s2m, s5, s6

            f1 = calc_f1(p['mag'], p['dist_rup'],
                         m1, a1, a2, a3, a4, a5, a6, a8, a17, c4)

            # Style of faulting
            f7 = a11 * np.clip(p['mag'] - 4, 0, 1) * flag_rv
            f8 = a12 * np.clip(p['mag'] - 4, 0, 1) * flag_nm

            # Site term
            if period <= 0.5:
                v_1 = 1500
            elif period < 3:
                v_1 = np.exp(-0.35 * np.log(period / 0.5) + np.log(1500))
            else:
                v_1 = 800

            vs_ratio = min(v_s30, v_1) / v_lin
            if vs_ratio < 1.:
                f5 = (a10 * np.log(vs_ratio) -
                      b * np.log(psa_ref + c) +
                      b * np.log(psa_ref + c * vs_ratio ** self.COEFF_N)
                      )
            else:
                f5 = (a10 + b * self.COEFF_N) * np.log(vs_ratio)

            if p['on_hanging_wall']:
                # Hanging-wall term
                f4 = calc_f4(p['dist_jb'], p['dist_x'], p['dist_y0'],
                             p['width'], p['dip'], p['depth_tor'], p['mag'],
                             a13)
            else:
                f4 = 0

            # Depth to top of rupture term
            f6 = a15 * np.clip(p['depth_tor'] / 20, 0, 1)

            # Basin term
            if v_s30 == self.V_REF or p['depth_1_0'] is None:
                # No basin response
                f10 = 0
            else:
                depth_1_0_ref = self.calc_depth_1_0(v_s30, p['region'])
                ln_depth_ratio = np.log(
                    (p['depth_1_0'] + 0.01) / (depth_1_0_ref + 0.01))
                slope = np.interp(v_s30,
                                  [150, 250, 400, 700],
                                  [a43, a44, a45, a46],
                                  a43, a46)
                f10 = slope * ln_depth_ratio

            # Aftershock term
            if p['is_aftershock']:
                f11 = a14 * np.clip(1 - (p['dist_crjb'] - 5) / 10, 0, 1)
            else:
                f11 = 0

            region = calc_region(p['region'], v_s30, p['dist_rup'],
                                 vs_ratio, a25, a28, a29, a31, a36, a37,
                                 a38, a39, a40, a41, a42)

            return f1 + f4 + f5 + f6 + f7 + f8 + f10 + f11 + region

        def calc_ln_std(psa_ref, period, m1, v_lin, b, c, c4, a1, a2,
                        a3, a4, a5, a6, a8, a10, a11, a12, a13, a14, a15,
                        a17, a43, a44, a45, a46, a25, a28, a29, a31, a36,
                        a37, a38, a39, a40, a41, a42, s1e, s2e, s3, s4, s1m,
                        s2m, s5, s6):
            del (period, m1, c4, a1, a2, a3, a4, a5, a6, a8, a10, a11, a12,
                 a13, a14, a15, a17, a43, a44, a45, a46, a25, a28, a29, a31,
                 a36, a37, a38, a39, a40, a41, a42)

            if p['region'] == 'japan':
                phi_al = s5 + (s6 - s5) * np.clip((p['dist_rup'] - 30) / 50,
                                                  0, 1)
            else:
                transition = np.clip((p['mag'] - 4) / 2, 0, 1)
                if p['vs_source'] == 'measured':
                    phi_al = s1m + (s2m - s1m) * transition
                else:
                    phi_al = s1e + (s2e - s1e) * transition

            tau_al = s3 + (s4 - s3) * np.clip((p['mag'] - 5) / 2, 0, 1)
            tau_b = tau_al

            # Remove period independent site amplification uncertainty of 0.4
            phi_amp = 0.4
            phi_b = np.sqrt(phi_al ** 2 - phi_amp ** 2)

            # The partial derivative of the amplification with respect to
            # the reference intensity
            if p['v_s30'] >= v_lin:
                deriv = 0
            else:
                deriv = ((-b * psa_ref) / (psa_ref + c) +
                         (b * psa_ref) /
                         (psa_ref + c * (p['v_s30'] / v_lin) ** self.COEFF_N))

            tau = tau_b * (1 + deriv)
            phi = np.sqrt(phi_b ** 2 * (1 + deriv) ** 2 + phi_amp ** 2)

            return np.sqrt(phi ** 2 + tau ** 2)

        psa_ref = np.exp(
            [calc_ln_resp(self.V_REF, None, *c) for c in self.COEFF])

        self._ln_resp = np.array([calc_ln_resp(p['v_s30'], psa_ref_i, *c)
                                  for psa_ref_i, c in zip(psa_ref, self.COEFF)])
        self._ln_std = np.array([calc_ln_std(psa_ref_i, *c)
                                 for psa_ref_i, c in zip(psa_ref, self.COEFF)])

    @staticmethod
    def calc_width(mag, dip):
        """Compute the fault width based on equation in NGW2 spreadsheet.

        Inputs:
            mag: float
                magnitude
            dip: float
                Dip of the fault (degrees).
        Returns:
            width: float
                estimated fault width (km).
        """
        return min(
            18 / np.sin(np.radians(dip)),
            10 ** (-1.75 + 0.45 * mag)
        )

    @staticmethod
    def calc_depth_tor(mag):
        """Calculate the depth to top of rupture (km).

        Inputs:
            mag: float
                magnitude
        Returns:
            depth_tor: float
                Estimated depth (m).
        """
        return np.interp(mag, [5., 7.2], [7.8, 0])

    @staticmethod
    def calc_depth_1_0(v_s30, region):
        """Calculate an estimate of the depth to 1 km/sec, :math:`Z_{1.0}`
        based on :math:`V_{s30}` and region.

        This is based on equations 18 and 19 in the Spectra paper,
        and differs from the equations in the CY14 paper.

        Inputs:
            v_s30: float
                time-averaged shear-wave velocity over the top 30 m of the
                site, :math:`V_{s30}` (m/s)
            region: str
                region for the :math:`V_{s30}-Z_{1.0}` correlation. Potential regions:
                    california (default)
                    japan
        Returns:
            depth_1_0: float
                depth (m) to a shear-wave velocity of 1,000 m/sec,
                :math:`Z_{1.0}`.

        """
        if region in ['japan']:
            # Japan
            power = 2
            v_ref = 412
            slope = -5.23 / power
        else:
            # Global
            power = 4
            v_ref = 610
            slope = -7.67 / power

        return np.exp(slope * np.log((v_s30 ** power + v_ref ** power) /
                                     (1360. ** power + v_ref ** power))) / 1000
