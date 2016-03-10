#!/usr/bin/env python3
# encoding: utf-8

"""Chiou and Youngs (2014) ground motion model."""

from __future__ import division

import logging

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class ChiouYoungs2014(model.Model):
    """Chiou and Youngs (2014) :cite:`chiou14` model.
    """
    NAME = 'Chiou and Youngs (2014)'
    ABBREV = 'CY14'

    # Reference velocity (m/s)
    V_REF = 1130.

    # Load the coefficients for the model
    COEFF = model.load_data_file('chiou_youngs_2014.csv', 2)

    PERIODS = COEFF['period']

    INDICES_PSA = np.arange(24)
    INDEX_PGA = 24
    INDEX_PGV = 25

    PARAMS = [
        model.NumericParameter('dist_rup', True, 0, 300),
        model.NumericParameter('dist_x', True),
        model.NumericParameter('dist_jb', True),
        model.NumericParameter('mag', True, 3.5, 8.5),
        model.NumericParameter('v_s30', True, 180, 1500),
        model.NumericParameter('depth_tor', False, 0, 20),
        model.NumericParameter('depth_1_0', False),
        model.NumericParameter('dpp_centered', False, default=0),
        model.NumericParameter('dip', True),
        model.CategoricalParameter('mechanism', False,
                                   ['U', 'SS', 'NS', 'RS'], 'U'),
        model.CategoricalParameter('on_hanging_wall', False,
                                   [True, False], False),
        model.CategoricalParameter('region', False,
                                   ['california', 'china', 'italy', 'japan'],
                                   'california'),
        model.CategoricalParameter('vs_source', False,
                                   ['measured', 'inferred'], 'measured'),
    ]

    def __init__(self, **kwds):
        """Compute the response predicted the Chiou and Youngs (2014) ground
        motion model.

        Keyword Args:
            depth_1_0 (Optional[float]): depth to the 1.0 km∕s shear-wave
                velocity horizon beneath the site, :math:`Z_{1.0}` in (km).
                Used to estimate `depth_2_5`.

            depth_tor (Optional[float]): depth to the top of the rupture
                plane (:math:`Z_{tor}`, km). If *None*, then  the average
                model is used.

            dist_jb (float): Joyner-Boore distance to the rupture plane
                (:math:`R_\\text{JB}`, km)

            dist_rup (float): closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)

            dist_x (float): site coordinate measured perpendicular to the
                fault strike from the fault line with the down-dip direction
                being positive (:math:`R_x`, km).

            dip (float): fault dip angle (:math:`\phi`, deg).

            dpp_centered (Optional[float]): direct point parameter for
                directivity effect (see Chiou and Spudich (2014)) centered on the
                earthquake-specific average DPP for California. If *None*,
                then value of 0 used to provide the average directivity.

            mag (float): moment magnitude of the event (:math:`M_w`)

            mechanism (str): fault mechanism. Valid options: "U", "SS", "NS", "RS".

            on_hanging_wall (Optional[bool]): If the site is located on the
                hanging wall of the fault. If *None*, then *False* is assumed.

            region (Optional[str]): region. Valid options: "california",
                "china", "italy", "japan". If *None*, then "california" is used as
                a default value.

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

            vs_source (Optional[str]): source of the `v_s30` value.  Valid
                options: "measured", "inferred"

        """
        super(ChiouYoungs2014, self).__init__(**kwds)
        p = self.params

        flag_nm = 1 if p['mechanism'] == 'NS' else 0
        flag_rv = 1 if p['mechanism'] == 'RS' else 0
        flag_meas = 1 if p['vs_source'] == 'measured' else 0
        flag_hw = int(p['on_hanging_wall'])

        # Difference between requested and model centered values
        diff_depth_tor = (p['depth_tor'] -
                          self.calc_depth_tor(p['mag'], p['mechanism']))
        diff_depth_1_0 = 1000 * (p['depth_1_0'] -
                                 self.calc_depth_1_0(p['v_s30'], p['region']))

        def calc_ln_resp_ref(
                period, c_1, c_1a, c_1b, c_1c, c_1d, c_2, c_3, c_4, c_4a, c_5,
                c_6, c_7, c_7b, c_8, c_8a, c_8b, c_9, c_9a, c_9b, c_11, c_11b,
                c_hm, c_m, c_n, c_rb, c_gamma1, c_gamma2, c_gamma3, gamma_ji,
                gamma_c, phi_1, phi_1jp, phi_2, phi_3, phi_4, phi_5, phi_5jp,
                phi_6, phi_6jp, tau_1, tau_2, sigma_1, sigma_2, sigma_2jp,
                sigma_3):
            # Compute the response for the reference velocity

            del (period, phi_1, phi_1jp, phi_2, phi_3, phi_4,
                 phi_5, phi_5jp, phi_6, phi_6jp, tau_1, tau_2, sigma_1, sigma_2,
                 sigma_2jp, sigma_3)

            cosh_mag = np.cosh(2 * max(p['mag'] - 4.5, 0))

            ln_resp = c_1

            # Reverse fault term
            ln_resp += (c_1a + c_1c / cosh_mag) * flag_rv

            # Normal fault term
            ln_resp += (c_1b + c_1d / cosh_mag) * flag_nm

            # Magnitude scaling
            ln_resp += c_2 * (p['mag'] - 6)
            ln_resp += (c_2 - c_3) / c_n * np.log(
                1 + np.exp(c_n * (c_m - p['mag'])))

            # Top of rupture term relative to model average
            ln_resp += (c_7 + c_7b / cosh_mag) * diff_depth_tor

            # Dip angle term
            ln_resp += (c_11 + c_11b / cosh_mag) * np.cos(
                np.radians(p['dip'])) ** 2

            # Distance terms
            ln_resp += c_4 * np.log(
                p['dist_rup'] + c_5 * np.cosh(c_6 * max(p['mag'] - c_hm, 0)))
            ln_resp += (c_4a - c_4) * np.log(
                np.sqrt(p['dist_rup'] ** 2 + c_rb ** 2))

            # Regional adjustment
            if p['region'] in ['japan', 'italy'] and (6 < p['mag'] < 6.9):
                scale = gamma_ji
            elif p['region'] in ['china']:
                scale = gamma_c
            else:
                scale = 1.
            ln_resp += (scale *
                        (c_gamma1 + c_gamma2 / np.cosh(
                            max(p['mag'] - c_gamma3, 0))) * p['dist_rup'])

            # Directivity term
            ln_resp += (c_8 *
                        max(1 - max(p['dist_rup'] - 40, 0) / 30, 0) *
                        min(max(p['mag'] - 5.5, 0) / 0.8, 1) *
                        np.exp(-c_8a * (p['mag'] - c_8b) ** 2) *
                        p['dpp_centered'])

            # Hanging wall term
            ln_resp += (c_9 * flag_hw * np.cos(np.radians(p['dip'])) *
                        (c_9a + (1 - c_9a) * np.tanh(p['dist_x'] / c_9b)) *
                        (1 - np.sqrt(p['dist_jb'] ** 2 + p['depth_tor'] ** 2) /
                         (p['dist_rup'] + 1)))

            return ln_resp

        def calc_ln_resp(
                ln_resp_ref, period, c_1, c_1a, c_1b, c_1c, c_1d, c_2, c_3,
                c_4, c_4a, c_5, c_6, c_7, c_7b, c_8, c_8a, c_8b, c_9, c_9a,
                c_9b, c_11, c_11b, c_hm, c_m, c_n, c_rb, c_gamma1, c_gamma2,
                c_gamma3, gamma_ji, gamma_c, phi_1, phi_1jp, phi_2, phi_3,
                phi_4, phi_5, phi_5jp, phi_6, phi_6jp, tau_1, tau_2,
                sigma_1, sigma_2, sigma_2jp, sigma_3):
            # Calculate response for the site condition

            del (period, c_1, c_1a, c_1b, c_1c, c_1d, c_2, c_3, c_4, c_4a, c_5,
                c_6, c_7, c_7b, c_8, c_8a, c_8b, c_9, c_9a, c_9b, c_11, c_11b,
                c_hm, c_m, c_n, c_rb, c_gamma1, c_gamma2, c_gamma3, gamma_ji,
                gamma_c, tau_1, tau_2, sigma_1, sigma_2, sigma_2jp, sigma_3)

            if p['region'] in ['japan']:
                phi_1 = phi_1jp
                phi_5 = phi_5jp
                phi_6 = phi_6jp

            ln_resp = ln_resp_ref
            ln_resp += phi_1 * min(np.log(p['v_s30'] / 1130.), 0)

            ln_resp += (phi_2 * (np.exp(phi_3 * (min(p['v_s30'], 1130.) - 360.)) -
                                 np.exp(phi_3 * (1130. - 360.))) *
                        np.log((np.exp(ln_resp_ref) + phi_4) / phi_4))

            ln_resp += phi_5 * (1 - np.exp(-diff_depth_1_0 / phi_6))

            return ln_resp

        def calc_ln_std(
                ln_resp_ref, period, c_1, c_1a, c_1b, c_1c, c_1d, c_2, c_3,
                c_4, c_4a, c_5, c_6, c_7, c_7b, c_8, c_8a, c_8b, c_9, c_9a,
                c_9b, c_11, c_11b, c_hm, c_m, c_n, c_rb, c_gamma1, c_gamma2,
                c_gamma3, gamma_ji, gamma_c, phi_1, phi_1jp, phi_2, phi_3,
                phi_4, phi_5, phi_5jp, phi_6, phi_6jp, tau_1, tau_2,
                sigma_1, sigma_2, sigma_2jp, sigma_3):

            # Calculate the standard deviation

            del (period, c_1, c_1a, c_1b, c_1c, c_1d, c_2, c_3, c_4, c_4a, c_5,
                 c_6, c_7, c_7b, c_8, c_8a, c_8b, c_9, c_9a, c_9b, c_11, c_11b,
                 c_hm, c_m, c_n, c_rb, c_gamma1, c_gamma2, c_gamma3, gamma_ji,
                 gamma_c, phi_1, phi_1jp, phi_5, phi_5jp, phi_6, phi_6jp)

            if p['region'] in ['japan']:
                sigma_2 = sigma_2jp

            clipped_mag = np.clip(p['mag'], 5., 6.5) - 5.
            resp_ref = np.exp(ln_resp_ref)

            tau = tau_1 + (tau_2 - tau_1) / 1.5 * clipped_mag

            nl_0 = (phi_2 * (np.exp(phi_3 * (min(p['v_s30'], 1130.) - 360.)) -
                             np.exp(phi_3 * (1130. - 360.))) *
                    (resp_ref / (resp_ref + phi_4)))

            phi_nl = (
                (sigma_1 + (sigma_2 - sigma_1) / 1.5 * clipped_mag) *
                np.sqrt(sigma_3 * (1 - flag_meas) + 0.7 * flag_meas +
                        (1 + nl_0) ** 2))

            ln_std = np.sqrt((1 + nl_0) ** 2 * tau ** 2 + phi_nl ** 2)

            return ln_std

        ln_resp_ref = [calc_ln_resp_ref(*c) for c in self.COEFF]

        self._ln_resp = np.array([calc_ln_resp(lrr, *c)
                                  for (lrr, c) in zip(ln_resp_ref, self.COEFF)])
        self._ln_std = np.array([calc_ln_std(lrr, *c)
                                 for (lrr, c) in zip(ln_resp_ref, self.COEFF)])

    def _check_inputs(self):
        super(ChiouYoungs2014, self)._check_inputs()

        if self.params['mechanism'] in ['RS', 'NS']:
            _min, _max = 3.5, 8.0
        else:
            _min, _max = 3.5, 8.5

        if not (_min <= self.params['mag'] <= _max):
            logging.warning(
                'Magnitude (%g) exceeds recommended bounds (%g to %g)'
                ' for a %s earthquake!',
                    self.params['mag'], _min, _max, self.params['mechanism']
                )

        if self.params.get('depth_tor', None) is None:
            self.params['depth_tor'] = self.calc_depth_tor(
                self.params['mag'], self.params['mechanism'])

        if self.params.get('depth_1_0', None) is None:
            # Calculate depth (m) and convert to (km)
            self.params['depth_1_0'] = self.calc_depth_1_0(
                self.params['v_s30'], self.params['region'])

    @staticmethod
    def calc_depth_1_0(v_s30, region):
        """Calculate an estimate of the depth to 1 km/sec (:math:`Z_{1.0}`)
        based on :math:`V_{s30}` and region.

        Args:
            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

            region (str): basin region. Valid options: "california", "japan"

        Returns:
            float: estimated depth to a shear-wave velocity of 1 km/sec (km)
        """
        if region in ['japan']:
            # Japan
            power = 2
            v_ref = 412.39
            slope = -5.23 / power
        else:
            # Global
            power = 4
            v_ref = 570.94
            slope = -7.15 / power

        return np.exp(slope * np.log((v_s30 ** power + v_ref ** power) /
                                     (1360. ** power + v_ref ** power))) / 1000

    @staticmethod
    def calc_depth_tor(mag, mechanism):
        """Calculate an estimate of the depth to top of rupture (km).

        Args:
            mag (float): moment magnitude of the event (:math:`M_w`)

            mechanism (str): fault mechanism. Valid options: "U", "SS", "NS", "RS".

        Returns:
            float: estimated depth to top of rupture (km)
        """
        if mechanism == 'RS':
            # Reverse and reverse-oblique faulting
            foo = 2.704 - 1.226 * max(mag - 5.849, 0)
        else:
            # Combined strike-slip and normal faulting
            foo = 2.673 - 1.136 * max(mag - 4.970, 0)

        return max(foo, 0) ** 2
