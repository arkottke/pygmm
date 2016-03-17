#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division

import logging

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class ChiouYoungs2014(model.Model):
    """Chiou and Youngs (2014, :cite:`chiou14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.
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
        """Initialize the model.

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
                directivity effect (see Chiou and Spudich (2014)) centered on
                the earthquake-specific average DPP for California. If *None*,
                then value of 0 used to provide the average directivity.

            mag (float): moment magnitude of the event (:math:`M_w`)

            mechanism (str): fault mechanism. Valid options: "U", "SS", "NS",
                "RS".

            on_hanging_wall (Optional[bool]): If the site is located on the
                hanging wall of the fault. If *None*, then *False* is assumed.

            region (Optional[str]): region. Valid options: "california",
                "china", "italy", "japan". If *None*, then "california" is used
                as a default value.

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

            vs_source (Optional[str]): source of the `v_s30` value.  Valid
                options: "measured", "inferred"

        """
        super(ChiouYoungs2014, self).__init__(**kwds)
        ln_resp_ref = self._calc_ln_resp_ref()
        self._ln_resp = self._calc_ln_resp_site(ln_resp_ref)
        self._ln_std = self._calc_ln_std(np.exp(ln_resp_ref))

    def _calc_ln_resp_ref(self):
        """Calculate the natural logarithm of the response at the reference
        site condition.

        Returns:
            :class:`np.array`: Natural logarithm of the response.
        """
        c = self.COEFF
        p = self.params

        cosh_mag = np.cosh(2 * max(p['mag'] - 4.5, 0))
        ln_resp = np.array(c.c_1)

        if p['mechanism'] == 'RS':
            # Reverse fault term
            ln_resp += (c.c_1a + c.c_1c / cosh_mag)
        elif p['mechanism'] == 'NS':
            # Normal fault term
            ln_resp += (c.c_1b + c.c_1d / cosh_mag)

        # Magnitude scaling
        ln_resp += c.c_2 * (p['mag'] - 6)
        ln_resp += (c.c_2 - c.c_3) / c.c_n * np.log(
            1 + np.exp(c.c_n * (c.c_m - p['mag'])))

        # Top of rupture term relative to model average
        diff_depth_tor = (p['depth_tor'] -
                          self.calc_depth_tor(p['mag'], p['mechanism']))
        ln_resp += (c.c_7 + c.c_7b / cosh_mag) * diff_depth_tor

        # Dip angle term
        ln_resp += (c.c_11 + c.c_11b / cosh_mag) * np.cos(
            np.radians(p['dip'])) ** 2

        # Distance terms
        ln_resp += c.c_4 * np.log(
            p['dist_rup'] +
            c.c_5 * np.cosh(c.c_6 * np.maximum(p['mag'] - c.c_hm, 0)))
        ln_resp += (c.c_4a - c.c_4) * np.log(
            np.sqrt(p['dist_rup'] ** 2 + c.c_rb ** 2))

        # Regional adjustment
        if p['region'] in ['japan', 'italy'] and (6 < p['mag'] < 6.9):
            scale = c.gamma_ji
        elif p['region'] in ['china']:
            scale = c.gamma_c
        else:
            scale = 1.
        ln_resp += (scale *
                    (c.c_gamma1 + c.c_gamma2 / np.cosh(
                        np.maximum(p['mag'] - c.c_gamma3, 0))) * p['dist_rup'])

        # Directivity term
        ln_resp += (c.c_8 *
                    max(1 - max(p['dist_rup'] - 40, 0) / 30, 0) *
                    min(max(p['mag'] - 5.5, 0) / 0.8, 1) *
                    np.exp(-c.c_8a * (p['mag'] - c.c_8b) ** 2) *
                    p['dpp_centered'])

        # Hanging wall term
        if p['on_hanging_wall']:
            ln_resp += (
                c.c_9 * np.cos(np.radians(p['dip'])) *
                (c.c_9a + (1 - c.c_9a) * np.tanh(p['dist_x'] / c.c_9b)) *
                (1 - np.sqrt(p['dist_jb'] ** 2 + p['depth_tor'] ** 2) /
                 (p['dist_rup'] + 1))
            )

        return ln_resp

    def _calc_ln_resp_site(self, ln_resp_ref):
        """Calculate the natural logarithm of the response including site
        effects.

        Args:
            ln_resp_ref (:class:`np.array`): Natural logarithm of the
            response at the reference site condition at each of the periods
            specified by the model coefficients.

        Returns:
            :class:`np.array`: Natural log of the response.
        """
        c = self.COEFF
        p = self.params
        if p['region'] in ['japan']:
            phi_1 = c.phi_1jp
            phi_5 = c.phi_5jp
            phi_6 = c.phi_6jp
        else:
            phi_1 = c.phi_1
            phi_5 = c.phi_5
            phi_6 = c.phi_6

        ln_resp = np.array(ln_resp_ref)
        ln_resp += phi_1 * min(np.log(p['v_s30'] / 1130.), 0)

        ln_resp += (
            c.phi_2 * (np.exp(c.phi_3 * (min(p['v_s30'], 1130.) - 360.)) -
                       np.exp(c.phi_3 * (1130. - 360.))) *
            np.log((np.exp(ln_resp_ref) + c.phi_4) / c.phi_4)
        )

        diff_depth_1_0 = 1000 * (p['depth_1_0'] -
                                 self.calc_depth_1_0(p['v_s30'], p['region']))
        ln_resp += phi_5 * (1 - np.exp(-diff_depth_1_0 / phi_6))

        return ln_resp

    def _calc_ln_std(self, resp_ref):
        """Calculate the logarithmic standard deviation.

        Args:
            resp_ref (:class:`np.array`): Response at the reference site
            condition at each of the periods specified by the model
            coefficients.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        c = self.COEFF
        p = self.params

        if p['region'] in ['japan']:
            sigma_2 = c.sigma_2jp
        else:
            sigma_2 = c.sigma_2

        clipped_mag = np.clip(p['mag'], 5., 6.5) - 5.
        tau = c.tau_1 + (c.tau_2 - c.tau_1) / 1.5 * clipped_mag

        nl_0 = (c.phi_2 * (np.exp(c.phi_3 * (min(p['v_s30'], 1130.) - 360.)) -
                           np.exp(c.phi_3 * (1130. - 360.))) *
                (resp_ref / (resp_ref + c.phi_4)))

        flag_meas = 1 if p['vs_source'] == 'measured' else 0
        phi_nl = (
            (c.sigma_1 + (sigma_2 - c.sigma_1) / 1.5 * clipped_mag) *
            np.sqrt(c.sigma_3 * (1 - flag_meas) + 0.7 * flag_meas +
                    (1 + nl_0) ** 2))

        ln_std = np.sqrt((1 + nl_0) ** 2 * tau ** 2 + phi_nl ** 2)

        return ln_std

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

            mechanism (str): fault mechanism. Valid options: "U", "SS", "NS",
                "RS".

        Returns:
            float: estimated depth to top of rupture (km)
        """
        if mechanism == 'RS':
            # Reverse and reverse-oblique faulting
            depth_tor = 2.704 - 1.226 * max(mag - 5.849, 0)
        else:
            # Combined strike-slip and normal faulting
            depth_tor = 2.673 - 1.136 * max(mag - 4.970, 0)

        return max(depth_tor, 0) ** 2
