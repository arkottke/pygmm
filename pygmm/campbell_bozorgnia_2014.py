#!/usr/bin/env python
# encoding: utf-8

"""Model for the Campbell and Bozorgnia (2014) ground motion model."""

from __future__ import division

import logging

import numpy as np

from . import model
from .chiou_youngs_2014 import ChiouYoungs2014 as CY14

__author__ = 'Albert Kottke'


class CampbellBozorgnia2014(model.Model):
    """Campbell and Bozorgnia (2014, :cite:`campbell14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.
    """

    NAME = 'Campbell & Bozorgnia (2014)'
    ABBREV = 'CB14'

    # Reference velocity (m/sec)
    V_REF = 1100.

    # Load the coefficients for the model
    COEFF = model.load_data_file('campbell_bozorgnia_2014.csv', 2)

    PERIODS = COEFF['period']

    # Period independent model coefficients
    COEFF_C = 1.88
    COEFF_N = 1.18
    COEEF_H_4 = 1

    INDICES_PSA = np.arange(21)
    INDEX_PGA = -2
    INDEX_PGV = -1

    PARAMS = [
        model.NumericParameter('depth_1_0', False),
        model.NumericParameter('depth_2_5', False, 0, 10),
        model.NumericParameter('depth_bor', False),
        model.NumericParameter('depth_bot', False, default=15.),
        model.NumericParameter('depth_hyp', False, 0, 20),
        model.NumericParameter('depth_tor', False, 0, 20),
        model.NumericParameter('dip', True, 15, 90),
        model.NumericParameter('dist_jb', True),
        model.NumericParameter('dist_rup', True, None, 300),
        model.NumericParameter('dist_x', True),
        model.NumericParameter('mag', True, 3.3, 8.5),
        model.NumericParameter('v_s30', True, 150, 1500),
        model.NumericParameter('width', False),

        model.CategoricalParameter(
            'region', False,
            ['global', 'california', 'japan', 'italy', 'china'], 'global'),
        model.CategoricalParameter('mechanism', True, ['SS', 'NS', 'RS']),
    ]

    def _check_inputs(self, **kwds):
        super(CampbellBozorgnia2014, self)._check_inputs(**kwds)
        p = self.params

        for mech, limit in [('SS', 8.5), ('RS', 8.0), ('NS', 7.5)]:
            if mech == p['mechanism'] and p['mag'] > limit:
                logging.warning(
                    'Magnitude of %g is greater than the recommended limit of'
                    '%g for %s style faults',
                    p['mag'], limit, mech
                )

        if p['depth_2_5'] is None:
            p['depth_2_5'] = self.calc_depth_2_5(
                p['v_s30'], p['region'], p['depth_1_0'])

        if p['depth_tor'] is None:
            p['depth_tor'] = CY14.calc_depth_tor(p['mag'], p['mechanism'])

        if p['width'] is None:
            p['width'] = CampbellBozorgnia2014.calc_width(
                p['mag'], p['dip'], p['depth_tor'], p['depth_bot'])

        if p['depth_bor'] is None:
            p['depth_bor'] = self.calc_depth_bor(
                p['depth_tor'], p['dip'], p['width'])

        if p['depth_hyp'] is None:
            p['depth_hyp'] = CampbellBozorgnia2014.calc_depth_hyp(
                p['mag'], p['dip'], p['depth_tor'], p['depth_bor'])

    def __init__(self, **kwds):
        """Compute the response predicted the Campbell and Bozorgnia (2014)
        ground motion model.

        Keyword Args:
            depth_1_0 (Optional[float]): depth to the 1.0 km∕s shear-wave
                velocity horizon beneath the site, :math:`Z_{1.0}` in (km).
                Used to estimate `depth_2_5`.

            depth_2_5 (Optional[float]): depth to the 2.5 km∕s shear-wave
                velocity horizon beneath the site, :math:`Z_{2.5}` in (km).
                If *None*, then it is computed from `depth_1_0` or `v_s30`
                and the `region` parameter.

            depth_tor (Optional[float]): depth to the top of the rupture
                plane (:math:`Z_{tor}`, km). If *None*, then  the average
                model is used.

            depth_bor (Optional[float]): depth to the bottom of the rupture
                plane (:math:`Z_{bor}`, km). If *None*, then  the average
                model is used.

            depth_bot (Optional[float]): depth to bottom of seismogenic crust
                (km). Used to calculate fault width if none is specified. If
                *None*, then a value of 15 km is used.

            depth_hyp (Optional[float]): depth of the hypocenter (km). If
                *None*, then the model average is used.

            dip (float): fault dip angle (:math:`\phi`, deg).

            dist_jb (float): Joyner-Boore distance to the rupture plane
                (:math:`R_\\text{JB}`, km)

            dist_rup (float): closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)

            dist_x (float): site coordinate measured perpendicular to the
                fault strike from the fault line with the down-dip direction
                being positive (:math:`R_x`, km).

            mag (float): moment magnitude of the event (:math:`M_w`)

            mechanism (str): fault mechanism. Valid values: "SS", "NS", "RS".

            region (Optional[str]): region. Valid values: "california",
                "china", "italy", "japan". If *None*, then "california" is
                used as a default value.

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

            width (Optional[float]): Down-dip width of the fault. If *None*,
                then the model average is used.
        """
        super(CampbellBozorgnia2014, self).__init__(**kwds)
        p = self.params

        pga_ref = np.exp(
            self._calc_ln_resp(np.nan, self.V_REF)[self.INDEX_PGA])

        self._ln_resp = self._calc_ln_resp(pga_ref, p['v_s30'])
        self._ln_std = self._calc_ln_std(pga_ref)

    def _calc_ln_resp(self, pga_ref, v_s30):
        """Calculate the natural logarithm of the response.

        Args:
            pga_ref (float): peak ground acceleration (g) at the reference
                condition. If :class:`np.nan`, then no site term is applied.

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

        Returns:
            :class:`np.array`: Natural log of the response.
        """
        p = self.params
        c = self.COEFF

        # Magnitude term
        f_mag = c.c_0 + c.c_1 * p['mag']
        for min_mag, slope in ([4.5, c.c_2], [5.5, c.c_3], [6.5, c.c_4]):
            if min_mag < p['mag']:
                f_mag += slope * (p['mag'] - min_mag)
            else:
                break

        # Geometric attenuation term
        f_dis = (c.c_5 + c.c_6 * p['mag']) * np.log(np.sqrt(
            p['dist_rup'] ** 2 + c.c_7 ** 2
        ))

        # Style of faulting term
        taper = np.clip(p['mag'] - 4.5, 0, 1)
        if p['mechanism'] == 'RS':
            f_flt = c.c_8 * taper
        elif p['mechanism'] == 'NS':
            f_flt = c.c_9 * taper
        else:
            f_flt = 0

        # Hanging-wall term
        R_1 = p['width'] * np.cos(np.radians(p['dip']))
        R_2 = 62 * p['mag'] - 350
        if p['dist_x'] < 0:
            f_hngRx = 0
        elif p['dist_x'] <= R_1:
            ratio = p['dist_x'] / R_1
            f_hngRx = c.h_1 + c.h_2 * ratio + c.h_3 * ratio ** 2
        else:
            ratio = (p['dist_x'] - R_1) / (R_2 - R_1)
            f_hngRx = np.maximum(0, c.h_4 + c.h_5 * ratio + c.h_6 * ratio ** 2)

        if p['dist_rup'] == 0:
            f_hngRrup = 1
        else:
            f_hngRrup = (p['dist_rup'] - p['dist_jb']) / p['dist_rup']

        if p['mag'] <= 5.5:
            f_hngM = 0
        else:
            f_hngM = \
                np.minimum(p['mag'] - 5.5, 1) * (1 + c.a_2 * (p['mag'] - 6.5))

        f_hngZ = 0 if p['depth_tor'] > 16.66 else 1 - 0.06 * p['depth_tor']
        f_hngDip = (90 - p['dip']) / 45

        f_hng = c.c_10 * f_hngRx * f_hngRrup * f_hngM * f_hngZ * f_hngDip

        # Site term
        f_site = np.zeros_like(c.period)
        vs_ratio = v_s30 / c.k_1
        mask = (v_s30 <= c.k_1)
        f_site[mask] = (
            c.c_11 * np.log(vs_ratio) +
            c.k_2 * (np.log(pga_ref +
                            self.COEFF_C * vs_ratio ** self.COEFF_N) -
                     np.log(pga_ref + self.COEFF_C))
            )[mask]
        f_site[~mask] = (
            (c.c_11 + c.k_2 * self.COEFF_N) * np.log(vs_ratio)
        )[~mask]

        if p['region'] == 'japan':
            # Apply regional correction for Japan
            if v_s30 <= 200:
                f_site += (
                    (c.c_12 + c.k_2 * self.COEFF_N) *
                    (np.log(vs_ratio) - np.log(200 / c.k_1))
                )
            else:
                f_site += (c.c_13 + c.k_2 * self.COEFF_N) * np.log(vs_ratio)

        # Basin response term
        if np.isnan(pga_ref):
            # Use model to compute depth_2_5 for the reference velocity case
            depth_2_5 = self.calc_depth_2_5(v_s30, p['region'])
        else:
            depth_2_5 = p['depth_2_5']

        if depth_2_5 <= 1:
            f_sed = c.c_14 * (depth_2_5 - 1)
            if p['region'] == 'japan':
                f_sed += c.c_15 * (depth_2_5 - 1)
        elif depth_2_5 <= 3:
            f_sed = 0
        else:
            f_sed = (c.c_16 * c.k_3 * np.exp(-0.75) *
                     (1 - np.exp(-0.25 * (depth_2_5 - 3))))

        # Hypocentral depth term
        f_hypH = np.clip(p['depth_hyp'] - 7, 0, 13)
        f_hypM = c.c_17 + (c.c_18 - c.c_17) * np.clip(p['mag'] - 5.5, 0, 1)
        f_hyp = f_hypH * f_hypM

        # Fault dip term
        f_dip = c.c_19 * p['dip'] * np.clip(5.5 - p['mag'], 0, 1)

        # Anaelastic attenuation term
        if p['region'] in ['japan', 'italy']:
            dc_20 = c.dc_20jp
        elif p['region'] == ['china']:
            dc_20 = c.dc_20ch
        else:
            dc_20 = c.dc_20ca

        f_atn = (c.c_20 + dc_20) * max(p['dist_rup'] - 80, 0)

        ln_resp = (f_mag + f_dis + f_flt + f_hng + f_site + f_sed + f_hyp +
                   f_dip + f_atn)
        return ln_resp

    def _calc_ln_std(self, pga_ref):
        """Calculate the logarithmic standard deviation.

        Args:
            pga_ref (float): peak ground acceleration (g) at the reference
                condition.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        p = self.params
        c = self.COEFF

        tau_lnY = c.tau_2 + (c.tau_1 - c.tau_2) * np.clip(5.5 - p['mag'], 0, 1)
        phi_lnY = c.phi_2 + (c.phi_1 - c.phi_2) * np.clip(5.5 - p['mag'], 0, 1)

        vs_ratio = p['v_s30'] / c.k_1
        alpha = np.zeros_like(c.period)
        mask = p['v_s30'] < c.k_1
        alpha[mask] = (
            c.k_2 * pga_ref * (
                (pga_ref + self.COEFF_C * vs_ratio ** self.COEFF_N) ** (-1) -
                (pga_ref + self.COEFF_C) ** -1)
        )[mask]

        tau_lnPGA = tau_lnY[self.INDEX_PGA]
        tau = np.sqrt(tau_lnY ** 2 + alpha ** 2 * tau_lnPGA ** 2 +
                      2 * alpha * c.rho_lnPGAlnY * tau_lnY * tau_lnPGA)

        phi_lnPGA = phi_lnY[self.INDEX_PGA]
        phi_lnAF_PGA = self.COEFF['phi_lnAF'][self.INDEX_PGA]
        phi_lnPGA_B = np.sqrt(phi_lnPGA ** 2 - phi_lnAF_PGA ** 2)
        phi_lnY_B = np.sqrt(phi_lnY ** 2 - c.phi_lnAF ** 2)

        phi = np.sqrt(phi_lnY_B ** 2 + c.phi_lnAF ** 2 +
                      alpha ** 2 * (phi_lnPGA ** 2 - phi_lnAF_PGA ** 2) +
                      2 * alpha * c.rho_lnPGAlnY * phi_lnY_B * phi_lnPGA_B)

        ln_std = np.sqrt(phi ** 2 + tau ** 2)

        return ln_std

    @staticmethod
    def calc_depth_2_5(v_s30, region='global', depth_1_0=None):
        """Calculate the depth to a shear-wave velocity of 2.5 km/sec
        (:math:`Z_{2.5}`).

        Provide either `v_s30` or `depth_1_0`.

        Args:
            v_s30 (Optional[float]): time-averaged shear-wave velocity over
                the top 30 m of the site (:math:`V_{s30}`, m/s).

        Keyword Args:
            region (Optional[str]): region of the basin model. Valid values:
                "california", "japan".

            depth_1_0 (Optional[float]): depth to the 1.0 km∕s shear-wave
                velocity horizon beneath the site, :math:`Z_{1.0}` in (km).

        Returns:
            float: estimated depth to a shear-wave velocity of 2.5 km/sec
            (km).
        """
        if v_s30:
            param = v_s30
            if region == 'japan':
                # From equation 6.10 on page 63
                intercept = 5.359
                slope = 1.102
            else:
                # From equation 6.9 on page 63
                intercept = 7.089
                slope = 1.144

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 6.510
            # bar = 1.181
        elif depth_1_0:
            param = depth_1_0
            if region == 'japan':
                # From equation 6.13 on page 64
                intercept = 0.408
                slope = 1.745
            else:
                # From equation 6.12 on page 64
                intercept = 1.392
                slope = 1.798

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 0.748
            # bar = 2.128
        else:
            raise NotImplementedError

        return np.exp(intercept - slope * np.log(param))

    @staticmethod
    def calc_depth_hyp(mag, dip, depth_tor, depth_bor):
        """Estimate the depth to hypocenter.

        Args:
            mag (float): moment magnitude of the event (:math:`M_w`)

            dip (float): fault dip angle (:math:`\phi`, deg).

            depth_tor (float): depth to the top of the rupture
                plane (:math:`Z_{tor}`, km).

            depth_bor (float): depth to the bottom of the rupture
                plane (:math:`Z_{bor}`, km).

        Returns:
            float: estimated hypocenter depth (km)
        """
        # Equations 35, 36, and 37 of journal article
        ln_dZ = min(
            min(-4.317 + 0.984 * mag, 2.325) +
            min(0.0445 * (dip - 40), 0),
            np.log(0.9 * (depth_bor - depth_tor))
        )

        depth_hyp = depth_tor + np.exp(ln_dZ)

        return depth_hyp

    @staticmethod
    def calc_width(mag, dip, depth_tor, depth_bot=15.0):
        """Estimate the fault width using Equation (39) of CB14.

        Args:
            mag (float): moment magnitude of the event (:math:`M_w`)

            dip (float): fault dip angle (:math:`\phi`, deg).

            depth_tor (float): depth to the top of the rupture
                plane (:math:`Z_{tor}`, km).

        Keyword Args:
            depth_bot (Optional[float]): depth to bottom of seismogenic crust
                (km). Used to calculate fault width if none is specified. If
                *None*, then a value of 15 km is used.

        Returns:
            float: estimated fault width (km)
        """
        return min(
            np.sqrt(10 ** ((mag - 4.07) / 0.98)),
            (depth_bot - depth_tor) / np.sin(np.radians(dip))
        )

    @staticmethod
    def calc_depth_bor(depth_tor, dip, width):
        """Compute the depth to bottom of the rupture (km).

        Args:
            dip (float): fault dip angle (:math:`\phi`, deg).

            depth_tor (float): depth to the top of the rupture
                plane (:math:`Z_{tor}`, km).

            width (float): Down-dip width of the fault.

        Returns:
            float: depth to bottom of the fault rupture (km)
        """
        return depth_tor + width * np.sin(np.radians(dip))
