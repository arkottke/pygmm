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
    """Model for the Campbell and Bozorgnia (2014) ground motion model.

    Citation:
        Campbell, K. W., & Bozorgnia, Y. (2014). NGA-West2 ground motion
        model for the average horizontal components of PGA, PGV, and 5%
        damped linear acceleration response spectra. Earthquake Spectra,
        30(3), 1087-1115.
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
                    'Magnitude of {mag:g} is greater than the recommended limit ≥'
                    '{:g} for {} style faults'.format(limit, mech, **p)
                )

        if p['depth_2_5'] is None:
            p['depth_2_5'] = self.calc_depth_2_5(p['v_s30'], p['region'], p['depth_1_0'])

        if p['depth_tor'] is None:
            p['depth_tor'] = CY14.calc_depth_tor(p['mag'], p['mechanism'])

        if p['width'] is None:
            p['width'] = CampbellBozorgnia2014.calc_width(
                p['mag'], p['dip'], p['depth_tor'], p['depth_bot'])

        if p['depth_bor'] is None:
            p['depth_bor'] = self.calc_depth_bor(p['depth_tor'], p['dip'], p['width'])

        if p['depth_hyp'] is None:
            p['depth_hyp'] = CampbellBozorgnia2014.calc_depth_hyp(
                p['mag'], p['dip'], p['depth_tor'], p['depth_bor'])


    def __init__(self, **kwds):
        """Compute the response predicted the Campbell and Bozorgnia (2014)
        ground motion model.

        Inputs:
            depth_1_0: float, default: None
                depth to the 1.0 km∕s shear wave velocity horizon beneath the
                site, :math:`Z_{1.0}`. Used to compute :math:`Z_{2.5}`,
                if it is not explicitly provided.
            depth_2_5: float, default: None
                depth to the 2.5 km∕s shear wave velocity horizon beneath the
                site (a.k.a. sediment depth). If `None`, then it is computed
                from `depth_1_0`, or `v_s30`. If `depth_2_5` is to be
                estimated from the `v_s30`.
            depth_tor: float, default: None
                Depth to the top of the rupture (km). If not provided,
                average model is used.
            depth_bor: float, default: None
                (optional) Depth to bottom of the rupture (km).
            depth_bot: float, default: 15.0
                (optional) depth to bottom of seismogenic crust (km). Used to
                calculate fault width if none is specified.
            depth_hyp: float, default: None
                Depth of the hypocenter (km).
            dip: float
                Fault dip angle (deg)
            dist_jb: float
                Joyner-Boore distance to the rupture plane, :math:`R_{JB}` (km)
            dist_rup: float
                closest distance to the rupture (km)
            dist_x: float
                site coordinate (km) measured perpendicular to the fault strike
                from the fault line with the down-dip direction being positive
                (see Figure 3.12 in Chiou and Youngs (2008).
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
            region: str, default: california
                region. Potential regions:
                    california
                    japan_italy
                    eastern_china
            v_s30: float
                time-averaged shear-wave velocity over the top 30 m of the
                site, :math:`V_{s30}` (m/s)
            width: float, default: None
                (optional) Down-dip width of the fault. If not provided,
                it is estimated.
        """
        super(CampbellBozorgnia2014, self).__init__(**kwds)
        p = self.params

        # Japan specific scaling
        S_J = 1 if 'japan' in p['region'] else 0

        flag_nm = flag_rv = 0
        if p['mechanism'] == 'NS':
            flag_nm = 1
        elif p['mechanism'] == 'RS':
            flag_rv = 1

        def calc_ln_resp(pga_ref, v_s30, period, c_0, c_1, c_2, c_3, c_4,
                         c_5, c_6, c_7, c_8, c_9, c_10, c_11, c_12, c_13,
                         c_14, c_15, c_16, c_17, c_18, c_19, k_1, k_2, k_3,
                         a_2, h_1, h_2, h_3, h_4, h_5, h_6, c_20, dc_20ca,
                         dc_20jp, dc_20ch, tau_1, tau_2, phi_1, phi_2,
                         phi_lnaf, phi_c, sigma_s, sigma_arbs, sigma_l,
                         sigma_arbl, rho_lnpgalny):
            # Calculate response
            del (period, tau_1, tau_2, phi_1, phi_2, phi_lnaf, phi_c,
                 sigma_s, sigma_arbs, sigma_l, sigma_arbl, rho_lnpgalny)
            # Magnitude term
            f_mag = c_0 + c_1 * p['mag']
            for min_mag, slope in ([4.5, c_2], [5.5, c_3], [6.5, c_4]):
                if min_mag < p['mag']:
                    f_mag += slope * (p['mag'] - min_mag)
                else:
                    break

            # Geometric attenuation term
            f_dis = (c_5 + c_6 * p['mag']) * np.log(np.sqrt(
                p['dist_rup'] ** 2 + c_7 ** 2
            ))

            # Style of faulting term
            f_fltF = c_8 * flag_rv + c_9 * flag_nm
            f_fltM = np.clip(p['mag'] - 4.5, 0, 1)

            f_flt = f_fltF * f_fltM

            # Hanging-wall term
            R_1 = p['width'] * np.cos(np.radians(p['dip']))
            R_2 = 62 * p['mag'] - 350
            if p['dist_x'] < 0:
                f_hngRx = 0
            elif p['dist_x'] <= R_1:
                ratio = p['dist_x'] / R_1
                f_hngRx = h_1 + h_2 * ratio + h_3 * ratio ** 2
            else:
                ratio = (p['dist_x'] - R_1) / (R_2 - R_1)
                f_hngRx = max(0, h_4 + h_5 * ratio + h_6 * ratio ** 2)

            if p['dist_rup'] == 0:
                f_hngRrup = 1
            else:
                f_hngRrup = (p['dist_rup'] - p['dist_jb']) / p['dist_rup']

            if p['mag'] <= 5.5:
                f_hngM = 0
            else:
                f_hngM = min(p['mag'] - 5.5, 1) * (1 + a_2 * (p['mag'] - 6.5))

            f_hngZ = 0 if p['depth_tor'] > 16.66 else 1 - 0.06 * p['depth_tor']

            f_hngDip = (90 - p['dip']) / 45

            f_hng = c_10 * f_hngRx * f_hngRrup * f_hngM * f_hngZ * f_hngDip

            # Site term
            v_s30_ratio = v_s30 / k_1
            if v_s30 <= k_1:
                f_siteG = c_11 * np.log(v_s30_ratio) + k_2 * (
                    np.log(pga_ref + self.COEFF_C *
                           v_s30_ratio ** self.COEFF_N) -
                    np.log(pga_ref + self.COEFF_C)
                )
            else:
                f_siteG = (c_11 + k_2 * self.COEFF_N) * np.log(v_s30_ratio)

            if v_s30 <= 200:
                f_siteJ = (c_12 + k_2 * self.COEFF_N) * (np.log(v_s30_ratio) -
                                                         np.log(200 / k_1))
            else:
                f_siteJ = (c_13 + k_2 * self.COEFF_N) * np.log(v_s30_ratio)

            f_site = f_siteG + S_J * f_siteJ

            # Basin response term
            if pga_ref is None:
                # Use model to compute depth_2_5 for the reference velocity case
                depth_2_5 = self.calc_depth_2_5(v_s30, p['region'])
            else:
                depth_2_5 = p['depth_2_5']

            if depth_2_5 <= 1:
                f_sed = (c_14 + c_15 * S_J) * (depth_2_5 - 1)
            elif depth_2_5 <= 3:
                f_sed = 0
            else:
                f_sed = (c_16 * k_3 * np.exp(-0.75) *
                         (1 - np.exp(-0.25 * (depth_2_5 - 3))))

            # Hypocentral depth term
            f_hypH = np.clip(p['depth_hyp'] - 7, 0, 13)
            f_hypM = c_17 + (c_18 - c_17) * np.clip(p['mag'] - 5.5, 0, 1)

            f_hyp = f_hypH * f_hypM

            # Fault dip term
            f_dip = c_19 * p['dip'] * np.clip(5.5 - p['mag'], 0, 1)

            # Anaelastic attenuation term
            if p['region'] in ['japan', 'italy']:
                dc_20 = dc_20jp
            elif p['region'] == ['china']:
                dc_20 = dc_20ch
            else:
                dc_20 = dc_20ca

            f_atn = (c_20 + dc_20) * max(p['dist_rup'] - 80, 0)

            ln_resp = (f_mag + f_dis + f_flt + f_hng + f_site + f_sed + f_hyp +
                       f_dip + f_atn)
            return ln_resp

        def calc_tau_lnY(tau_1, tau_2):
            # Equation 27
            return tau_2 + (tau_1 - tau_2) * np.clip(5.5 - p['mag'], 0, 1)

        def calc_phi_lnY(phi_1, phi_2):
            # Equation 28
            return phi_2 + (phi_1 - phi_2) * np.clip(5.5 - p['mag'], 0, 1)

        def calc_ln_std(pga_ref, tau_lnPGA, phi_lnPGA, period, c_0, c_1, c_2,
                        c_3, c_4, c_5, c_6, c_7, c_8, c_9, c_10, c_11, c_12,
                        c_13, c_14, c_15, c_16, c_17, c_18, c_19, k_1, k_2, k_3,
                        a_2, h_1, h_2, h_3, h_4, h_5, h_6, c_20, dc_20ca,
                        dc_20jp, dc_20ch, tau_1, tau_2, phi_1, phi_2,
                        phi_lnAF, phi_c, sigma_s, sigma_arbs, sigma_l,
                        sigma_arbl, rho_lnPGAlny):
            # Compute standard deviation
            del (period, c_0, c_1, c_2, c_3, c_4, c_5, c_6, c_7, c_8, c_9,
                 c_10, c_11, c_12, c_13, c_14, c_15, c_16, c_17, c_18, c_19,
                 k_3, a_2, h_1, h_2, h_3, h_4, h_5, h_6, c_20, dc_20ca,
                 dc_20jp, dc_20ch, phi_c, sigma_s, sigma_arbs, sigma_l,
                 sigma_arbl)
            tau_lnY = calc_tau_lnY(tau_1, tau_2)
            phi_lnY = calc_phi_lnY(phi_1, phi_2)

            v_s30_ratio = p['v_s30'] / k_1
            if p['v_s30'] < k_1:
                alpha = k_2 * pga_ref * (
                    (pga_ref +
                     self.COEFF_C * v_s30_ratio ** self.COEFF_N) ** (-1) -
                    (pga_ref + self.COEFF_C) ** -1)
            else:
                alpha = 0

            tau = np.sqrt(tau_lnY ** 2 + alpha ** 2 * tau_lnPGA ** 2 +
                          2 * alpha * rho_lnPGAlny * tau_lnY * tau_lnPGA)

            phi_lnAF_PGA = self.COEFF['phi_lnaf'][self.INDEX_PGA]
            phi_lnPGA_B = np.sqrt(phi_lnPGA ** 2 - phi_lnAF_PGA ** 2)
            phi_lnY_B = np.sqrt(phi_lnY ** 2 - phi_lnAF ** 2)

            phi = np.sqrt(phi_lnY_B ** 2 + phi_lnAF ** 2 +
                          alpha ** 2 * (phi_lnPGA ** 2 - phi_lnAF_PGA ** 2) +
                          2 * alpha * rho_lnPGAlny * phi_lnY_B * phi_lnPGA_B)

            ln_std = np.sqrt(phi ** 2 + tau ** 2)

            return ln_std

        pga_ref = np.exp(
            calc_ln_resp(None, self.V_REF, *self.COEFF[self.INDEX_PGA]))
        self._ln_resp = np.array(
            [calc_ln_resp(pga_ref, p['v_s30'], *c) for c in self.COEFF])

        tau_lnPGA = calc_tau_lnY(self.COEFF.tau_1[self.INDEX_PGA],
                                 self.COEFF.tau_2[self.INDEX_PGA])
        phi_lnPGA = calc_phi_lnY(self.COEFF.phi_1[self.INDEX_PGA],
                                 self.COEFF.phi_2[self.INDEX_PGA])

        self._ln_std = np.array(
            [calc_ln_std(pga_ref, tau_lnPGA, phi_lnPGA, *c)
             for c in self.COEFF])

    @staticmethod
    def calc_depth_2_5(v_s30, region='global', depth_1_0=None):
        """Calculate the depth to a shear-wave velocity of 2.5 km/sec.

        Inputs:
            v_s30: float
                time-averaged shear-wave velocity over the top 30 m of the site
                (m/s)
            region: str, default: 'global'
                region for the :math:`V_s30`-:math:`Z_{2.5}` correlation. Potential regions:
                    global:     Data from California and Japan
                    california: Data only from California
                    japan:      Data only from Japan
            depth_1_0: float, default: None
                depth (m) to a shearw-wave velocity of 1,000 (m/sec). Only
                used if :math:`V_{s30}` is not provided.

        Returns:
            depth_2_5: float
                Estimated depth (km) to a shear-wave velocity of 2,500 (m/sec)
        """
        if v_s30:
            param = v_s30
            if region == 'japan':
                # From equation 6.10 on page 63
                foo = 5.359
                bar = 1.102
            else:
                # From equation 6.9 on page 63
                foo = 7.089
                bar = 1.144

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 6.510
            # bar = 1.181
        elif depth_1_0:
            param = depth_1_0
            if region == 'japan':
                # From equation 6.13 on page 64
                foo = 0.408
                bar = 1.745
            else:
                # From equation 6.12 on page 64
                foo = 1.392
                bar = 1.798

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 0.748
            # bar = 2.128
        else:
            raise NotImplementedError

        return np.exp(foo - bar * np.log(param))

    @staticmethod
    def calc_depth_hyp(mag, dip, depth_tor, depth_bor):
        """Calculate the depth to hypocenter.

        Inputs:
            mag: float
                magnitude
            dip: float
                Dip of the fault (degrees).
            depth_tor: float
                Depth to top of rupture (km).
            depth_bor: float
                Depth to bottom of seismogenic crust (km).

        Returns:
            depth_hyp: float
                Estimated hypocenter depth (m).
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
        """Compute the fault width based on Equation (39) of CB14.

        Inputs:
            mag: float
                magnitude
            dip: float
                Dip of the fault (degrees).
            depth_tor: float
                Depth to top of rupture (km).
            depth_bot: float, default: 15.0
                Depth to bottom of seismogenic crust (km).

        Returns:
            width: float
                estimated fault width (km).
        """
        return min(
            np.sqrt(10 ** ((mag - 4.07) / 0.98)),
            (depth_bot - depth_tor) / np.sin(np.radians(dip))
        )

    @staticmethod
    def calc_depth_bor(depth_tor, dip, width=None):
        """Compute the depth to bottom of the rupture (km).

        Inputs:
            depth_tor: float
                Depth to top of rupture (km).
            dip: float
                Dip of the fault (degrees).
            width: float
                Fault width (km).

        Returns:
            depth_bor: float
                estimated depth to bottom of the fault rupture (km).
        """
        return depth_tor + width * np.sin(np.radians(dip))
