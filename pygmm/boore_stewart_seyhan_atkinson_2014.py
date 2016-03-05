#!/usr/bin/env python
# encoding: utf-8

"""Boore, Stewart, Seyhan, and Atkinson (2014) ground motion model."""

from __future__ import division

import logging

import numpy as np

from . import model
from .chiou_youngs_2014 import ChiouYoungs2014 as CY14

__author__ = 'Albert Kottke'


class BooreStewartSeyhanAtkinson2014(model.Model):
    """Boore, Stewart, Seyhan, and Atkinson (2014) ground motion model.

    Developed for the California and other active tectonic environments.

    Citation:
        Boore, D. M., Stewart, J. P., Seyhan, E., & Atkinson, G. M. (2014).
        NGA-West2 Equations for Predicting PGA, PGV, and 5% Damped PSA for
        Shallow Crustal Earthquakes. Earthquake Spectra, 30(3), 1057-1085.
    """
    NAME = 'Boore, Stewart, Seyhan, and Atkinson (2014)'
    ABBREV = 'BSSA14'

    # Reference shear-wave velocity in m/sec
    V_REF = 760.

    # Load the coefficients for the model
    COEFF = model.load_data_file('boore_stewart_seyhan_atkinson-2014.csv', 2)
    PERIODS = COEFF['period']

    INDEX_PGV = 0
    INDEX_PGA = 1
    INDICES_PSA = np.arange(2, 107)

    LIMITS = dict(
        mag=(3.0, 8.5),
        dist_jb=(0., 300.),
        v_s30=(150., 1500.),
    )

    PARAMS = [
        model.NumericParameter('mag', True, 3, 8.5),
        model.NumericParameter('depth_1_0', False),
        model.NumericParameter('dist_jb', True, None, 300.),
        model.NumericParameter('v_s30', True, 150., 1500.),

        model.CategoricalParameter('mechanism', False,
                                   ['U', 'SS', 'NS', 'RS'], 'U'),
        model.CategoricalParameter(
            'region', False,
            ['global', 'california', 'china', 'italy', 'japan', 'new_zealand',
             'taiwan', 'turkey'],
            'global'),
    ]

    def __init__(self, **kwds):
        """Compute the response predicted the Boore, Stewart, Seyhan,
        and Atkinson (2014) ground motion model.

        Inputs:
            dist_jb: float
                Joyner-Boore distance to the rupture plane (km)
            mag: float
                moment magnitude of the event
            mechanism: str, default: 'U'
                (optional) fault mechanism with the following options:

                    +--------------+--------------+
                    | Abbreviation | Name         |
                    +==============+==============+
                    | U            | Unspecified  |
                    +--------------+--------------+
                    | SS           | Strike-slip  |
                    +--------------+--------------+
                    | NS           | Normal slip  |
                    +--------------+--------------+
                    | RS           | Reverse slip |
                    +--------------+--------------+

            region: str, default: 'global'
                (optional) region for distance attenuation and basin model.
                The BSSA14 model defines the following distance attenuation
                models:
                    global:         Global; California and Taiwan
                    china_turkey:   China and Turkey
                    italy_japan:    Italy and Japan
                and the following basin region models:
                    global:         Global / California
                    japan:          Japan
                These are simplified into one regional parameter with the
                following possibilities:

                    +-------------+--------------+------------+
                    | Region      | Attenuation  | Basin      |
                    +=============+==============+============+
                    | global      | global       | global     |
                    +-------------+--------------+------------+
                    | california  | global       | global     |
                    +-------------+--------------+------------+
                    | china       | china_turkey | global     |
                    +-------------+--------------+------------+
                    | italy       | italy_japan  | global     |
                    +-------------+--------------+------------+
                    | japan       | italy_japan  | japan      |
                    +-------------+--------------+------------+
                    | new zealand | italy_japan  | global     |
                    +-------------+--------------+------------+
                    | taiwan      | global       | global     |
                    +-------------+--------------+------------+
                    | turkey      | china_turkey | global     |
                    +-------------+--------------+------------+

                If *None* is specified, then 'global' parameters are used.

            v_s30: float
                time-averaged shear-wave velocity over the top 30 m of the
                site (m/s).
            depth_1_0: float, default: None
                (optional) depth to the 1.0 km∕s shear wave velocity horizon
                beneath the site (:math:`Z_{1.0}`). If *None* is specified,
                then no adjustment is applied.
        """
        super(BooreStewartSeyhanAtkinson2014, self).__init__(**kwds)
        p = self.params

        U = SS = NS = RS = 0
        if p['mechanism'] == 'U':
            U = 1
        elif p['mechanism'] == 'SS':
            SS = 1
        elif p['mechanism'] == 'NS':
            NS = 1
        elif p['mechanism'] == 'RS':
            RS = 1
        else:
            raise NotImplementedError

        def calc_ln_resp(pga_ref, period, e_0, e_1, e_2, e_3, e_4,
                         e_5, e_6, M_h, c_1, c_2, c_3, M_ref, R_ref, h,
                         dc_3global, dc_3ct, dc_3ij, c, V_c, V_ref, f_1, f_3,
                         f_4, f_5, f_6, f_7, R_1, R_2, dphi_R, dphi_V, V_1,
                         V_2, phi_1, phi_2, tau_1, tau_2):
            del (R_1, R_2, dphi_R, dphi_V, V_1, V_2, phi_1, phi_2, tau_1,
                 tau_2)
            # Compute the event term
            event = e_0 * U + e_1 * SS + e_2 * NS + e_3 * RS
            if p['mag'] <= M_h:
                event += e_4 * (p['mag'] - M_h) + e_5 * (p['mag'] - M_h) ** 2
            else:
                event += e_6 * (p['mag'] - M_h)

            # Compute the distance terms
            R = np.sqrt(p['dist_jb'] ** 2 + h ** 2)

            if p['region'] in ['global', 'california', 'new_zealand', 'taiwan']:
                dc_3 = dc_3global
            elif p['region'] in ['china', 'turkey']:
                dc_3 = dc_3ct
            elif p['region'] in ['italy', 'japan']:
                dc_3 = dc_3ij
            else:
                raise NotImplementedError

            path = ((c_1 + c_2 * (p['mag'] - M_ref)) * np.log(R / R_ref) +
                    (c_3 + dc_3) * (R - R_ref))

            if pga_ref is not None:
                # Compute the site term
                f_lin = c * np.log(min(p['v_s30'], V_c) / V_ref)

                # Add the nonlinearity to the site term
                f_2 = f_4 * (np.exp(f_5 * (min(p['v_s30'], 760) - 360.)) -
                             np.exp(f_5 * (760. - 360.)))
                f_nl = f_1 + f_2 * np.log((pga_ref + f_3) / f_3)

                # Add the basin effect to the site term
                if period >= 0.65:
                    # Compute the average from the Chiou and Youngs (2014)
                    # model convert from m to km.
                    ln_mz1 = np.log(
                        CY14.calc_depth_1_0(p['v_s30'], p['region']))

                    if p.get('depth_1_0', None) is not None:
                        dz1 = p['depth_1_0'] - np.exp(ln_mz1)
                    else:
                        dz1 = 0.
                    F_dz1 = min(f_6 * dz1, f_7)
                else:
                    F_dz1 = 0

                site = f_lin + f_nl + F_dz1
            else:
                site = 0

            return event + path + site

        def calc_ln_std(period, e_0, e_1, e_2, e_3, e_4, e_5, e_6, M_h, c_1,
                        c_2, c_3, M_ref, R_ref, h, dc_3global, dc_3ct,
                        dc_3ij, c, V_c, V_ref, f_1, f_3, f_4, f_5, f_6, f_7,
                        R_1, R_2, dphi_R, dphi_V, V_1, V_2, phi_1, phi_2,
                        tau_1, tau_2):
            del (period, e_0, e_1, e_2, e_3, e_4, e_5, e_6, M_h, c_1, c_2,
                 c_3, M_ref, R_ref, h, dc_3global, dc_3ct, dc_3ij, c, V_c,
                 V_ref, f_1, f_3, f_4, f_5, f_6, f_7)
            # Uncertainty model
            tau = tau_1 + (tau_2 - tau_1) * (np.clip(p['mag'], 4.5, 5.5) - 4.5)
            phi = phi_1 + (phi_2 - phi_1) * (np.clip(p['mag'], 4.5, 5.5) - 4.5)

            # Modify phi for Vs30
            phi -= dphi_V * np.clip(np.log(V_2 / p['v_s30']) /
                                    np.log(V_2 / V_1), 0, 1)

            # Modify phi for R
            phi += dphi_R * np.clip(np.log(p['dist_jb'] / R_1) /
                                    np.log(R_2 / R_1), 0, 1)

            return np.sqrt(phi ** 2 + tau ** 2)

        pga_ref = np.exp(
            calc_ln_resp(None, *self.COEFF[self.INDEX_PGA]))
        self._ln_resp = np.array(
            [calc_ln_resp(pga_ref, *c) for c in self.COEFF])
        self._ln_std = np.array(
            [calc_ln_std(*c) for c in self.COEFF])

    def _check_inputs(self):
        super(BooreStewartSeyhanAtkinson2014, self)._check_inputs()

        # Mechanism specific limits
        if self.params['mechanism'] == 'SS':
            _min, _max = 3., 8.5
            if not (_min <= self.params['mag'] <= _max):
                logging.warning(
                    'Magnitude ({}) exceeds recommended bounds ({} to {})'\
                    ' for a strike-slip earthquake!'.format(
                        self.params['mag'], _min, _max))
        elif self.params['mechanism'] == 'NS':
            _min, _max = 3., 7.0
            if not (_min <= self.params['mag'] <= _max):
                logging.warning(
                    'Magnitude ({}) exceeds recommended bounds ({} to {})'\
                    ' for a normal-slip earthquake!'.format(
                        self.params['mag'], _min, _max))
