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
    """Boore, Stewart, Seyhan, and Atkinson (2014, :cite:`boore14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    The BSSA14 model defines the following distance attenuation models:

        +--------------+-------------------------------+
        | Name         | Description                   |
        +--------------+-------------------------------+
        | global       | Global; California and Taiwan |
        +--------------+-------------------------------+
        | china_turkey | China and Turkey              |
        +--------------+-------------------------------+
        | italy_japan  | Italy and Japan               |
        +--------------+-------------------------------+

    and the following basin region models:

        +--------+---------------------+
        | Name   | Description         |
        +========+=====================+
        | global | Global / California |
        +--------+---------------------+
        | japan  | Japan               |
        +--------+---------------------+

    These are simplified into one regional parameter with the following
    possibilities:

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
        """Initialize the model.

        Keyword Args:
            mag (float): moment magnitude of the event (:math:`M_w`)

            depth_1_0 (Optional[float]): depth to the 1.0 km∕s shear-wave
                velocity horizon beneath the site, :math:`Z_{1.0}` in (km).
                If *None* is specified, then no adjustment is applied.

            dist_jb (float): Joyner-Boore distance to the rupture plane
                (:math:`R_\\text{JB}`, km)

            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).

            mechanism (str): fault mechanism. Valid options: "U", "SS", "NS",
                "RS"

            region (Optional[str]): region for distance attenuation and basin
                model.  Valid options: "global", "california", "china",
                "italy", "japan", "new_zealand", "taiwan", "turkey".  If
                *None* is specified, then "global" is used as default.
        """
        super(BooreStewartSeyhanAtkinson2014, self).__init__(**kwds)
        pga_ref = np.exp(self._calc_ln_resp(np.nan)[self.INDEX_PGA])
        self._ln_resp = self._calc_ln_resp(pga_ref)
        self._ln_std = self._calc_ln_std()

    def _check_inputs(self):
        super(BooreStewartSeyhanAtkinson2014, self)._check_inputs()

        # Mechanism specific limits
        if self.params['mechanism'] == 'SS':
            _min, _max = 3., 8.5
            if not (_min <= self.params['mag'] <= _max):
                logging.warning(
                    'Magnitude (%g) exceeds recommended bounds (%g to %g)'
                    ' for a strike-slip earthquake!',
                    self.params['mag'], _min, _max
                )
        elif self.params['mechanism'] == 'NS':
            _min, _max = 3., 7.0
            if not (_min <= self.params['mag'] <= _max):
                logging.warning(
                    'Magnitude (%g) exceeds recommended bounds (%g to %g)'
                    ' for a normal-slip earthquake!',
                    self.params['mag'], _min, _max
                )

    def _calc_ln_resp(self, pga_ref):
        """Calculate the natural logarithm of the response.

        Args:
            pga_ref (float): peak ground acceleration (g) at the reference
                condition. If :class:`np.nan`, then no site term is applied.

        Returns:
            :class:`np.array`: Natural log of the response.
        """
        p = self.params
        c = self.COEFF

        # Compute the event term
        ########################
        if p['mechanism'] == 'SS':
            event = np.array(c.e_1)
        elif p['mechanism'] == 'NS':
            event = np.array(c.e_2)
        elif p['mechanism'] == 'RS':
            event = np.array(c.e_3)
        else:
            # Unspecified
            event = np.array(c.e_0)

        mask = p['mag'] <= c.M_h
        event[mask] += (
            c.e_4 * (p['mag'] - c.M_h) +
            c.e_5 * (p['mag'] - c.M_h) ** 2
        )[mask]
        event[~mask] += (c.e_6 * (p['mag'] - c.M_h))[~mask]

        # Compute the distance terms
        ############################
        if p['region'] in ['china', 'turkey']:
            dc_3 = c.dc_3ct
        elif p['region'] in ['italy', 'japan']:
            dc_3 = c.dc_3ij
        else:
            # p['region'] in 'global', 'california', 'new_zealand', 'taiwan'
            dc_3 = c.dc_3global

        dist = np.sqrt(p['dist_jb'] ** 2 + c.h ** 2)
        path = (
            (c.c_1 +
             c.c_2 * (p['mag'] - c.M_ref)) * np.log(dist / c.R_ref) +
            (c.c_3 + dc_3) * (dist - c.R_ref)
        )

        if np.isnan(pga_ref):
            # Reference condition. No site effect
            site = 0
        else:
            # Compute the site term
            f_lin = c.c * np.log(np.minimum(p['v_s30'], c.V_c) / c.V_ref)

            # Add the nonlinearity to the site term
            f_2 = c.f_4 * (np.exp(c.f_5 * (min(p['v_s30'], 760) - 360.)) -
                           np.exp(c.f_5 * (760. - 360.)))
            f_nl = c.f_1 + f_2 * np.log((pga_ref + c.f_3) / c.f_3)

            # Add the basin effect to the site term
            F_dz1 = np.zeros_like(c.period)

            # Compute the average from the Chiou and Youngs (2014)
            # model convert from m to km.
            ln_mz1 = np.log(
                CY14.calc_depth_1_0(p['v_s30'], p['region']))

            if p.get('depth_1_0', None) is not None:
                delta_depth_1_0 = p['depth_1_0'] - np.exp(ln_mz1)
            else:
                delta_depth_1_0 = 0.

            mask = (c.period >= 0.65)
            F_dz1[mask] = np.minimum(c.f_6 * delta_depth_1_0, c.f_7)[mask]

            site = f_lin + f_nl + F_dz1

        ln_resp = event + path + site
        return ln_resp

    def _calc_ln_std(self):
        """Calculate the logarithmic standard deviation.

        Returns:
            :class:`np.array`: Logarithmic standard deviation.
        """
        c = self.COEFF
        p = self.params

        # Uncertainty model
        tau = c.tau_1 + (c.tau_2 - c.tau_1) * \
                        (np.clip(p['mag'], 4.5, 5.5) - 4.5)
        phi = c.phi_1 + (c.phi_2 - c.phi_1) * \
                        (np.clip(p['mag'], 4.5, 5.5) - 4.5)

        # Modify phi for Vs30
        phi -= c.dphi_V * np.clip(np.log(c.V_2 / p['v_s30']) /
                                  np.log(c.V_2 / c.V_1), 0, 1)

        # Modify phi for R
        phi += c.dphi_R * np.clip(np.log(p['dist_jb'] / c.R_1) /
                                  np.log(c.R_2 / c.R_1), 0, 1)

        ln_std = np.sqrt(phi ** 2 + tau ** 2)

        return ln_std
