#!/usr/bin/env python3
# encoding: utf-8

from __future__ import division

import collections

import numpy as np

from . import model

__author__ = 'Albert Kottke'


class AkkarSandikkayaBommer2014(model.Model):
    """Akkar, Sandikkaya, & Bommer (2014, :cite:`akkar14`) model.

    The model is specified for three different distance metrics. However,
    the implementation uses only one distance metric. They are used in
    the following order:

        1. `dist_jb`

        2. `dist_hyp`

        3. `dist_epi`

    This order was selected based on evaluation of the total standard
    deviation. To compute the response for differing metrics, call the
    model multiple times with different keywords.
    """
    NAME = 'Akkar, Sandikkaya, & Bommer (2014)'
    ABBREV = 'ASB14'

    # Reference velocity (m/sec)
    V_REF = 750.

    # Load the coefficients for the model
    COEFF = collections.OrderedDict(
        (k, model.load_data_file('akkar-sandikkaya-bommer-2014-%s.csv' % k, 2))
        for k in ['dist_jb', 'dist_hyp', 'dist_epi'])
    PERIODS = np.array(COEFF['dist_jb'].period)

    INDICES_PSA = np.arange(2, 64)
    INDEX_PGA = 0
    INDEX_PGV = 1
    PARAMS = [
        model.NumericParameter('dist_jb', False, 0, 200),
        model.NumericParameter('dist_epi', False, 0, 200),
        model.NumericParameter('dist_hyp', False, 0, 200),
        model.NumericParameter('mag', True, 4, 8),
        model.NumericParameter('v_s30', True, 150, 1200),
        model.CategoricalParameter('mechanism', True, ['SS', 'NS', 'RS']),
    ]

    def __init__(self, scenario):
        """Initialize the model.

<<<<<<< HEAD
        The model is specified for three different distance metrics. However,
        the implementation uses only one distance metric. They are used in
        the following order:

            1. `dist_jb`

            2. `dist_hyp`

            3. `dist_epi`

        This order was selected based on evaluation of the total standard
        deviation. To compute the response for differing metrics, call the
        model multiple times with different keywords.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
=======
        One distance metric must be provided. If multiple are provided, only
        one will be used.

        Parameters
        ----------
        dist_jb : float, optional
            Joyner-Boore distance to the rupture plane
            (:math:`R_\\text{JB}`, km)
        dist_epi : float, optional
            Epicentral distance to the rupture plane
            (:math:`R_\\text{epi}`, km)
        dist_hyp : float, optional
            Hypocentral distance to the rupture plane
            (:math:`R_\\text{hyp}`, km).
        mag : float
            moment magnitude of the event (:math:`M_w`)
        mechanism : str
            fault mechanism. Valid options: "SS", "NS", "RS",
            and "U". See :ref:`Mechanism` for more information.
        v_s30 : float
            time-averaged shear-wave velocity over the top 30 m
            of the site (:math:`V_{s30}`, m/s).
>>>>>>> 463f156a57779d7fb9def11b795e00bc38ad0dd8
        """
        super(AkkarSandikkayaBommer2014, self).__init__(scenario)

        s = self._scenario
        for k in self.COEFF:
            if s[k] is not None:
                dist = s[k]
                c = self.COEFF[k]
                break
        else:
            raise NotImplementedError(
                "Must provide at least one distance metric.")

        # Compute the reference response
        ln_resp_ref = (
<<<<<<< HEAD
            c.a_1 + c.a_3 * (8.5 - s.mag) ** 2 +
            (c.a_4 + c.a_5 * (s.mag - c.c_1)) *
            np.log(np.sqrt(dist ** 2 + c.a_6 ** 2))
        )
=======
            c.a_1 + c.a_3 * (8.5 - p['mag']) ** 2 +
            (c.a_4 + c.a_5 *
             (p['mag'] - c.c_1)) * np.log(np.sqrt(dist ** 2 + c.a_6 ** 2)))
>>>>>>> 463f156a57779d7fb9def11b795e00bc38ad0dd8

        mask = (s.mag <= c.c_1)
        ln_resp_ref[mask] += (c.a_2 * (s.mag - c.c_1))[mask]
        ln_resp_ref[~mask] += (c.a_7 * (s.mag - c.c_1))[~mask]

        if s.mechanism == 'NS':
            ln_resp_ref += c.a_8
        elif s.mechanism == 'RS':
            ln_resp_ref += c.a_9

        pga_ref = np.exp(ln_resp_ref[self.INDEX_PGA])

        # Compute the nonlinear site term
<<<<<<< HEAD
        if s.v_s30 <= self.V_REF:
            vs_ratio = s.v_s30 / self.V_REF
            site = (c.b_1 * np.log(vs_ratio) +
                    c.b_2 * np.log((pga_ref + c.c * vs_ratio ** c.n) /
                                   ((pga_ref + c.c) * vs_ratio ** c.n))
                    )
=======
        if p['v_s30'] <= self.V_REF:
            vs_ratio = p['v_s30'] / self.V_REF
            site = (c.b_1 * np.log(vs_ratio) + c.b_2 * np.log(
                (pga_ref + c.c * vs_ratio ** c.n) / (
                    (pga_ref + c.c) * vs_ratio ** c.n)))
>>>>>>> 463f156a57779d7fb9def11b795e00bc38ad0dd8
        else:
            site = c.b_1 * np.log(np.minimum(s.v_s30, c.v_con) / self.V_REF)

        self._ln_resp = ln_resp_ref + site
        self._ln_std = np.array(c.sd_total)
