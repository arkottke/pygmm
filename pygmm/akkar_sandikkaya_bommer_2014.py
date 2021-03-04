# -*- coding: utf-8 -*-
"""Akkar, Sandikkaya, & Bommer (2014, :cite:`akkar14`) model."""
import collections

import numpy as np

from . import model

__author__ = "Albert Kottke"


class AkkarSandikkayaBommer2014(model.GroundMotionModel):
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

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Akkar, Sandikkaya, & Bommer (2014)"
    ABBREV = "ASB14"

    # Reference velocity (m/sec)
    V_REF = 750.0

    # Load the coefficients for the model
    COEFF = collections.OrderedDict(
        (k, model.load_data_file("akkar-sandikkaya-bommer-2014-%s.csv" % k, 2))
        for k in ["dist_jb", "dist_hyp", "dist_epi"]
    )
    PERIODS = np.array(COEFF["dist_jb"].period)

    INDICES_PSA = np.arange(2, 64)
    INDEX_PGA = 0
    INDEX_PGV = 1
    PARAMS = [
        model.NumericParameter("dist_jb", False, 0, 200),
        model.NumericParameter("dist_epi", False, 0, 200),
        model.NumericParameter("dist_hyp", False, 0, 200),
        model.NumericParameter("mag", True, 4, 8),
        model.NumericParameter("v_s30", True, 150, 1200),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super(AkkarSandikkayaBommer2014, self).__init__(scenario)

        s = self._scenario
        for k in self.COEFF:
            if s[k] is not None:
                dist = s[k]
                c = self.COEFF[k]
                break
        else:
            raise NotImplementedError("Must provide at least one distance metric.")

        # Compute the reference response
        ln_resp_ref = (
            c.a_1
            + c.a_3 * (8.5 - s.mag) ** 2
            + (c.a_4 + c.a_5 * (s.mag - c.c_1))
            * np.log(np.sqrt(dist ** 2 + c.a_6 ** 2))
        )
        mask = s.mag <= c.c_1
        ln_resp_ref[mask] += (c.a_2 * (s.mag - c.c_1))[mask]
        ln_resp_ref[~mask] += (c.a_7 * (s.mag - c.c_1))[~mask]

        if s.mechanism == "NS":
            ln_resp_ref += c.a_8
        elif s.mechanism == "RS":
            ln_resp_ref += c.a_9

        pga_ref = np.exp(ln_resp_ref[self.INDEX_PGA])

        # Compute the nonlinear site term
        if s.v_s30 <= self.V_REF:
            vs_ratio = s.v_s30 / self.V_REF
            site = c.b_1 * np.log(vs_ratio) + c.b_2 * np.log(
                (pga_ref + c.c * vs_ratio ** c.n) / ((pga_ref + c.c) * vs_ratio ** c.n)
            )
        else:
            site = c.b_1 * np.log(np.minimum(s.v_s30, c.v_con) / self.V_REF)

        self._ln_resp = ln_resp_ref + site
        self._ln_std = np.array(c.sd_total)
