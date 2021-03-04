# -*- coding: utf-8 -*-
"""Derras, Bard and Cotton (2014, :cite:`derras14`) model."""
import json
import os

import numpy as np

from . import model

__author__ = "Albert Kottke"


class DerrasBardCotton2014(model.GroundMotionModel):
    """Derras, Bard and Cotton (2014, :cite:`derras14`) model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Derras, Bard & Cotton (2014)"
    ABBREV = "DBC13"

    # Load the coefficients for the model
    COEFF = json.load(
        open(
            os.path.join(
                os.path.dirname(__file__), "data", "derras_bard_cotton_2014.json"
            )
        )
    )
    GRAVITY = 9.80665
    PERIODS = np.array(COEFF["period"])

    INDICES_PSA = np.arange(2, 64)
    INDEX_PGA = 1
    INDEX_PGV = 0
    PARAMS = [
        model.NumericParameter("dist_jb", True, 5, 200),
        model.NumericParameter("mag", True, 4, 7),
        model.NumericParameter("v_s30", True, 200, 800),
        model.NumericParameter("depth_hyp", True, 0, 25),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)
        c = self.COEFF
        # Values modified during the calculation
        s = dict(self._scenario)

        for k in ["v_s30", "dist_jb"]:
            s["log10_" + k] = np.log10(s[k])
        # Translate to mechanism integer
        s["mechanism"] = dict(NS=1, RS=3, SS=4)[s["mechanism"]]

        # Create the normalized parameter matrix
        keys = ["log10_dist_jb", "mag", "log10_v_s30", "depth_hyp", "mechanism"]
        values = np.array([s[k] for k in keys])
        limits = np.rec.array([c["min_max"][k] for k in keys], names="min,max")
        p_n = np.array(
            2 * (values - limits["min"]) / (limits["max"] - limits["min"]) - 1
        ).T

        # Compute the normalized response
        b_1 = np.array(c["b_1"]).T
        b_2 = np.array(c["b_2"]).T
        w_1 = np.array(c["w_1"])
        w_2 = np.array(c["w_2"])

        log10_resp_n = b_2 + w_2 @ np.tanh(b_1 + w_1 @ p_n)

        # Convert from normalized values
        log10_resp_limits = np.rec.array(c["min_max"]["log10_resp"], names="min,max")
        log10_resp = (
            0.5
            * (log10_resp_n + 1)
            * (log10_resp_limits["max"] - log10_resp_limits["min"])
            + log10_resp_limits["min"]
        )
        # Convert from m/sec and m/sec² into cm/sec and g
        scale = np.log10(np.r_[0.01, self.GRAVITY * np.ones(self.PERIODS.size - 1)])

        self._ln_resp = np.log(10 ** (log10_resp - scale))
        self._ln_std = np.log(10 ** np.array(c["log10_std"]["total"]))
