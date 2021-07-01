"""Gülerce and Abrahamson (2011, :cite:`gulerce11`) V/H model."""
from __future__ import division

import numpy as np

from . import model
from .abrahamson_silva_kamai_2014 import AbrahamsonSilvaKamai2014 as ASK14

__author__ = "Albert Kottke"


class GulerceAbrahamson2011(model.Model):
    """Gülerce and Abrahamson (2011, :cite:`gulerce11`) V/H model.

    This model was developed for active tectonic regions.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Gülerce & Abrahamson (2011)"
    ABBREV = "GA11"

    # Reference velocity (m/s)
    V_REF = 1100.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("gulerce_abrahamson_2011-coefficients.csv", 2)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDEX_PGV = 1
    INDICES_PSA = 2 + np.arange(27)

    PARAMS = [
        model.NumericParameter("dist_rup", True, None, 200),
        model.NumericParameter("mag", True, 5, 8.5),
        model.NumericParameter("v_s30", True, 450, 1200),
        model.NumericParameter("pga_ref", True, 0.001, 10),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def __init__(self, *args, **kwargs):
        """Initialize the model."""
        super(GulerceAbrahamson2011, self).__init__(*args, **kwargs)

        if "pga_ref" not in self._scenario:
            s = model.Scenario(self._scenario)
            s["v_s30"] = self.V_REF
            m = ASK14(s)
            self._scenario["pga_ref"] = m.pga

        self._ln_ratio = self._calc_ln_ratio()
        self._ln_std = self._calc_ln_std()

    @property
    def ratio(self):
        return np.exp(self._ln_ratio)

    @property
    def ln_std(self):
        return self._ln_std

    def _calc_ln_ratio(self):
        s = self._scenario
        c = self.COEFF

        dist = np.sqrt(s.dist_rup ** 2 + c.c4 ** 2)

        # Magnitude and distance scaling
        f1 = (
            c.a1
            + c.a8 * (8.5 - s.mag) ** 2
            + (c.a2 + c.a3 * (s.mag - c.c1)) * np.log(dist)
            + np.select(
                [(s.mag < c.c1), True], [c.a4 * (s.mag - c.c1), c.a5 * (s.mag - c.c1)]
            )
        )
        # Site model
        v_1 = np.select(
            [
                self.PERIODS < 0,
                self.PERIODS <= 0.50,
                self.PERIODS <= 1.00,
                self.PERIODS <= 2.00,
                self.PERIODS <= 10.0,
            ],
            # Use minimum on period to prevent runtime warning.
            [
                862,
                1500.0,
                np.exp(8 - 0.795 * np.log(np.maximum(0.001, self.PERIODS) / 0.21)),
                np.exp(6.76 - 0.297 * np.log(np.maximum(0.001, self.PERIODS))),
                700,
            ],
        )
        v_s30_lim = np.minimum(s.v_s30, v_1)
        vs_ratio = v_s30_lim / c.v_lin
        f5 = c.a10 * np.log(vs_ratio) - np.select(
            [v_s30_lim < c.v_lin, True],
            [
                -c.b * np.log(s.pga_ref + c.c)
                + c.b * np.log(s.pga_ref + c.c * vs_ratio ** c.n),
                c.b * c.n * np.log(vs_ratio),
            ],
        )

        ln_ratio = f1 + f5
        if s.mechanism == "RS":
            ln_ratio += c.a6
        elif s.mechanism == "NS":
            ln_ratio += c.a7

        return ln_ratio

    def _calc_ln_std_within(self):
        c = self.COEFF
        mag = self._scenario.mag
        return np.select(
            [mag < 5, mag <= 7, True],
            [c.s1, c.s1 + 0.5 * (c.s2 - c.s1) * (mag - 5), c.s2],
        )

    def _calc_ln_std_between(self):
        c = self.COEFF
        mag = self._scenario.mag
        return np.select(
            [mag < 5, mag <= 7, True],
            [c.s3, c.s3 + 0.5 * (c.s4 - c.s3) * (mag - 5), c.s4],
        )

    def _calc_ln_std(self):
        between = self._calc_ln_std_between()
        within = self._calc_ln_std_within()

        return np.sqrt(between ** 2 + within ** 2)
