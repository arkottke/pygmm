"""Abrahamson and Silva (1996, :cite:`abrahamson96`) duration model."""
from __future__ import division

import numpy as np

from . import model

__author__ = "Albert Kottke"


class AbrahamsonSilva1996(model.Model):
    """Abrahamson and Silva (1996, :cite:`abrahamson96`) duration model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Abrahamson Silva (1996)"
    ABBREV = "AS96"

    PARAMS = [
        model.NumericParameter("mag", True, 4, 7.5),
        model.NumericParameter("dist_rup", True, 0, 250),
        # FIXME add site_cond to Scenario
        model.CategoricalParameter("site_cond", True, ["soil", "rock"]),
    ]

    def __init__(self, scenario):
        super(AbrahamsonSilva1996, self).__init__(scenario)

        s = self._scenario

        stress_drop = np.exp(5.204 + 0.851 * (s.mag - 6))
        moment = 10 ** (1.5 * s.mag + 16.05)

        self._ln_dur = np.log(
            (stress_drop / moment) ** (-1 / 3) / (4.9e6 * 3.2)
            + 0.805 * (1 if s.site_cond == "soil" else 0)
            + 0.063 * max(s.dist_rup - 10, 0)
        ) + self.calc_ln_dur_incr(0.75)
        self._std_err = 0.55

    @property
    def duration(self):
        return np.exp(self._ln_dur)

    @property
    def std_err(self):
        return self._std_err

    @staticmethod
    def calc_ln_dur_incr(nias):
        # Mask out inappropriate values
        nias = np.asarray(nias)
        mask = (nias < 0.10) | (0.95 < nias)
        nias[mask] = np.nan

        # Compute the increment due to the difference in the normalized Arias
        # intensity
        ln_i_ratio = np.log((nias - 0.05) / (1 - nias))
        ln_dur_incr = -0.532 + 0.552 * ln_i_ratio - 0.0262 * ln_i_ratio ** 2

        return ln_dur_incr

    def interp(self, nias, stds=None):
        # Mask out inappropriate values
        nias = np.asarray(nias)
        mask = (nias < 0.10) | (0.95 < nias)
        nias[mask] = np.nan

        ln_dur = self._ln_dur + self.calc_ln_dur_incr(nias)

        if stds is not None:
            std_errs = np.interp(
                nias,
                [
                    0.10,
                    0.15,
                    0.20,
                    0.25,
                    0.30,
                    0.35,
                    0.40,
                    0.45,
                    0.50,
                    0.55,
                    0.60,
                    0.65,
                    0.70,
                    0.75,
                    0.80,
                    0.85,
                    0.90,
                    0.95,
                ],
                [
                    0.843,
                    0.759,
                    0.713,
                    0.691,
                    0.674,
                    0.660,
                    0.646,
                    0.636,
                    0.628,
                    0.616,
                    0.605,
                    0.594,
                    0.582,
                    0.565,
                    0.545,
                    0.528,
                    0.510,
                    0.493,
                ],
            )
            ln_dur = ln_dur + np.array(stds)[:, np.newaxis] * std_errs

        return np.exp(ln_dur)
