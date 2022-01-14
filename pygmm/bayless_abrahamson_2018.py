# -*- coding: utf-8 -*-
"""Bayless and Abrahamson (2018, :cite:`bayless19`) model."""
from typing import Optional
from typing import Union

import numpy as np

from . import model

# Based on code from Artie Rodgers
__author__ = "Albert Kottke"


class BaylessAbrahamson18(model.Model):
    """Bayless and Abrahamson (2018, :cite:`bayless19`) model.

    This Fourier-amplitude spectra model was developed for active tectonic regions.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Bayless and Abrahamson (2018)"
    ABBREV = "BA18"

    # Reference velocity (m/s)
    V_REF = 1000.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("bayless_abrahamson_2018.csv", 2)

    FREQS = COEFF["freq_hz"]

    IDX_5Hz = 174
    IDX_MAX = 238

    PARAMS = [
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("mag", True, 3.0, 8.0),
        model.NumericParameter("v_s30", True, 180, 1500),
        model.NumericParameter("depth_tor", True, 0, 20),
        model.NumericParameter("depth_1_0", False),
        model.CategoricalParameter("mechanism", False, ["SS", "NS", "RS"], "SS"),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)

        # Taek the reference intensity at 5 Hz
        ln_eas_ref = self._calc_ln_eas(None)
        self._ln_eas = self._calc_ln_eas(ln_eas_ref)
        self._ln_std = self._calc_ln_std()

    @property
    def freqs(self):
        return self.FREQS

    @property
    def ln_eas(self):
        return self._ln_eas

    @property
    def eas(self):
        return np.exp(self._ln_eas)

    @property
    def ln_std(self):
        return self._ln_std

    def _check_inputs(self) -> None:
        super(BaylessAbrahamson18, self)._check_inputs()
        s = self._scenario

        if s["depth_1_0"] is None:
            s["depth_1_0"] = self.calc_depth_1_0(s["v_s30"])

    def _calc_ln_eas(self, ln_eas_ref: Optional[float]) -> Union[float, np.ndarray]:
        """Compute the effective amplitude."""
        s = self._scenario

        # Constants
        c4a = -0.5
        mbreak = 6.0

        if ln_eas_ref is None:
            v_s30 = self.V_REF
            # Only perform the 5 Hz calculation
            c = self.COEFF[self.IDX_MAX]
        else:
            v_s30 = s.v_s30
            # Only use the coefficients only the usable frequency range
            c = self.COEFF[: (self.IDX_MAX + 1)]

        c11 = np.select(
            [v_s30 <= 200, v_s30 <= 300, v_s30 <= 500], [c.c11a, c.c11b, c.c11c], c.c11d
        )

        ln_eas = (
            c.c1
            + c.c2 * (s.mag - mbreak)
            + ((c.c2 - c.c3) / c.cn) * np.log(1 + np.exp(c.cn * (c.cM - s.mag)))
            + c.c4
            * np.log(s.dist_rup + c.c5 * np.cosh(c.c6 * np.maximum(s.mag - c.chm, 0)))
            + (c4a - c.c4) * np.log(np.sqrt(s.dist_rup ** 2 + 50 ** 2))
            + c.c7 * s.dist_rup
            + c.c8 * np.log(np.minimum(v_s30, 1000) / self.V_REF)
            + c.c9 * s.depth_tor
            + c11
            * np.log(
                (np.minimum(s.depth_1_0, 2) + 0.01)
                / (self.calc_depth_1_0(v_s30) + 0.01)
            )
        )

        if s.mechanism == "NS":
            ln_eas += c.c10

        if ln_eas_ref is not None:
            ln_eas += self._calc_site_response(ln_eas_ref)

            # Extrapolate to 100 Hz using the kappa
            kappa = np.exp(-0.4 * np.log(v_s30 / 760) - 3.5)
            # Diminuition operator
            freq_max = self.FREQS[self.IDX_MAX]
            dimin = np.exp(
                -np.pi * kappa * (self.FREQS[(self.IDX_MAX + 1) :] - freq_max)
            )
            ln_eas = np.r_[ln_eas, ln_eas[self.IDX_MAX] + np.log(dimin)]

        return ln_eas

    def _calc_site_response(self, ln_eas_ref: float) -> np.ndarray:
        v_s30 = self.scenario.v_s30
        # Only use the coefficients only the usable frequency range
        c = self.COEFF[: (self.IDX_MAX + 1)]

        I_R = np.exp(1.238 + 0.846 * ln_eas_ref)

        f_sl = c.c8 * np.log(min(v_s30, 1000) / 1000)

        f2 = c.f4 * (
            np.exp(c.f5 * (min(v_s30, self.V_REF) - 360))
            - np.exp(c.f5 * (self.V_REF - 360))
        )

        f_nl = f2 + np.log((I_R + c.f3) / c.f3)
        f_s = f_sl + f_nl

        return f_s

    def _calc_ln_std(self) -> np.ndarray:
        """Compute the total standard deviation."""
        s = self._scenario
        c = self.COEFF

        def interp(left, right):
            return np.clip(left + (right - left) / 2 * (s.mag - 4), left, right)

        tau = interp(c.s1, c.s2)
        phi_s2s = interp(c.s3, c.s4)
        phi_ss = interp(c.s5, c.s6)

        ln_std = np.sqrt(tau ** 2 + phi_s2s ** 2 + phi_ss ** 2 + c.c1a ** 2)
        return ln_std

    @staticmethod
    def calc_depth_1_0(v_s30: float) -> float:
        """Compute the generic depth to 1 km/sec based on Vs30.

        Parameters
        ---------
        v_s30 : float
            site condition. Set `v_s30` to the reference
            velocity (e.g., 1180 m/s) for the reference response.

        Returns
        -------
        depth_1_0 : float
            depth to a shear-wave velocity of 1,000 m/sec
            (:math:`Z_{1.0}`, km).
        """
        power = 4
        v_ref = 610
        slope = -7.67 / power
        return (
            np.exp(
                slope
                * np.log(
                    (v_s30 ** power + v_ref ** power)
                    / (1360.0 ** power + v_ref ** power)
                )
            )
            / 1000
        )
