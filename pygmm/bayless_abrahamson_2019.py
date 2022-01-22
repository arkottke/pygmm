# -*- coding: utf-8 -*-
"""Bayless and Abrahamson (2019, :cite:`bayless19`) model."""
from typing import Optional
from typing import Union

import matplotlib.pyplot as plt
import numpy as np

from . import model

# Based on code from Artie Rodgers
__author__ = "Albert Kottke"


class BaylessAbrahamson2019(model.Model):
    """Bayless and Abrahamson (2019, :cite:`bayless19`) model.

    This Fourier-amplitude spectra model was developed for active tectonic regions.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Bayless and Abrahamson (2019)"
    ABBREV = "BA19"

    # Reference velocity (m/s)
    V_REF = 760.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("bayless_abrahamson_2019.csv", 2)

    FREQS = COEFF["freq_hz"]

    IDX_5Hz = 170
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
        super(BaylessAbrahamson2019, self)._check_inputs()
        s = self._scenario

        if s["depth_1_0"] is None:
            s["depth_1_0"] = self.calc_depth_1_0(s["v_s30"])

    def _calc_ln_eas(self, ln_eas_ref: Optional[float]) -> Union[float, np.ndarray]:
        """Compute the effective amplitude."""
        s = self._scenario

        if ln_eas_ref is None:
            # Compute the motion at the reference condition
            v_s30 = self.V_REF
            depth_1_0 = self.calc_depth_1_0(self.V_REF)
            # Only perform the 5 Hz calculation
            c = self.COEFF[self.IDX_5Hz]
        else:
            # Use site-specific properties
            v_s30 = s.v_s30
            depth_1_0 = s.depth_1_0
            # Only use the coefficients only the usable frequency range
            c = self.COEFF[: (self.IDX_MAX + 1)]

        f_M = (
            c.c1
            + c.c2 * (s.mag - 6)
            + ((c.c2 - c.c3) / c.cn) * np.log(1 + np.exp(c.cn * (c.cM - s.mag)))
        )

        f_P = (
            c.c4
            * np.log(s.dist_rup + c.c5 * np.cosh(c.c6 * np.maximum(s.mag - c.chm, 0)))
            + (-0.5 - c.c4) * np.log(np.sqrt(s.dist_rup ** 2 + 50 ** 2))
            + c.c7 * s.dist_rup
        )

        f_Ztor = c.c9 * min(s.depth_tor, 20)

        f_NM = c.c10 if s.mechanism == "NS" else 0

        c11 = np.select(
            [v_s30 <= 200, v_s30 <= 300, v_s30 <= 500], [c.c11a, c.c11b, c.c11c], c.c11d
        )
        f_Z1 = c11 * np.log(
            (min(depth_1_0, 2) + 0.01) / (self.calc_depth_1_0(v_s30) + 0.01)
        )

        ln_eas = f_M + f_P + f_Ztor + f_NM + f_Z1

        if ln_eas_ref is not None:
            # Add the site term
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

        f_nl = f2 * np.log((I_R + c.f3) / c.f3)

        # Use the minimum value at higher frequencies as described in the text and
        # referenced to Bayless and Abrahamson (2018b)
        i = np.argmin(f_nl)
        f_nl[(i + 1) :] = f_nl[i]
        width = 7
        kernel = np.ones(width) / width
        f_nl = np.convolve(f_nl, kernel, mode="same")
        # End needs to be forced down
        f_nl[-width:] = f_nl[i]

        f_s = f_sl + f_nl

        if v_s30 < 800:
            freqs = self.FREQS[: (self.IDX_MAX + 1)]
            fig, ax = plt.subplots()
            ax.plot(freqs, f_nl)

            ax.set(xscale="log")

            fig.savefig("test_fnl.png", dpi=200)

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
