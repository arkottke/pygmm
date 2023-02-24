"""Afshari and Stewart (2016, :cite:`afshari2016`) duration model."""
from __future__ import division

import numpy as np

from . import model


class AfshariStewart2016(model.Model):
    """Afshari and Stewart (2016, :cite:`afshari2016`) duration model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Afshari and Stewart (2016)"
    ABBREV = "AS16"

    PARAMS = [
        model.NumericParameter("mag", True, 3, 7.9),
        model.NumericParameter("dist_rup", True, 0, 200),
        model.NumericParameter("v_s30", True, 200, 1000),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
        model.NumericParameter("depth_1_0", False),
    ]

    def __init__(self, scenario):
        super(AfshariStewart2016, self).__init__(scenario)

        # scenario
        s = self._scenario

        # duration D_{5-75}, D_{5-95}, D_{20-80},

        # source coefficients
        mag_1 = np.array([5.35, 5.2, 5.2])
        mag_2 = np.array([7.15, 7.4, 7.4])
        mag_star = 6
        if s.mechanism is None:
            b_0 = np.array([1.280, 2.182, 0.8822])
            b_1 = np.array([5.576, 3.628, 6.182])
        elif s.mechanism == "NS":
            b_0 = np.array([1.555, 2.541, 1.409])
            b_1 = np.array([4.992, 3.170, 4.778])
        elif s.mechanism == "RS":
            b_0 = np.array([0.7806, 1.612, 0.7729])
            b_1 = np.array([7.061, 4.536, 6.188])
        elif s.mechanism == "SS":
            b_0 = np.array([1.279, 2.302, 0.8804])
            b_1 = np.array([5.578, 3.467, 6.188])
        b_2 = np.array([0.9011, 0.9443, 0.7414])
        b_3 = np.array([-1.684, -3.911, -3.164])
        # stress drop indices for duration measures

        stress_drop = np.exp(
            b_1
            + b_2 * (np.minimum(s.mag, mag_2) - mag_star)
            + b_3 * np.maximum((s.mag - mag_2), 0)
        )
        # corner frequency
        f_0 = self._corner_freq(stress_drop)

        # source duration
        F_E = (s.mag <= mag_1) * b_0 + (s.mag > mag_1) * (1 / f_0)

        # path coeficients
        c_1 = np.array([0.1159, 0.3165, 0.0646])
        c_2 = np.array([0.1065, 0.2539, 0.0865])
        c_3 = np.array([0.0682, 0.0932, 0.0373])
        c_4 = np.array([-0.2246, -0.3183, -0.4237])
        c_5 = np.array([0.0006, 0.0006, 0.0005])
        R_1 = 10
        R_2 = 50
        # path duration
        F_P = (
            c_1 * np.minimum(s.dist_rup, R_1)
            + c_2 * np.maximum((np.minimum(s.dist_rup, R_2) - R_1), 0)
            + c_3 * np.maximum((s.dist_rup - R_2), 0)
        )

        # site coefficients
        V_1 = 600
        V_ref = np.array([368.2, 369.9, 369.6])
        dz_1ref = 200
        # depth to bedrock duration
        if s.depth_1_0 is not None:
            dz_1 = s.depth_1_0 - self.calc_depth_1_0(s.v_s30, s.mechanism)
            F_dz1 = c_5 * (dz_1 if dz_1 <= dz_1ref else dz_1ref)
        else:
            F_dz1 = 0
        # site duration
        F_S = c_4 * np.log(np.minimum(s.v_s30, V_1) / V_ref) + F_dz1

        # aleatory coefficients
        tau_1 = np.array([0.28, 0.25, 0.30])
        tau_2 = np.array([0.25, 0.19, 0.19])
        phi_1 = np.array([0.54, 0.43, 0.56])
        phi_2 = np.array([0.41, 0.35, 0.45])
        # aleatory variability
        tau = tau_1 + (tau_2 - tau_1) * (
            np.minimum(np.maximum(s.mag, 6.5), 7.0) - 6.5
        ) / (7.0 - 6.5)
        phi = phi_1 + (phi_2 - phi_1) * (
            np.minimum(np.maximum(s.mag, 5.5), 5.75) - 5.5
        ) / (5.75 - 5.50)

        # total druation
        self._ln_dur = np.log(F_E + F_P) + F_S

        # aleatory standard deviation
        self._std_err = np.sqrt(tau**2 + phi**2)

    @staticmethod
    def calc_depth_1_0(v_s30: float, region: str = "california") -> float:
        if region in ["japan"]:
            # Japan
            power = 2
            v_ref = 412.39
            slope = -5.23 / power
        else:
            # Global
            power = 4
            v_ref = 570.94
            slope = -7.15 / power

        return (
            np.exp(
                slope
                * np.log(
                    (v_s30**power + v_ref**power)
                    / (1360.0**power + v_ref**power)
                )
            )
            / 1000
        )

    @property
    def duration(self):
        """Return the durations as a `np.recarray`. Values can be accessed with
        'D_5t75', 'D_5t95', or 'D_20t80'."""
        return AfshariStewart2016._as_recarray(np.exp(self._ln_dur))

    @property
    def std_err(self):
        """Return the standard errors of the durations as a `np.recarray`. Values can be
        accessed with 'D_5t75', 'D_5t95', or 'D_20t80'."""
        return AfshariStewart2016._as_recarray(self._std_err[:3])

    @staticmethod
    def _as_recarray(values):
        return np.rec.fromarrays(
            values,
            names=[
                "D_5t75",
                "D_5t95",
                "D_20t80",
            ],
        )

    def _corner_freq(self, stress_drop):
        """Corner frequency."""
        moment = 10 ** (1.5 * self.scenario.mag + 16.05)
        return (4.9e6 * 3.2) * (stress_drop / moment) ** (1 / 3)
