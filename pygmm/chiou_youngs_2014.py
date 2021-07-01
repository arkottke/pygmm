# -*- coding: utf-8 -*-
"""Chiou and Youngs (2014, :cite:`chiou14`) model."""
import logging
from typing import Optional

import numpy as np

from . import model
from .types import ArrayLike

__author__ = "Albert Kottke"


class ChiouYoungs2014(model.GroundMotionModel):
    """Chiou and Youngs (2014, :cite:`chiou14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Chiou and Youngs (2014)"
    ABBREV = "CY14"

    # Reference velocity (m/s)
    V_REF = 1130.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("chiou_youngs_2014.csv", 2)

    PERIODS = COEFF["period"]

    INDICES_PSA = np.arange(24)
    INDEX_PGA = 24
    INDEX_PGV = 25

    PARAMS = [
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("dist_x", True),
        model.NumericParameter("dist_jb", True),
        model.NumericParameter("mag", True, 3.5, 8.5),
        model.NumericParameter("v_s30", True, 180, 1500),
        model.NumericParameter("depth_tor", False, 0, 20),
        model.NumericParameter("depth_1_0", False),
        model.NumericParameter("dpp_centered", False, default=0),
        model.NumericParameter("dip", True),
        model.CategoricalParameter("mechanism", False, ["U", "SS", "NS", "RS"], "U"),
        model.CategoricalParameter("on_hanging_wall", False, [True, False], False),
        model.CategoricalParameter(
            "region", False, ["california", "china", "italy", "japan"], "california"
        ),
        model.CategoricalParameter(
            "vs_source", False, ["measured", "inferred"], "measured"
        ),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)
        ln_resp_ref = self._calc_ln_resp_ref()
        self._ln_resp = self._calc_ln_resp_site(ln_resp_ref)
        self._ln_std = self._calc_ln_std(np.exp(ln_resp_ref))

    def _calc_ln_resp_ref(self) -> np.ndarray:
        """Calculate the natural logarithm of the reference response.

        Returns
        -------
        ln_resp_ref : class:`np.array`:
            natural log of the response

        """
        c = self.COEFF
        s = self._scenario

        cosh_mag = np.cosh(2 * max(s.mag - 4.5, 0))
        ln_resp = np.array(c.c_1)

        if s.mechanism == "RS":
            # Reverse fault term
            ln_resp += c.c_1a + c.c_1c / cosh_mag
        elif s.mechanism == "NS":
            # Normal fault term
            ln_resp += c.c_1b + c.c_1d / cosh_mag

        # Magnitude scaling
        ln_resp += c.c_2 * (s.mag - 6)
        ln_resp += (c.c_2 - c.c_3) / c.c_n * np.log(1 + np.exp(c.c_n * (c.c_m - s.mag)))

        # Top of rupture term relative to model average
        diff_depth_tor = s.depth_tor - self.calc_depth_tor(s.mag, s.mechanism)
        ln_resp += (c.c_7 + c.c_7b / cosh_mag) * diff_depth_tor

        # Dip angle term
        ln_resp += (c.c_11 + c.c_11b / cosh_mag) * np.cos(np.radians(s.dip)) ** 2

        # Distance terms
        ln_resp += c.c_4 * np.log(
            s.dist_rup + c.c_5 * np.cosh(c.c_6 * np.maximum(s.mag - c.c_hm, 0))
        )
        ln_resp += (c.c_4a - c.c_4) * np.log(np.sqrt(s.dist_rup ** 2 + c.c_rb ** 2))

        # Regional adjustment
        if s.region in ["japan", "italy"] and (6 < s.mag < 6.9):
            scale = c.gamma_ji
        elif s.region in ["china"]:
            scale = c.gamma_c
        else:
            scale = 1.0
        ln_resp += (
            scale
            * (c.c_gamma1 + c.c_gamma2 / np.cosh(np.maximum(s.mag - c.c_gamma3, 0)))
            * s.dist_rup
        )

        # Directivity term
        ln_resp += (
            c.c_8
            * max(1 - max(s.dist_rup - 40, 0) / 30, 0)
            * min(max(s.mag - 5.5, 0) / 0.8, 1)
            * np.exp(-c.c_8a * (s.mag - c.c_8b) ** 2)
            * s.dpp_centered
        )

        # Hanging wall term
        if s.on_hanging_wall:
            ln_resp += (
                c.c_9
                * np.cos(np.radians(s.dip))
                * (c.c_9a + (1 - c.c_9a) * np.tanh(s.dist_x / c.c_9b))
                * (1 - np.sqrt(s.dist_jb ** 2 + s.depth_tor ** 2) / (s.dist_rup + 1))
            )

        return ln_resp

    def _calc_ln_resp_site(self, ln_resp_ref: np.ndarray) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        ln_resp_ref : :class:`np.array`
            natural logarithm of the response at the reference site condition
            at each of the periods specified by the model coefficients.

        Returns
        -------
        ln_resp_rsite : class:`np.ndarray`:
            natural log of the response including the site effects

        """
        s = self._scenario

        site_term = self.calc_site_term(
            np.exp(ln_resp_ref), s.v_s30, s.depth_1_0, s.region
        )
        return ln_resp_ref + site_term

    @classmethod
    def calc_site_term(
        cls, resp_ref: ArrayLike, v_s30: float, depth_1_0: float, region: Optional[str]
    ) -> ArrayLike:
        """Calculate the site term, which includes site and basin effects.

        Parameters
        ----------
        resp_ref : :class:`np.array`
            natural logarithm of the response at the reference site condition
            at each of the periods specified by the model coefficients.
        v_s30 : float
            site condition. Set `v_s30` to the reference
            velocity (e.g., 1180 m/s) for the reference response.
        depth_1_0 : float
            depth to the 1.0 km∕s shear-wave velocity horizon beneath the site,
            :math:`Z_{1.0}` in (km).
        region : str, optional
            region of basin model. Valid options: 'california', 'japan'. If
            *None*, then 'california' is used as the default value.

        Returns
        -------
        site_term: :class:`np.ndarray`
            site term that is applied to the natural log response.
        """

        c = cls.COEFF

        if region in ["japan"]:
            phi_1 = c.phi_1jp
            phi_5 = c.phi_5jp
            phi_6 = c.phi_6jp
        else:
            phi_1 = c.phi_1
            phi_5 = c.phi_5
            phi_6 = c.phi_6

        ln_resp_ref = np.log(resp_ref)

        # Linear scaling
        site_term = phi_1 * min(np.log(v_s30 / 1130.0), 0)
        # Nonlinear scaling
        site_term += (
            c.phi_2
            * (
                np.exp(c.phi_3 * (min(v_s30, 1130.0) - 360.0))
                - np.exp(c.phi_3 * (1130.0 - 360.0))
            )
            * np.log((np.exp(ln_resp_ref) + c.phi_4) / c.phi_4)
        )
        # Basin model
        diff_depth_1_0 = 1000 * (depth_1_0 - cls.calc_depth_1_0(v_s30, region))
        site_term += phi_5 * (1 - np.exp(-diff_depth_1_0 / phi_6))

        return site_term

    def _calc_ln_std(self, resp_ref: np.ndarray) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Parameters
        ----------
        ln_resp_ref : :class:`np.array`
            natural logarithm of the response at the reference site condition
            at each of the periods specified by the model coefficients.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        c = self.COEFF
        s = self._scenario

        if s.region in ["japan"]:
            sigma_2 = c.sigma_2jp
        else:
            sigma_2 = c.sigma_2

        clipped_mag = np.clip(s.mag, 5.0, 6.5) - 5.0
        tau = c.tau_1 + (c.tau_2 - c.tau_1) / 1.5 * clipped_mag

        nl_0 = (
            c.phi_2
            * (
                np.exp(c.phi_3 * (min(s.v_s30, 1130.0) - 360.0))
                - np.exp(c.phi_3 * (1130.0 - 360.0))
            )
            * (resp_ref / (resp_ref + c.phi_4))
        )

        flag_meas = 1 if s.vs_source == "measured" else 0
        phi_nl = (c.sigma_1 + (sigma_2 - c.sigma_1) / 1.5 * clipped_mag) * np.sqrt(
            c.sigma_3 * (1 - flag_meas) + 0.7 * flag_meas + (1 + nl_0) ** 2
        )

        ln_std = np.sqrt((1 + nl_0) ** 2 * tau ** 2 + phi_nl ** 2)
        return ln_std

    def _check_inputs(self) -> None:
        """Check the inputs."""
        super(ChiouYoungs2014, self)._check_inputs()
        s = self._scenario

        if s.mechanism in ["RS", "NS"]:
            _min, _max = 3.5, 8.0
        else:
            _min, _max = 3.5, 8.5

        if not (_min <= s.mag <= _max):
            logging.warning(
                "Magnitude (%g) exceeds recommended bounds (%g to %g)"
                " for a %s earthquake!",
                s.mag,
                _min,
                _max,
                s.mechanism,
            )

        if s.get("depth_tor", None) is None:
            s["depth_tor"] = self.calc_depth_tor(s.mag, s.mechanism)

        if s.get("depth_1_0", None) is None:
            # Calculate depth (m) and convert to (km)
            s["depth_1_0"] = self.calc_depth_1_0(s.v_s30, s.region)

    @staticmethod
    def calc_depth_1_0(v_s30: float, region: str) -> float:
        """Calculate the depth to 1 km/sec (:math:`Z_{1.0}`).

        Parameters
        ----------
        v_s30 : float
            time-averaged shear-wave velocity over the top 30 m
            of the site (:math:`V_{s30}`, m/s).
        region : str
            basin region. Valid options: "california", "japan"

        Returns
        -------
        depth_1_0 : float
            estimated depth to a shear-wave velocity of 1 km/sec (km)

        """
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
                    (v_s30 ** power + v_ref ** power)
                    / (1360.0 ** power + v_ref ** power)
                )
            )
            / 1000
        )

    @staticmethod
    def calc_depth_tor(mag: float, mechanism: str) -> float:
        """Calculate an estimate of the depth to top of rupture (km).

        Parameters
        ----------
        mag : float
            moment magnitude of the event (:math:`M_w`)
        mechanism : str
            fault mechanism. Valid options: "U", "SS", "NS",
            "RS".

        Returns
        -------
        depth_tor : float
            estimated depth to top of rupture (km)

        """
        if mechanism == "RS":
            # Reverse and reverse-oblique faulting
            depth_tor = 2.704 - 1.226 * max(mag - 5.849, 0)
        else:
            # Combined strike-slip and normal faulting
            depth_tor = 2.673 - 1.136 * max(mag - 4.970, 0)

        return max(depth_tor, 0) ** 2
