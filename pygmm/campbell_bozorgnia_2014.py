# -*- coding: utf-8 -*-
"""Model for the Campbell and Bozorgnia (2014) ground motion model."""
import logging
from typing import Optional

import numpy as np

from . import model
from .chiou_youngs_2014 import ChiouYoungs2014 as CY14
from .types import ArrayLike

__author__ = "Albert Kottke"


class CampbellBozorgnia2014(model.GroundMotionModel):
    """Campbell and Bozorgnia (2014, :cite:`campbell14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.
    """

    NAME = "Campbell & Bozorgnia (2014)"
    ABBREV = "CB14"

    # Reference velocity (m/sec)
    V_REF = 1100.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("campbell_bozorgnia_2014.csv", 2)

    PERIODS = COEFF["period"]

    # Period independent model coefficients
    COEFF_C = 1.88
    COEFF_N = 1.18
    COEEF_H_4 = 1

    INDICES_PSA = np.arange(21)
    INDEX_PGA = -2
    INDEX_PGV = -1

    PARAMS = [
        model.NumericParameter("depth_1_0", False),
        model.NumericParameter("depth_2_5", False, 0, 10),
        model.NumericParameter("depth_bor", False),
        model.NumericParameter("depth_bot", False, default=15.0),
        model.NumericParameter("depth_hyp", False, 0, 20),
        model.NumericParameter("depth_tor", False, 0, 20),
        model.NumericParameter("dip", True, 15, 90),
        model.NumericParameter("dist_jb", True),
        model.NumericParameter("dist_rup", True, None, 300),
        model.NumericParameter("dist_x", True),
        model.NumericParameter("mag", True, 3.3, 8.5),
        model.NumericParameter("v_s30", True, 150, 1500),
        model.NumericParameter("width", False),
        model.CategoricalParameter(
            "region",
            False,
            ["global", "california", "japan", "italy", "china"],
            "global",
        ),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def _check_inputs(self) -> None:
        """Check the inputs."""
        super()._check_inputs()
        s = self._scenario

        for mech, limit in [("SS", 8.5), ("RS", 8.0), ("NS", 7.5)]:
            if mech == s.mechanism and s.mag > limit:
                logging.warning(
                    "Magnitude of %g is greater than the recommended limit of"
                    "%g for %s style faults",
                    s.mag,
                    limit,
                    mech,
                )
        if s.depth_2_5 is None:
            s.depth_2_5 = self.calc_depth_2_5(s.v_s30, s.region, s.depth_1_0)

        if s.depth_tor is None:
            s.depth_tor = CY14.calc_depth_tor(s.mag, s.mechanism)

        if s.width is None:
            s.width = CampbellBozorgnia2014.calc_width(
                s.mag, s.dip, s.depth_tor, s.depth_bot
            )

        if s.depth_bor is None:
            s.depth_bor = self.calc_depth_bor(s.depth_tor, s.dip, s.width)

        if s.depth_hyp is None:
            s.depth_hyp = CampbellBozorgnia2014.calc_depth_hyp(
                s.mag, s.dip, s.depth_tor, s.depth_bor
            )

    def __init__(self, scenario: model.Scenario):
        """Initialize the model.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
        """
        super().__init__(scenario)

        pga_ref = np.exp(self._calc_ln_resp(np.nan, self.V_REF)[self.INDEX_PGA])
        self._ln_resp = self._calc_ln_resp(pga_ref, self._scenario.v_s30)
        self._ln_std = self._calc_ln_std(pga_ref)

    def _calc_ln_resp(self, pga_ref: float, v_s30: float) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        pga_ref : float
            peak ground acceleration (g) at the reference
            condition. If :class:`np.nan`, then no site term is applied.
        v_s30 : float
            time-averaged shear-wave velocity over the top 30 m
            of the site (:math:`V_{s30}`, m/s).

        Returns
        -------

            class:`np.array`: Natural log of the response.

        """
        c = self.COEFF
        s = self._scenario

        # Magnitude term
        f_mag = c.c_0 + c.c_1 * s.mag
        for min_mag, slope in ([4.5, c.c_2], [5.5, c.c_3], [6.5, c.c_4]):
            if min_mag < s.mag:
                f_mag += slope * (s.mag - min_mag)
            else:
                break

        # Geometric attenuation term
        f_dis = (c.c_5 + c.c_6 * s.mag) * np.log(np.sqrt(s.dist_rup ** 2 + c.c_7 ** 2))

        # Style of faulting term
        taper = np.clip(s.mag - 4.5, 0, 1)
        if s.mechanism == "RS":
            f_flt = c.c_8 * taper
        elif s.mechanism == "NS":
            f_flt = c.c_9 * taper
        else:
            f_flt = 0

        # Hanging-wall term
        R_1 = s.width * np.cos(np.radians(s.dip))
        R_2 = 62 * s.mag - 350
        if s.dist_x < 0:
            f_hngRx = 0
        elif s.dist_x <= R_1:
            ratio = s.dist_x / R_1
            f_hngRx = c.h_1 + c.h_2 * ratio + c.h_3 * ratio ** 2
        else:
            ratio = (s.dist_x - R_1) / (R_2 - R_1)
            f_hngRx = np.maximum(0, c.h_4 + c.h_5 * ratio + c.h_6 * ratio ** 2)

        if s.dist_rup == 0:
            f_hngRrup = 1
        else:
            f_hngRrup = (s.dist_rup - s.dist_jb) / s.dist_rup

        if s.mag <= 5.5:
            f_hngM = 0
        else:
            f_hngM = np.minimum(s.mag - 5.5, 1) * (1 + c.a_2 * (s.mag - 6.5))

        f_hngZ = 0 if s.depth_tor > 16.66 else 1 - 0.06 * s.depth_tor
        f_hngDip = (90 - s.dip) / 45

        f_hng = c.c_10 * f_hngRx * f_hngRrup * f_hngM * f_hngZ * f_hngDip

        site_term = self.calc_site_term(pga_ref, v_s30, s.depth_2_5, s.region)

        # Hypocentral depth term
        f_hypH = np.clip(s.depth_hyp - 7, 0, 13)
        f_hypM = c.c_17 + (c.c_18 - c.c_17) * np.clip(s.mag - 5.5, 0, 1)
        f_hyp = f_hypH * f_hypM

        # Fault dip term
        f_dip = c.c_19 * s.dip * np.clip(5.5 - s.mag, 0, 1)

        # Anaelastic attenuation term
        if s.region in ["japan", "italy"]:
            dc_20 = c.dc_20jp
        elif s.region == ["china"]:
            dc_20 = c.dc_20ch
        else:
            dc_20 = c.dc_20ca

        f_atn = (c.c_20 + dc_20) * max(s.dist_rup - 80, 0)

        ln_resp = f_mag + f_dis + f_flt + f_hng + site_term + f_hyp + f_dip + f_atn
        return ln_resp

    @classmethod
    def calc_site_term(
        cls,
        pga_ref: float,
        v_s30: float,
        depth_2_5: float,
        region: str = "california",
    ) -> ArrayLike:
        """Calculate the site term, which includes site and basin effects.

        Parameters
        ----------
        resp_ref :  array_like, optional
            response at the reference condition
        v_s30 : float
            site condition. Set `v_s30` to the reference
            velocity (e.g., 1180 m/s) for the reference response.
        depth_2_5 : float
            depth to the 2.5 km∕s shear-wave velocity horizon beneath the site,
            :math:`Z_{2.5}` in (km).
        region : str, optional
            region of basin model. Valid options:
            "global", "california", "japan", "italy", "china". If
            *None*, then "global" is used as the default value.

        Returns
        -------
        site_term: :class:`np.ndarray`
            site term that is applied to the natural log response.
        """

        c = cls.COEFF

        # Site term
        f_site = np.zeros_like(c.period)
        vs_ratio = v_s30 / c.k_1
        mask = v_s30 <= c.k_1
        f_site[mask] = (
            c.c_11 * np.log(vs_ratio)
            + c.k_2
            * (
                np.log(pga_ref + cls.COEFF_C * vs_ratio ** cls.COEFF_N)
                - np.log(pga_ref + cls.COEFF_C)
            )
        )[mask]
        f_site[~mask] = ((c.c_11 + c.k_2 * cls.COEFF_N) * np.log(vs_ratio))[~mask]

        if region == "japan":
            # Apply regional correction for Japan
            if v_s30 <= 200:
                f_site += (c.c_12 + c.k_2 * cls.COEFF_N) * (
                    np.log(vs_ratio) - np.log(200 / c.k_1)
                )
            else:
                f_site += (c.c_13 + c.k_2 * cls.COEFF_N) * np.log(vs_ratio)

        # Basin response term
        if np.isnan(pga_ref):
            # Use model to compute depth_2_5 for the reference velocity case
            depth_2_5 = cls.calc_depth_2_5(v_s30, region)
        else:
            depth_2_5 = depth_2_5

        if depth_2_5 <= 1:
            f_sed = c.c_14 * (depth_2_5 - 1)
            if region == "japan":
                f_sed += c.c_15 * (depth_2_5 - 1)
        elif depth_2_5 <= 3:
            f_sed = 0
        else:
            f_sed = (
                c.c_16 * c.k_3 * np.exp(-0.75) * (1 - np.exp(-0.25 * (depth_2_5 - 3)))
            )

        return f_site + f_sed

    def _calc_ln_std(self, pga_ref: ArrayLike) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Parameters
        ----------
        pga_ref : float
            peak ground acceleration (g) at the reference
            condition.

        Returns
        -------

            class:`np.array`: Logarithmic standard deviation.

        """
        c = self.COEFF
        s = self._scenario

        tau_lnY = c.tau_2 + (c.tau_1 - c.tau_2) * np.clip(5.5 - s.mag, 0, 1)
        phi_lnY = c.phi_2 + (c.phi_1 - c.phi_2) * np.clip(5.5 - s.mag, 0, 1)

        vs_ratio = s.v_s30 / c.k_1
        alpha = np.zeros_like(c.period)
        mask = s.v_s30 < c.k_1
        alpha[mask] = (
            c.k_2
            * pga_ref
            * (
                (pga_ref + self.COEFF_C * vs_ratio ** self.COEFF_N) ** (-1)
                - (pga_ref + self.COEFF_C) ** -1
            )
        )[mask]

        tau_lnPGA = tau_lnY[self.INDEX_PGA]
        tau = np.sqrt(
            tau_lnY ** 2
            + alpha ** 2 * tau_lnPGA ** 2
            + 2 * alpha * c.rho_lnPGAlnY * tau_lnY * tau_lnPGA
        )

        phi_lnPGA = phi_lnY[self.INDEX_PGA]
        phi_lnAF_PGA = self.COEFF["phi_lnAF"][self.INDEX_PGA]
        phi_lnPGA_B = np.sqrt(phi_lnPGA ** 2 - phi_lnAF_PGA ** 2)
        phi_lnY_B = np.sqrt(phi_lnY ** 2 - c.phi_lnAF ** 2)

        phi = np.sqrt(
            phi_lnY_B ** 2
            + c.phi_lnAF ** 2
            + alpha ** 2 * (phi_lnPGA ** 2 - phi_lnAF_PGA ** 2)
            + 2 * alpha * c.rho_lnPGAlnY * phi_lnY_B * phi_lnPGA_B
        )

        ln_std = np.sqrt(phi ** 2 + tau ** 2)

        return ln_std

    @staticmethod
    def calc_depth_2_5(
        v_s30: float, region: str = "global", depth_1_0: Optional[float] = None
    ) -> float:
        """Calculate the depth to a shear-wave velocity of 2.5 km/sec
        (:math:`Z_{2.5}`).

        Provide either `v_s30` or `depth_1_0`.

        Parameters
        ----------
        v_s30 : Optional[float]
            time-averaged shear-wave velocity over
            the top 30 m of the site (:math:`V_{s30}`, m/s).
            Keyword Args:
        region : Optional[str]
            region of the basin model. Valid values:
            "california", "japan". (Default value = 'global')
        depth_1_0 : Optional[float]
            depth to the 1.0 km∕s shear-wave
            velocity horizon beneath the site, :math:`Z_{1.0}` in (km).
            (Default value = None)

        Returns
        -------
        float
            estimated depth to a shear-wave velocity of 2.5 km/sec
        float
            estimated depth to a shear-wave velocity of 2.5 km/sec
            (km).

        """
        if v_s30:
            param = v_s30
            if region == "japan":
                # From equation 6.10 on page 63
                intercept = 5.359
                slope = 1.102
            else:
                # From equation 6.9 on page 63
                intercept = 7.089
                slope = 1.144

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 6.510
            # bar = 1.181
        elif depth_1_0:
            param = depth_1_0
            if region == "japan":
                # From equation 6.13 on page 64
                intercept = 0.408
                slope = 1.745
            else:
                # From equation 6.12 on page 64
                intercept = 1.392
                slope = 1.798

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 0.748
            # bar = 2.128
        else:
            raise NotImplementedError

        return np.exp(intercept - slope * np.log(param))

    @staticmethod
    def calc_depth_hyp(
        mag: float, dip: float, depth_tor: float, depth_bor: float
    ) -> float:
        """Estimate the depth to hypocenter.

        Parameters
        ----------
        mag : float
            moment magnitude of the event (:math:`M_w`)
        dip : float
            fault dip angle (:math:`\phi`, deg).
        depth_tor : float
            depth to the top of the rupture
            plane (:math:`Z_{tor}`, km).
        depth_bor : float
            depth to the bottom of the rupture
            plane (:math:`Z_{bor}`, km).

        Returns
        -------
        float
            estimated hypocenter depth (km)

        """
        # Equations 35, 36, and 37 of journal article
        ln_dZ = min(
            min(-4.317 + 0.984 * mag, 2.325) + min(0.0445 * (dip - 40), 0),
            np.log(0.9 * (depth_bor - depth_tor)),
        )

        depth_hyp = depth_tor + np.exp(ln_dZ)

        return depth_hyp

    @staticmethod
    def calc_width(
        mag: float, dip: float, depth_tor: float, depth_bot: float = 15.0
    ) -> float:
        """Estimate the fault width using Equation (39) of CB14.

        Parameters
        ----------
        mag : float
            moment magnitude of the event (:math:`M_w`)
        dip : float
            fault dip angle (:math:`\phi`, deg).
        depth_tor : float
            depth to the top of the rupture
            plane (:math:`Z_{tor}`, km).
            Keyword Args:
        depth_bot : Optional[float]
            depth to bottom of seismogenic crust
            (km). Used to calculate fault width if none is specified. If
            *None*, then a value of 15 km is used. (Default value = 15.0)

        Returns
        -------
        float
            estimated fault width (km)
        """
        return min(
            np.sqrt(10 ** ((mag - 4.07) / 0.98)),
            (depth_bot - depth_tor) / np.sin(np.radians(dip)),
        )

    @staticmethod
    def calc_depth_bor(depth_tor: float, dip: float, width: float) -> float:
        """Compute the depth to bottom of the rupture (km).

        Parameters
        ----------
        dip : float
            fault dip angle (:math:`\phi`, deg).
        depth_tor : float
            depth to the top of the rupture
            plane (:math:`Z_{tor}`, km).
        width : float
            Down-dip width of the fault.

        Returns
        -------
        float
            depth to bottom of the fault rupture (km)
        """
        return depth_tor + width * np.sin(np.radians(dip))
