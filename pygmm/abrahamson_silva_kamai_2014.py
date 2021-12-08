# -*- coding: utf-8 -*-
"""Abrahamson, Silva, and Kamai (2014, :cite:`abrahamson14`) model."""
import numpy as np
from scipy.interpolate import interp1d

from . import model
from .types import ArrayLike

__author__ = "Albert Kottke"


class AbrahamsonSilvaKamai2014(model.GroundMotionModel):
    """Abrahamson, Silva, and Kamai (2014, :cite:`abrahamson14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Abrahamson, Silva, & Kamai (2014)"
    ABBREV = "ASK14"

    # Reference velocity (m/sec)
    V_REF = 1180.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("abrahamson_silva_kamai_2014.csv", 2)

    PERIODS = COEFF["period"]

    INDICES_PSA = np.arange(22)
    INDEX_PGA = -2
    INDEX_PGV = -1

    PARAMS = [
        model.NumericParameter("dist_rup", True, None, 300),
        model.NumericParameter("dist_jb", True),
        model.NumericParameter("mag", True, 3, 8.5),
        model.NumericParameter("v_s30", True, 180, 1000),
        model.NumericParameter("depth_1_0", False),
        model.NumericParameter("depth_tor", False),
        model.NumericParameter("dip", True),
        model.NumericParameter("dist_crjb", False, default=15),
        model.NumericParameter("dist_x", False),
        model.NumericParameter("dist_y0", False),
        model.NumericParameter("width", False),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
        model.CategoricalParameter(
            "region",
            False,
            ["global", "california", "china", "italy", "japan", "taiwan"],
            "global",
        ),
        model.CategoricalParameter(
            "vs_source", False, ["measured", "inferred"], "measured"
        ),
        model.CategoricalParameter("is_aftershock", False, [True, False], False),
        model.CategoricalParameter("on_hanging_wall", False, [True, False], False),
    ]

    def _check_inputs(self) -> None:
        """Check the inputs."""
        super(AbrahamsonSilvaKamai2014, self)._check_inputs()
        s = self._scenario
        if s["width"] is None:
            s["width"] = self.calc_width(s.mag, s.dip)

        if s["depth_tor"] is None:
            s["depth_tor"] = self.calc_depth_tor(s.mag)

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super(AbrahamsonSilvaKamai2014, self).__init__(scenario)
        # Compute the response at the reference velocity
        resp_ref = np.exp(self._calc_ln_resp(self.V_REF, np.nan))

        self._ln_resp = self._calc_ln_resp(self._scenario.v_s30, resp_ref)
        self._ln_std = self._calc_ln_std(resp_ref)

    def _calc_ln_resp(self, v_s30: float, resp_ref: ArrayLike) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        v_s30 : float
            site condition. Set `v_s30` to the reference
            velocity (e.g., 1180 m/s) for the reference response.
        resp_ref :  array_like, optional
            response at the reference condition. Required if `v_s30` is not
            equal to reference velocity.

        Returns
        -------
        ln_resp: :class:`np.ndarray`
            natural log of the response.

        """
        c = self.COEFF
        s = self._scenario

        # Magnitude scaling
        f1 = self._calc_f1()

        if s.on_hanging_wall:
            # Hanging-wall term
            f4 = self._calc_f4()
        else:
            f4 = 0

        # Depth to top of rupture term
        f6 = c.a15 * np.clip(s.depth_tor / 20, 0, 1)
        # Style of faulting
        if s.mechanism == "RS":
            f7 = c.a11 * np.clip(s.mag - 4, 0, 1)
            f8 = 0
        elif s.mechanism == "NS":
            f7 = 0
            f8 = c.a12 * np.clip(s.mag - 4, 0, 1)
        else:
            f7, f8 = 0, 0

        site_term = self.calc_site_term(resp_ref, v_s30, s.depth_1_0, s.region)

        # Aftershock term
        if s.is_aftershock:
            f11 = c.a14 * np.clip(1 - (s.dist_crjb - 5) / 10, 0, 1)
        else:
            f11 = 0

        if s.region == "taiwan":
            v_1 = np.exp(-0.35 * np.log(np.clip(c.period, 0.5, 3) / 0.5) + np.log(1500))
            vs_ratio = np.minimum(v_s30, v_1) / c.v_lin
            freg = c.a31 * np.log(vs_ratio) + c.a25 * s.dist_rup
        elif s.region == "china":
            freg = c.a28 * s.dist_rup
        elif s.region == "japan":
            f13 = interp1d(
                [150, 250, 350, 450, 600, 850, 1150],
                np.c_[c.a36, c.a37, c.a38, c.a39, c.a40, c.a41, c.a42],
                copy=False,
                bounds_error=False,
                fill_value=(c.a36, c.a42),
            )(v_s30)
            freg = f13 + c.a29 * s.dist_rup
        else:
            freg = 0

        return f1 + f4 + f6 + f7 + f8 + f11 + freg + site_term

    @classmethod
    def calc_site_term(
        cls,
        resp_ref: ArrayLike,
        v_s30: float,
        depth_1_0: float,
        region: str = "california",
    ):
        """Calculate the site term, which includes site and basin effects.

        Parameters
        ----------
        resp_ref :  array_like, optional
            response at the reference condition
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

        v_1 = np.exp(-0.35 * np.log(np.clip(c.period, 0.5, 3) / 0.5) + np.log(1500))

        vs_ratio = np.minimum(v_s30, v_1) / c.v_lin
        # Linear site model
        f5 = (c.a10 + c.b * c.n) * np.log(vs_ratio)
        # Nonlinear model
        mask = vs_ratio < 1
        f5[mask] = (
            c.a10 * np.log(vs_ratio)
            - c.b * np.log(resp_ref + c.c)
            + c.b * np.log(resp_ref + c.c * vs_ratio ** c.n)
        )[mask]

        # Basin term
        if v_s30 == cls.V_REF or depth_1_0 is None:
            # No basin response
            f10 = 0
        else:
            # Ratio between site depth_1_0 and model center
            ln_depth_ratio = np.log(
                (depth_1_0 + 0.01) / (cls.calc_depth_1_0(v_s30, region) + 0.01)
            )
            slope = interp1d(
                [150, 250, 400, 700],
                np.c_[c.a43, c.a44, c.a45, c.a46],
                copy=False,
                bounds_error=False,
                fill_value=(c.a43, c.a46),
            )(v_s30)
            f10 = slope * ln_depth_ratio

        return f5 + f10

    def _calc_ln_std(self, psa_ref: ArrayLike) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Parameters
        ----------
        psa_ref : array_like
           spectral accelerations at the reference condition

        Returns
        -------
        ln_std :  :class:`np.ndarray`
            logarithmic standard deviation.
        """
        s = self._scenario
        c = self.COEFF

        if s.region == "japan":
            phi_al = c.s5 + (c.s6 - c.s5) * np.clip((s.dist_rup - 30) / 50, 0, 1)
        else:
            transition = np.clip((s.mag - 4) / 2, 0, 1)
            if s.vs_source == "measured":
                phi_al = c.s1m + (c.s2m - c.s1m) * transition
            else:
                phi_al = c.s1e + (c.s2e - c.s1e) * transition

        tau_al = c.s3 + (c.s4 - c.s3) * np.clip((s.mag - 5) / 2, 0, 1)
        tau_b = tau_al

        # Remove period independent site amplification uncertainty of 0.4
        phi_amp = 0.4
        phi_b = np.sqrt(phi_al ** 2 - phi_amp ** 2)

        # The partial derivative of the amplification with respect to
        # the reference intensity
        deriv = (-c.b * psa_ref) / (psa_ref + c.c) + (c.b * psa_ref) / (
            psa_ref + c.c * (s.v_s30 / c.v_lin) ** c.n
        )
        deriv[s.v_s30 >= c.v_lin] = 0
        tau = tau_b * (1 + deriv)
        phi = np.sqrt(phi_b ** 2 * (1 + deriv) ** 2 + phi_amp ** 2)

        ln_std = np.sqrt(phi ** 2 + tau ** 2)
        return ln_std

    @staticmethod
    def calc_width(mag: float, dip: float) -> float:
        """Compute the fault width based on equation in NGW2 spreadsheet.

        This equation is not provided in the paper.

        Parameters
        ----------
        mag : float
            moment magnitude of the event (:math:`M_w`)
        dip : float
            Fault dip angle (:math:`\phi`, deg)

        Returns
        -------
        width : float
            estimated fault width (:math:`W`, km)
        """
        return min(18 / np.sin(np.radians(dip)), 10 ** (-1.75 + 0.45 * mag))

    @staticmethod
    def calc_depth_tor(mag: float) -> float:
        """Calculate the depth to top of rupture (km).

        Parameters
        ----------
        mag : float
            moment magnitude of the event (:math:`M_w`)

        Returns
        -------
        depth_tor : float
            estimated depth to top of rupture (km)
        """
        return np.interp(mag, [5.0, 7.2], [7.8, 0])

    @staticmethod
    def calc_depth_1_0(v_s30: float, region: str = "california") -> float:
        """Estimate the depth to 1 km/sec horizon (:math:`Z_{1.0}`) based on
        :math:`V_{s30}` and region.

        This is based on equations 18 and 19 in the :cite:`abrahamson14`
        and differs from the equations in the :cite:`chiou14`.

        Parameters
        ----------
        v_s30 : float
            time-averaged shear-wave velocity over the top 30 m
            of the site (:math:`V_{s30}`, m/s).
            Keyword Args:
        region : str, optional
            region of basin model. Valid options: 'california', 'japan'. If
            *None*, then 'california' is used as the default value.

        Returns
        -------
        depth_1_0 : float
            depth to a shear-wave velocity of 1,000 m/sec
            (:math:`Z_{1.0}`, km).

        """
        if region in ["japan"]:
            # Japan
            power = 2
            v_ref = 412
            slope = -5.23 / power
        else:
            # Global
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

    def _calc_f1(self) -> np.ndarray:
        """Calculate the magnitude scaling parameter f1."""
        c = self.COEFF
        s = self._scenario

        # Magnitude dependent taper
        dist = np.sqrt(
            s.dist_rup ** 2 + (c.c4 - (c.c4 - 1) * np.clip(5 - s.mag, 0, 1)) ** 2
        )

        # Magnitude scaling
        # Need to copy c.a1 to that it isn't modified during the following
        # operations.
        f1 = np.array(c.a1)
        ma1 = s.mag <= c.m2
        f1[ma1] += (
            c.a4 * (c.m2 - c.m1)
            + c.a8 * (8.5 - c.m2) ** 2
            + c.a6 * (s.mag - c.m2)
            + c.a7 * (s.mag - c.m2) ** 2
            + (c.a2 + c.a3 * (c.m2 - c.m1)) * np.log(dist)
            + c.a17 * s.dist_rup
        )[ma1]

        f1[~ma1] += (
            c.a8 * (8.5 - s.mag) ** 2
            + (c.a2 + c.a3 * (s.mag - c.m1)) * np.log(dist)
            + c.a17 * s.dist_rup
        )[~ma1]

        ma2 = np.logical_and(~ma1, s.mag <= c.m1)
        f1[ma2] += (c.a4 * (s.mag - c.m1))[ma2]

        ma3 = np.logical_and(~ma1, s.mag > c.m1)
        f1[ma3] += (c.a5 * (s.mag - c.m1))[ma3]

        return f1

    def _calc_f4(self) -> np.ndarray:
        """Calculate the hanging-wall parameter f4."""
        # Move this into a decorator?
        c = self.COEFF
        s = self._scenario

        t1 = min(90 - s.dip, 60) / 45
        # Constant from page 1041
        a2hw = 0.2
        if s.mag <= 5.5:
            t2 = 0
        elif s.mag < 6.5:
            t2 = 1 + a2hw * (s.mag - 6.5) - (1 - a2hw) * (s.mag - 6.5) ** 2
        else:
            t2 = 1 + a2hw * (s.mag - 6.5)

        # Constants defined on page 1040
        r1 = s.width * np.cos(np.radians(s.dip))
        r2 = 3 * r1
        h1 = 0.25
        h2 = 1.5
        h3 = -0.75
        if s.dist_x < r1:
            t3 = h1 + h2 * (s.dist_x / r1) + h3 * (s.dist_x / r1) ** 2
        elif s.dist_x < r2:
            t3 = 1 - ((s.dist_x - r1) / (r2 - r1))
        else:
            t3 = 0

        t4 = np.clip(1 - s.depth_tor ** 2 / 100, 0, 1)

        if s.dist_y0 is None:
            t5 = np.clip(1 - s.dist_jb / 30, 0, 1)
        else:
            dist_y1 = s.dist_x * np.tan(np.radians(20))
            t5 = np.clip(1 - (s.dist_y0 - dist_y1) / 5, 0, 1)

        f4 = c.a13 * t1 * t2 * t3 * t4 * t5

        return f4
