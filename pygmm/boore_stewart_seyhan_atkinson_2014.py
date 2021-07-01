"""Boore, Stewart, Seyhan, and Atkinson (2014) ground motion model."""
import logging
from typing import Optional

import numpy as np

from . import model
from .chiou_youngs_2014 import ChiouYoungs2014 as CY14
from .types import ArrayLike

__author__ = "Albert Kottke"


class BooreStewartSeyhanAtkinson2014(model.GroundMotionModel):
    """Boore, Stewart, Seyhan, and Atkinson (2014, :cite:`boore14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    The BSSA14 model defines the following distance attenuation models:

        +--------------+-------------------------------+
        | Name         | Description                   |
        +--------------+-------------------------------+
        | global       | Global; California and Taiwan |
        +--------------+-------------------------------+
        | china_turkey | China and Turkey              |
        +--------------+-------------------------------+
        | italy_japan  | Italy and Japan               |
        +--------------+-------------------------------+

    and the following basin region models:

        +--------+---------------------+
        | Name   | Description         |
        +========+=====================+
        | global | Global / California |
        +--------+---------------------+
        | japan  | Japan               |
        +--------+---------------------+

    These are simplified into one regional parameter with the following
    possibilities:

        +-------------+--------------+------------+
        | Region      | Attenuation  | Basin      |
        +=============+==============+============+
        | global      | global       | global     |
        +-------------+--------------+------------+
        | california  | global       | global     |
        +-------------+--------------+------------+
        | china       | china_turkey | global     |
        +-------------+--------------+------------+
        | italy       | italy_japan  | global     |
        +-------------+--------------+------------+
        | japan       | italy_japan  | japan      |
        +-------------+--------------+------------+
        | new zealand | italy_japan  | global     |
        +-------------+--------------+------------+
        | taiwan      | global       | global     |
        +-------------+--------------+------------+
        | turkey      | china_turkey | global     |
        +-------------+--------------+------------+

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Boore, Stewart, Seyhan, and Atkinson (2014)"
    ABBREV = "BSSA14"

    # Reference shear-wave velocity in m/sec
    V_REF = 760.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("boore_stewart_seyhan_atkinson-2014.csv", 2)
    PERIODS = COEFF["period"]

    INDEX_PGV = 0
    INDEX_PGA = 1
    INDICES_PSA = np.arange(2, 107)

    LIMITS = dict(
        mag=(3.0, 8.5),
        dist_jb=(0.0, 300.0),
        v_s30=(150.0, 1500.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 3, 8.5),
        model.NumericParameter("depth_1_0", False),
        model.NumericParameter("dist_jb", True, None, 300.0),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.CategoricalParameter("mechanism", False, ["U", "SS", "NS", "RS"], "U"),
        model.CategoricalParameter(
            "region",
            False,
            [
                "global",
                "california",
                "china",
                "italy",
                "japan",
                "new_zealand",
                "taiwan",
                "turkey",
            ],
            "global",
        ),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
        """
        super(BooreStewartSeyhanAtkinson2014, self).__init__(scenario)
        pga_ref = np.exp(self._calc_ln_resp(np.nan)[self.INDEX_PGA])
        self._ln_resp = self._calc_ln_resp(pga_ref)
        self._ln_std = self._calc_ln_std()

    def _check_inputs(self) -> None:
        """Check the inputs."""
        super(BooreStewartSeyhanAtkinson2014, self)._check_inputs()
        s = self._scenario
        # Mechanism specific limits
        if s.mechanism == "SS":
            _min, _max = 3.0, 8.5
            if not (_min <= s.mag <= _max):
                logging.warning(
                    "Magnitude (%g) exceeds recommended bounds (%g to %g)"
                    " for a strike-slip earthquake!",
                    s.mag,
                    _min,
                    _max,
                )
        elif s.mechanism == "NS":
            _min, _max = 3.0, 7.0
            if not (_min <= s.mag <= _max):
                logging.warning(
                    "Magnitude (%g) exceeds recommended bounds (%g to %g)"
                    " for a normal-slip earthquake!",
                    s.mag,
                    _min,
                    _max,
                )

    def _calc_ln_resp(self, pga_ref: ArrayLike) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        pga_ref : float
            peak ground acceleration (g) at the reference
            condition. If :class:`np.nan`, then no site term is applied.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response

        """
        s = self._scenario
        c = self.COEFF

        # Compute the event term
        ########################
        if s.mechanism == "SS":
            event = np.array(c.e_1)
        elif s.mechanism == "NS":
            event = np.array(c.e_2)
        elif s.mechanism == "RS":
            event = np.array(c.e_3)
        else:
            # Unspecified
            event = np.array(c.e_0)

        mask = s.mag <= c.M_h
        event[mask] += (c.e_4 * (s.mag - c.M_h) + c.e_5 * (s.mag - c.M_h) ** 2)[mask]
        event[~mask] += (c.e_6 * (s.mag - c.M_h))[~mask]

        # Compute the distance terms
        ############################
        if s.region in ["china", "turkey"]:
            dc_3 = c.dc_3ct
        elif s.region in ["italy", "japan"]:
            dc_3 = c.dc_3ij
        else:
            # s.region in 'global', 'california', 'new_zealand', 'taiwan'
            dc_3 = c.dc_3global

        dist = np.sqrt(s.dist_jb ** 2 + c.h ** 2)
        path = (c.c_1 + c.c_2 * (s.mag - c.M_ref)) * np.log(dist / c.R_ref) + (
            c.c_3 + dc_3
        ) * (dist - c.R_ref)

        if np.isnan(pga_ref):
            # Reference condition. No site effect
            site = 0
        else:
            # Compute the site term
            site = self.calc_site_term(
                pga_ref, s.v_s30, s.depth_1_0, s.get("depth_1_0", None)
            )

        ln_resp = event + path + site
        return ln_resp

    @classmethod
    def calc_site_term(
        cls,
        pga_ref: float,
        v_s30: float,
        depth_1_0: Optional[float],
        region: str = "california",
    ) -> ArrayLike:
        """Calculate the site term, which includes site and basin effects.

        Parameters
        ----------
        pga_ref : float
            peak ground acceleration (g) at the reference
            condition. If :class:`np.nan`, then no site term is applied.
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

        f_lin = c.c * np.log(np.minimum(v_s30, c.V_c) / c.V_ref)

        # Add the nonlinearity to the site term
        f_2 = c.f_4 * (
            np.exp(c.f_5 * (min(v_s30, 760) - 360.0)) - np.exp(c.f_5 * (760.0 - 360.0))
        )
        f_nl = c.f_1 + f_2 * np.log((pga_ref + c.f_3) / c.f_3)

        # Add the basin effect to the site term
        F_dz1 = np.zeros_like(c.period)

        # Compute the average from the Chiou and Youngs (2014)
        # model convert from m to km.
        ln_mz1 = np.log(CY14.calc_depth_1_0(v_s30, region))

        if depth_1_0 is not None:
            delta_depth_1_0 = depth_1_0 - np.exp(ln_mz1)
        else:
            delta_depth_1_0 = 0.0

        mask = c.period >= 0.65
        F_dz1[mask] = np.minimum(c.f_6 * delta_depth_1_0, c.f_7)[mask]

        site = f_lin + f_nl + F_dz1

        return site

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        c = self.COEFF
        s = self._scenario

        # Uncertainty model
        tau = c.tau_1 + (c.tau_2 - c.tau_1) * (np.clip(s.mag, 4.5, 5.5) - 4.5)
        phi = c.phi_1 + (c.phi_2 - c.phi_1) * (np.clip(s.mag, 4.5, 5.5) - 4.5)

        # Modify phi for Vs30
        phi -= c.dphi_V * np.clip(np.log(c.V_2 / s.v_s30) / np.log(c.V_2 / c.V_1), 0, 1)

        # Modify phi for R
        phi += c.dphi_R * np.clip(
            # Maximum added for zero distance caes
            np.log(np.maximum(s.dist_jb, 0.1) / c.R_1) / np.log(c.R_2 / c.R_1),
            0,
            1,
        )

        ln_std = np.sqrt(phi ** 2 + tau ** 2)
        return ln_std
