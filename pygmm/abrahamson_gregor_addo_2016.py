"""Abrahamson, Gregor, and Addo (2016) ground motion moodel."""
from __future__ import division

import numpy as np

from . import model

__author__ = "Albert Kottke"


class AbrahamsonGregorAddo2016(model.GroundMotionModel):
    """Abrahamson, Gregor, and Addo (2016) ground motion moodel.

    This model was developed for subduction regions and is commonly referred to
    as the BCHydro model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Abrahamson, Gregor, & Addo (2016)"
    ABBREV = "AGA16"

    # Reference shear-wave velocity in m/sec
    V_REF = 1000.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("abrahamson_gregor_addo_2016.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 23)

    # FIXME
    LIMITS = dict(
        mag=(3.0, 8.5),
        dist_jb=(0.0, 300.0),
        v_s30=(150.0, 1500.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 3, 8.5),
        model.NumericParameter("dist_rup", False),
        model.NumericParameter("dist_hyp", False),
        model.NumericParameter("depth_hyp", False),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
        model.CategoricalParameter(
            "tectonic_region", False, ["forearc", "backarc", "unknown"], "unknown"
        ),
    ]

    def __init__(self, scenario, adjust_c1=None, adjust_c4=0, scale_atten=1.0):
        """Initialize the model.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
        """
        super(AbrahamsonGregorAddo2016, self).__init__(scenario)

        n = len(self.COEFF)
        if adjust_c1 is None:
            # Default adjustments
            if self.scenario.event_type == "intraslab":
                self._adjust_c1 = -0.3 * np.ones(n)
            else:
                self._adjust_c1 = np.interp(
                    np.log(np.maximum(0.01, self.COEFF.period)),
                    np.log([0.3, 0.5, 1, 2, 3]),
                    [0.2, 0.1, 0, -0.1, -0.2],
                    left=0.2,
                    right=-0.2,
                )
        else:
            if isinstance(adjust_c1, float):
                self._adjust_c1 = adjust_c1 * np.ones(n)
            else:
                self._adjust_c1 = np.asarray(adjust_c1)

        self._adjust_c4 = adjust_c4
        self._scale_atten = scale_atten

        pga_ref = np.exp(self._calc_ln_resp(np.nan)[self.INDEX_PGA])
        self._ln_resp = self._calc_ln_resp(pga_ref)
        self._ln_std = self._calc_ln_std()

    @property
    def adjust_c1(self):
        return self._adjust_c1

    @property
    def adjust_c4(self):
        return self._adjust_c4

    @property
    def scale_atten(self):
        return self._scale_atten

    def _calc_ln_resp(self, pga_ref):
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

        f_event = 1 if s.event_type == "intraslab" else 0
        dist = s.dist_hyp if s.event_type == "intraslab" else s.dist_rup

        path_atten = (c.t_2 + c.t_14 * f_event + c.t_3 * (s.mag - 7.8)) * np.log(
            dist + c.c_4 * np.exp(c.t_9 * (s.mag - 6))
        )

        ln_resp = (
            c.t_1
            + c.t_4 * self.adjust_c1
            + path_atten
            + self._scale_atten * c.t_6 * dist
            + c.t_10 * f_event
            + self._calc_f_mag(s.mag)
            + self._calc_f_site(pga_ref)
        )

        if s.tectonic_region == "backarc":
            ln_resp += self._calc_f_faba(dist)

        if s.event_type == "intraslab":
            ln_resp += self._calc_f_depth(s.depth_hyp)

        return ln_resp

    def _calc_f_mag(self, mag):
        """Calculate the magnitude scaling. Equation 2.

        Parameters
        ----------
        mag : float
            moment magnitude

        Returns
        -------
        f_mag : class:`np.array`:
            magnitude scaling

        """
        c = self.COEFF
        c_1 = 7.8
        coeff = np.select([mag <= (c_1 + self.adjust_c1), True], [c.t_4, c.t_5])
        f_mag = coeff * (mag - (c_1 + self.adjust_c1)) + c.t_13 * (10 - mag) ** 2
        return f_mag

    def _calc_f_depth(self, depth_hyp):
        """Calculate the depth scaling. Equation 3.

        Parameters
        ----------
        depth_hyp : float
            depth to hypocenter [km]

        Returns
        -------
        f_depth : class:`np.array`:
            depth scaling

        """
        c = self.COEFF
        f_depth = c.t_11 * (np.minimum(depth_hyp, 120) - 60)
        return f_depth

    def _calc_f_faba(self, dist):
        """Calculate the forearc/backarc scaling. Equation 4.

        Parameters
        ----------
        dist : float
            depth to hypocenter [km]

        Returns
        -------
        f_faba : class:`np.array`:
            forearc/backarc scaling

        """
        c = self.COEFF
        if self.scenario.event_type == "intraslab":
            f_faba = c.t_7 + c.t_8 * np.log(np.maximum(dist, 85) / 40)
        elif self.scenario.event_type == "interface":
            f_faba = c.t_15 + c.t_16 * np.log(np.maximum(dist, 100) / 40)
        else:
            raise NotImplementedError

        return f_faba

    def _calc_f_site(self, pga_ref):
        """Calculate the site amplification. Equation 5.

        Parameters
        ----------
        pga_ref : float, default: np.nan
            reference PGA (V_s30 = 1000 m/sec). If 'np.nan', then the ground
            motion at the reference condition is computed.

        Returns
        -------
        f_site : class:`np.array`:
            site scaling

        """
        c = self.COEFF
        v_s30 = 1000.0 if np.isnan(pga_ref) else self.scenario.v_s30
        vs_ratio = np.minimum(v_s30, 1000) / c.v_lin

        f_site = np.select(
            [v_s30 < c.v_lin, True],
            [
                c.t_12 * np.log(vs_ratio)
                - c.b * np.log(pga_ref + c.c)
                + c.b * np.log(pga_ref + c.c * vs_ratio ** c.n),
                (c.t_12 + c.b * c.n) * np.log(vs_ratio),
            ],
        )
        return f_site

    def _calc_ln_std(self):
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        c = self.COEFF

        ln_std = np.sqrt(c.phi ** 2 + c.tau ** 2)
        return ln_std
