"""Coppersmith and Bommer (2014) ground motion model, part of the PNNL Site
Wide Hazard."""
from __future__ import division

import numpy as np

from . import model
from .abrahamson_gregor_addo_2016 import AbrahamsonGregorAddo2016

__author__ = "Albert Kottke"


class CoppersmithBommer2014(AbrahamsonGregorAddo2016):
    """Coppersmith and Bommer (2014) ground motion model developed as part of
    the PNNL hazard study for the Hanford site.

    This model was developed for subduction regions and is a modified version
    of the BCHydro model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Coppersmith and Bommer (2014)"
    ABBREV = "CB14"

    # Load the coefficients for the model
    COEFF = model.load_data_file("coppersmith_bommer_2014.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDEX_PGV = None
    INDICES_PSA = np.arange(14)

    # FIXME
    LIMITS = dict(
        mag=(3.0, 8.5),
        dist_jb=(0.0, 300.0),
        v_s30=(150.0, 1500.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 3, 8.5),
        model.NumericParameter("dist_rup", False),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
    ]

    def __init__(self, scenario, scale_atten=1.0):
        """Initialize the model.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
        """
        super(CoppersmithBommer2014, self).__init__(
            scenario, adjust_c1=None, adjust_c4=0, scale_atten=scale_atten
        )

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

        path_atten = (c.t_2 + c.t_3 * (s.mag - 7.8)) * np.log(
            s.dist_rup + c.c_4 * np.exp(c.t_9 * (s.mag - 6))
        )

        ln_resp = (
            c.t_1
            + c.t_4 * self.adjust_c1
            + path_atten
            + self._scale_atten * c.t_6 * s.dist_rup
            + self._calc_f_mag(s.mag)
            + self._calc_f_site(pga_ref)
            + self._calc_f_faba(s.dist_rup)
        )

        return ln_resp

    def _calc_f_faba(self, dist):
        """Calculate the forearc/backarc scaling. Equation 4.

        Parameters
        ----------
        event_type : str
            type of subduction event, either: 'interface' or 'intraslab'

        dist : float
            depth to hypocenter [km]

        Returns
        -------
        f_faba : class:`np.array`:
            forearc/backarc scaling

        """
        c = self.COEFF
        f_faba = c.t_16 * np.log(np.maximum(dist, 40) / 40)

        return f_faba
