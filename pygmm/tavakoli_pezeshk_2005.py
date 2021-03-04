# -*- coding: utf-8 -*-
"""Tavakoli and Pezeshk (2005, :cite:`tavakoli05`) model."""
import numpy as np

from . import model

__author__ = "Albert Kottke"


class TavakoliPezeshk05(model.GroundMotionModel):
    """Tavakoli and Pezeshk (2005, :cite:`tavakoli05`) model.

    Developed for the Eastern North America with a reference velocity of 2880
    m/s.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Tavakoli and Pezeshk (2005)"
    ABBREV = "TP05"

    # Reference velocity (m/sec)
    V_REF = 2880.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("tavakoli_pezeshk_2005.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 14)

    PARAMS = [
        model.NumericParameter("dist_rup", True, None, 1000),
        model.NumericParameter("mag", True, 5.0, 8.2),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)
        self._ln_resp = self._calc_ln_resp()
        self._ln_std = self._calc_ln_std()

    def _calc_ln_resp(self) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response

        """
        s = self._scenario
        c = self.COEFF

        # Magnitude scaling
        f1 = c.c_1 + c.c_2 * s.mag + c.c_3 * (8.5 - s.mag) ** 2.5

        # Distance scaling
        f2 = c.c_9 * np.log(s.dist_rup + 4.5)

        if s.dist_rup > 70:
            f2 += c.c_10 * np.log(s.dist_rup / 70.0)

        if s.dist_rup > 130:
            f2 += c.c_11 * np.log(s.dist_rup / 130.0)

        # Calculate scaled, magnitude dependent distance R for use when
        # calculating f3
        dist = np.sqrt(
            s.dist_rup ** 2
            + (c.c_5 * np.exp(c.c_6 * s.mag + c.c_7 * (8.5 - s.mag) ** 2.5)) ** 2
        )
        f3 = (c.c_4 + c.c_13 * s.mag) * np.log(dist) + (c.c_8 + c.c_12 * s.mag) * dist

        # Compute the ground motion
        ln_resp = f1 + f2 + f3

        return ln_resp

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        s = self._scenario
        c = self.COEFF

        if s.mag < 7.2:
            ln_std = c.c_14 + c.c_15 * s.mag
        else:
            ln_std = c.c_16

        return ln_std
