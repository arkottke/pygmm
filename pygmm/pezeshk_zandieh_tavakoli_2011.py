# -*- coding: utf-8 -*-
"""Pezeshk, Zandieh, and Tavakoli (2011, :cite:`pezeshk11`) model."""
import numpy as np

from . import model

__author__ = "Albert Kottke"


class PezeshkZandiehTavakoli2011(model.GroundMotionModel):
    """Pezeshk, Zandieh, and Tavakoli (2011, :cite:`pezeshk11`) model.

    Developed for the Eastern North America with a reference velocity of 2000
    m/s.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Pezeshk et al. (2011)"
    ABBREV = "Pea11"

    # Reference shear-wave velocity (m/sec)
    V_REF = 2000.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("pezeshk_zandieh_tavakoli_2011.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 23)

    PARAMS = [
        model.NumericParameter("mag", True, 5, 8),
        model.NumericParameter("dist_rup", True, None, 1000),
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

        dist = np.sqrt(s.dist_rup ** 2 + c.c_11 ** 2)

        log10_resp = (
            c.c_1
            + c.c_2 * s.mag
            + c.c_3 * s.mag ** 2
            + (c.c_4 + c.c_5 * s.mag) * np.minimum(np.log10(dist), np.log10(70.0))
            + (c.c_6 + c.c_7 * s.mag)
            * np.maximum(np.minimum(np.log10(dist / 70.0), np.log10(140.0 / 70.0)), 0.0)
            + (c.c_8 + c.c_9 * s.mag) * np.maximum(np.log10(dist / 140.0), 0)
            + c.c_10 * dist
        )

        ln_resp = np.log(np.power(10, log10_resp))
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

        if s.mag <= 7.0:
            ln_std_mean = c.c_12 * s.mag + c.c_13
        else:
            ln_std_mean = -6.95e-3 * s.mag + c.c_14

        ln_std = np.sqrt(ln_std_mean ** 2 + c["sigma_reg"] ** 2)

        return ln_std
