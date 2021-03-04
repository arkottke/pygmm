# -*- coding: utf-8 -*-
"""Idriss (2014, :cite:`idriss14`) model."""
import numpy as np

from . import model

__author__ = "Albert Kottke"


class Idriss2014(model.GroundMotionModel):
    """Idriss (2014, :cite:`idriss14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Idriss (2014)"
    ABBREV = "I14"

    # Reference velocity (m/s)
    V_REF = 1200.0

    # Load the coefficients for the model
    COEFF = dict(
        small=model.load_data_file("idriss_2014-small.csv", 2),
        large=model.load_data_file("idriss_2014-large.csv", 2),
    )
    PERIODS = COEFF["small"]["period"]

    INDEX_PGA = 0
    INDICES_PSA = np.arange(22)

    PARAMS = [
        model.NumericParameter("dist_rup", True, None, 150),
        model.NumericParameter("mag", True, 5, None),
        model.NumericParameter("v_s30", True, 450, 1200),
        model.CategoricalParameter("mechanism", True, ["SS", "RS"], "SS"),
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
        c = self.COEFF["small"] if s.mag <= 6.75 else self.COEFF["large"]

        if s.mechanism == "RS":
            flag_mech = 1
        else:
            # SS/RS/U
            flag_mech = 0

        f_mag = c.alpha_1 + c.alpha_2 * s.mag + c.alpha_3 * (8.5 - s.mag) ** 2
        f_dst = (
            -(c.beta_1 + c.beta_2 * s.mag) * np.log(s.dist_rup + 10)
            + c.gamma * s.dist_rup
        )
        f_ste = c.epsilon * np.log(s.v_s30)
        f_mec = c.phi * flag_mech

        ln_resp = f_mag + f_dst + f_ste + f_mec

        return ln_resp

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        s = self._scenario
        ln_std = (
            1.18
            + 0.035 * np.log(np.clip(self.PERIODS, 0.05, 3.0))
            - 0.06 * np.clip(s.mag, 5.0, 7.5)
        )
        return ln_std
