# -*- coding: utf-8 -*-
"""Model for the Campbell (2003) ground motion model."""
import numpy as np

from . import model

__author__ = "Albert Kottke"


class Campbell2003(model.GroundMotionModel):
    """Campbell (2003, :cite:`campbell03`) model.

    This model was developed for the Eastern US.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario


    """

    NAME = "Campbell (2003)"
    ABBREV = "C03"

    # Reference velocity (m/sec)
    V_REF = 2800.0

    COEFF = model.load_data_file("campbell_2003.csv", 1)
    PERIODS = COEFF["period"]

    INDICES_PSA = np.arange(16)

    PARAMS = [
        model.NumericParameter("mag", True, 5.0, 8.2),
        model.NumericParameter("dist_rup", True, None, 1000.0),
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
        c = self.COEFF
        s = self._scenario

        f_1 = c.c_2 * s.mag + c.c_3 * (8.5 - s.mag) ** 2
        # Distance scaling
        f_2 = (
            c.c_4
            * np.log(np.sqrt(s.dist_rup ** 2 + (c.c_7 * np.exp(c.c_8 * s.mag)) ** 2))
            + (c.c_5 + c.c_6 * s.mag) * s.dist_rup
        )
        # Geometric attenuation
        r_1 = 70.0
        r_2 = 130.0
        if s.dist_rup <= r_1:
            f_3 = 0.0
        else:
            f_3 = c.c_9 * (np.log(s.dist_rup) - np.log(r_1))

            if r_2 < s.dist_rup:
                f_3 += c.c_10 * (np.log(s.dist_rup) - np.log(r_2))

        # Compute the ground motion
        ln_resp = c.c_1 + f_1 + f_2 + f_3

        return ln_resp

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        c = self.COEFF
        s = self._scenario

        if s.mag < 7.16:
            ln_std = c.c_11 + c.c_12 * s.mag
        else:
            ln_std = c.c_13

        return ln_std
