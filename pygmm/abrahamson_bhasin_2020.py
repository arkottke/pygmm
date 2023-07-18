# -*- coding: utf-8 -*-
"""Abrahamson and Bhasin (2020, :cite:`abrahamson20`) model."""
from typing import Optional

import numpy as np

from . import model

__author__ = "Albert Kottke"


class AbrahamsonBhasin2020(model.Model):
    """Abrahamson (2020, :cite:`abrahamson20`) model.

    Conditional PGV model

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Abrahamson and Bhasin (2020)"
    ABBREV = "AB20"

    # Reference velocity (m/s)
    V_REF = 425.0

    # Load the coefficients for the model
    COEFFS = {
        # From Table 3.2
        "psa(T=Tpgv)": model.Coefficients(
            a_1=5.39,
            a_2=0.799,
            a_3=0.654,
            a_4=0.479,
            a_5=-0.062,
            a_6=-0.359,
            a_7=-0.134,
            a_8=0.023,
            phi=0.29,
            tau=0.16,
        ),
        "pga": model.Coefficients(
            a_1=4.77,
            a_2=0.738,
            a_3=0.484,
            a_4=0.275,
            a_5=-0.036,
            a_6=-0.332,
            a_7=-0.44,
            a_8=0,
            phi_1=0.32,
            phi_2=0.42,
            tau_1=0.12,
            tau_2=0.26,
            M_1=5.0,
            M_2=7.0,
        ),
        "psa(T=1s)": model.Coefficients(
            a_1=4.80,
            a_2=0.82,
            a_3=0.55,
            a_4=0.27,
            a_5=0.054,
            a_6=-0.382,
            a_7=-0.21,
            a_8=0.0,
            phi_1=0.28,
            phi_2=0.38,
            tau_1=0.12,
            tau_2=0.17,
            M_1=5.0,
            M_2=7.0,
        ),
    }

    PARAMS = [
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("mag", True, 3.0, 8.0),
        model.NumericParameter("v_s30", True, 180, 1500),
    ]

    def __init__(
        self,
        scenario: model.Scenario,
        psa=Optional[float],
        pga=Optional[float],
        psa_1s=Optional[float],
    ):
        super().__init__(scenario)

        if psa is not None:
            method = "psa(T=Tpgv)"
            # FIXME: How to interpolate?
            raise NotImplementedError
        elif pga is not None:
            method = "pga"
            ln_gm = np.log(pga)
        else:
            method = "psa(T=1s)"
            ln_gm = np.log(psa_1s)

        C = self.COEFFS[method]
        s = scenario

        # Mean
        if s.mag < 5:
            f_1 = C.a_2
        elif s.mag > 7.5:
            f_1 = C.a_3
        else:
            f_1 = C.a_2 + (C.a_3 - C.a_2) * (s.mag - 5) / 2.5

        ln_mean = (
            C.a_1
            + f_1 * ln_gm
            + C.a_4 * (s.mag - 6)
            + C.a_5 * (8.5 - s.mag) ** 2
            + C.a_6 * np.log(s.dist_rup + 5 * np.exp(0.4 * (s.mag - 6)))
            + (C.a_7 + C.a_8 * (s.mag - 5)) * np.log(C.v_s30 / self.V_REF)
        )

        if method == "psa(T=Tpgv)":
            phi, tau = C.phi, C.tau
        else:

            def interp(var_1, var_2):
                if s.mag < C.mag_1:
                    val = var_1
                elif s.mag > C.mag_2:
                    val = var_2
                else:
                    val = var_1 + (var_2 - var_1) * (s.mag - C.mag_1) / (
                        C.mag_2 - C.mag_1
                    )
                return val

            phi = interp(C.phi_1, C.phi_2)
            tau = interp(C.tau_1, C.tau_2)

        ln_std = np.sqrt(phi**2 + tau**2)
        return ln_mean, ln_std

    @classmethod
    def ln_period_pgv(cls, mag: float) -> float:
        """Natural logarithm of the period for calculating PGV. Equation 3.6.

        Parameters
        ----------
        mag : float
            Magnitude

        Returns
        -------
        ln_period_pgv : float
            Best period for computing PGV.
        """
        return -4.09 + 0.66 * mag
