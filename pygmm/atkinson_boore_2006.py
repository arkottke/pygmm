# -*- coding: utf-8 -*-
"""Atkinson and Boore (2006, :cite:`atkinson06`) model."""
import numpy as np

from . import model

__author__ = "Albert Kottke"


class AtkinsonBoore2006(model.GroundMotionModel):
    """Atkinson and Boore (2006, :cite:`atkinson06`) model.

    Developed for the Eastern North America with a reference velocity of 760
    or 2000 m/s.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Atkinson and Boore (2006)"
    ABBREV = "AB06"

    # Load the coefficients for the model
    COEFF = dict(
        bc=model.load_data_file("atkinson_boore_2006-bc.csv", 2),
        rock=model.load_data_file("atkinson_boore_2006-rock.csv", 2),
    )

    PERIODS = COEFF["bc"]["period"]

    COEFF_SITE = model.load_data_file("atkinson_boore_2006-site.csv", 2)
    COEFF_SF = model.load_data_file("atkinson_boore_2006-sf.csv", 2)

    INDEX_PGD = 0
    INDEX_PGV = 1
    INDEX_PGA = 2
    INDICES_PSA = np.arange(3, 27)

    PARAMS = [
        model.NumericParameter("mag", True),
        model.NumericParameter("dist_rup", True),
        model.NumericParameter("v_s30", True),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super(AtkinsonBoore2006, self).__init__(scenario)
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
        c = self.COEFF["bc"] if s.v_s30 else self.COEFF["rock"]

        # Compute the response at the reference condition
        r0 = 10.0
        r1 = 70.0
        r2 = 140.0

        f0 = np.maximum(np.log10(r0 / s.dist_rup), 0)
        f1 = np.minimum(np.log10(s.dist_rup), np.log10(r1))
        f2 = np.maximum(np.log10(s.dist_rup / r2), 0)

        # Compute the log10 PSA in units of cm/sec/sec
        log10_resp = (
            c.c_1
            + c.c_2 * s.mag
            + c.c_3 * s.mag ** 2
            + (c.c_4 + c.c_5 * s.mag) * f1
            + (c.c_6 + c.c_7 * s.mag) * f2
            + (c.c_8 + c.c_9 * s.mag) * f0
            + c.c_10 * s.dist_rup
        )

        # Apply stress drop correction
        log10_resp += self._calc_stress_factor()

        if s.v_s30:
            # Compute the site amplification
            pga_bc = 10 ** log10_resp[self.INDEX_PGA]

            log10_site = self._calc_log10_site(pga_bc)

            log10_resp += log10_site

        # Convert from cm/sec/sec to gravity
        log10_resp -= np.log10(980.665)

        ln_resp = np.log(10 ** log10_resp)
        return ln_resp

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        ln_std = np.ones_like(self.PERIODS) * 0.30
        return ln_std

    def _calc_stress_factor(self) -> float:
        """Calculate the stress correction factor proposed by Atkinson and
        Boore (2011) :cite:`atkinson11`.

        Returns
        -------
        log10_stress_factor : class:`np.array`:
            log base 10 of the stress factor
        """
        s = self._scenario
        c = self.COEFF_SF

        stress_drop = 10.0 ** (3.45 - 0.2 * s.mag)
        v1 = c.delta + 0.05
        v2 = 0.05 + c.delta * np.maximum(s.mag - c.m_1, 0) / (c.m_h - c.m_1)

        log10_stress_factor = np.minimum(2.0, stress_drop / 140.0) * np.minimum(v1, v2)

        return np.interp(self.PERIODS, c.period, log10_stress_factor)

    def _calc_log10_site(self, pga_bc: float) -> np.ndarray:
        """Calculate the log10 of the site amplification.

        Parameters
        ----------
        pga_bc : float
            peak ground acceleration (PGA, g) at the B/C boundary.

        Returns
        -------
        log_10_site : :class:`np.ndarray`
            log base 10 of the  site amplification.
        """
        s = self._scenario
        c = self.COEFF_SITE
        VS_1 = 180.0
        VS_2 = 300.0
        VS_REF = 760.0

        if s.v_s30 <= VS_1:
            b_nl = c.b_1
        elif VS_1 < s.v_s30 <= VS_2:
            b_nl = (c.b_1 - c.b_2) * np.log(s.v_s30 / VS_2) / np.log(VS_1 / VS_2)
        elif VS_2 < s.v_s30 <= VS_REF:
            b_nl = c.b_2 * np.log(s.v_s30 / VS_REF) / np.log(VS_2 / VS_REF)
        else:
            # Vs30 > VS_REF
            b_nl = 0

        pga_bc = max(pga_bc, 60.0)

        log10_site = np.log10(
            np.exp(c.b_lin * np.log(s.v_s30 / VS_REF) + b_nl * np.log(pga_bc / 100.0))
        )
        return np.interp(self.PERIODS, c.period, log10_site)
