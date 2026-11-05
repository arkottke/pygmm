"""Tavakoli and Pezeshk (2005, :cite:p:`tavakoli05`) model."""

import numpy as np

from .. import model
from ..registry import register

__author__ = "Albert Kottke"


@register(tectonic="stable_continental")
class TavakoliPezeshk05(model.GroundMotionModel):
    """Tavakoli and Pezeshk (2005, :cite:p:`tavakoli05`) model.

    Developed for the Eastern North America with a reference velocity of 2880
    m/s.

    The model is vectorized. Each scenario value (``mag`` and ``dist_rup``) can
    be a scalar or an array, and the arrays are broadcast against each other.
    For a scalar scenario, the response and standard deviation have one value
    per period, as in other models. For arrays of N scenarios, the response and
    standard deviation have shape (N, periods), so, for example, ``pga`` has
    shape (N,) and ``spec_accels`` has shape (N, 13).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 13 periods). Computing only
        the needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 150.0]))
    >>> pygmm.TavakoliPezeshk05(s).pga.shape
    (2,)
    >>> pygmm.TavakoliPezeshk05(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.TavakoliPezeshk05(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

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

    def __init__(self, scenario: model.Scenario, ims=None):
        """Initialize the model."""
        super().__init__(scenario, ims)
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
        c = self._coeff_rows(self.COEFF)
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)

        # Magnitude scaling
        f1 = c.c_1 + c.c_2 * mag + c.c_3 * (8.5 - mag) ** 2.5

        # Distance scaling, with additional terms beyond 70 and 130 km. The
        # distance is limited to 70 km in those terms to avoid the logarithm of
        # zero, where the terms are not used.
        f2 = c.c_9 * np.log(dist_rup + 4.5)
        dist_far = np.maximum(dist_rup, 70.0)
        f2 = f2 + np.where(dist_rup > 70, c.c_10 * np.log(dist_far / 70.0), 0.0)
        f2 = f2 + np.where(dist_rup > 130, c.c_11 * np.log(dist_far / 130.0), 0.0)

        # Calculate scaled, magnitude dependent distance R for use when
        # calculating f3
        dist = np.sqrt(
            dist_rup**2
            + (c.c_5 * np.exp(c.c_6 * mag + c.c_7 * (8.5 - mag) ** 2.5)) ** 2
        )
        f3 = (c.c_4 + c.c_13 * mag) * np.log(dist) + (c.c_8 + c.c_12 * mag) * dist

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
        c = self._coeff_rows(self.COEFF)
        mag = model.as_column(self._scenario.mag)

        ln_std = np.where(mag < 7.2, c.c_14 + c.c_15 * mag, c.c_16)

        return ln_std
