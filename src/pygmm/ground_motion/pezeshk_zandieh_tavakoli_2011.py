"""Pezeshk, Zandieh, and Tavakoli (2011, :cite:p:`pezeshk11`) model."""

import numpy as np

from .. import model
from ..registry import register

__author__ = "Albert Kottke"


@register(tectonic="stable_continental")
class PezeshkZandiehTavakoli2011(model.GroundMotionModel):
    """Pezeshk, Zandieh, and Tavakoli (2011, :cite:p:`pezeshk11`) model.

    Developed for the Eastern North America with a reference velocity of 2000
    m/s.

    The model is vectorized. Each scenario value (``mag`` and ``dist_rup``) can
    be a scalar or an array, and the arrays are broadcast against each other.
    For a scalar scenario, the response and standard deviation have one value
    per period, as in other models. For arrays of N scenarios, the response and
    standard deviation have shape (N, periods), so, for example, ``pga`` has
    shape (N,) and ``spec_accels`` has shape (N, 22).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 22 periods). Computing only
        the needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(mag=np.array([6.0, 7.5]), dist_rup=np.array([10.0, 150.0]))
    >>> pygmm.PezeshkZandiehTavakoli2011(s).pga.shape
    (2,)
    >>> pygmm.PezeshkZandiehTavakoli2011(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.PezeshkZandiehTavakoli2011(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

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

        dist = np.sqrt(model.as_column(s.dist_rup) ** 2 + c.c_11**2)

        log10_resp = (
            c.c_1
            + c.c_2 * mag
            + c.c_3 * mag**2
            + (c.c_4 + c.c_5 * mag) * np.minimum(np.log10(dist), np.log10(70.0))
            + (c.c_6 + c.c_7 * mag)
            * np.maximum(np.minimum(np.log10(dist / 70.0), np.log10(140.0 / 70.0)), 0.0)
            + (c.c_8 + c.c_9 * mag) * np.maximum(np.log10(dist / 140.0), 0)
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
        c = self._coeff_rows(self.COEFF)
        mag = model.as_column(self._scenario.mag)

        # The standard deviations are in log10 units
        log10_std_mean = np.where(
            mag <= 7.0, c.c_12 * mag + c.c_13, -6.95e-3 * mag + c.c_14
        )

        log10_std = np.sqrt(log10_std_mean**2 + c["sigma_reg"] ** 2)
        ln_std = np.log(10) * log10_std

        return ln_std
