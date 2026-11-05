"""Model for the Campbell (2003) ground motion model."""

import numpy as np

from .. import model
from ..registry import register

__author__ = "Albert Kottke"


@register(tectonic="stable_continental")
class Campbell2003(model.GroundMotionModel):
    """Campbell (2003, :cite:p:`campbell03`) model.

    This model was developed for the Eastern US.

    The model is vectorized. Each scenario value (``mag`` and ``dist_rup``) can
    be a scalar or an array, and the arrays are broadcast against each other.
    For a scalar scenario, the response and standard deviation have one value
    per period, as in other models. For arrays of N scenarios, the response and
    standard deviation have shape (N, periods), so, for example,
    ``spec_accels`` has shape (N, 16).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: spectral periods such as "psa_1p000"
        (1.0 s) and/or "psa_all" (all 16 periods). The model does not provide
        PGA. Computing only the needed intensity measures is much faster for
        large vectorized scenarios. If *None* (default), all intensity
        measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 150.0]))
    >>> pygmm.Campbell2003(s).spec_accels.shape
    (2, 16)
    >>> m = pygmm.Campbell2003(s, ims=["psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

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
        c = self._coeff_rows(self.COEFF)
        s = self._scenario
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)

        f_1 = c.c_2 * mag + c.c_3 * (8.5 - mag) ** 2
        # Distance scaling
        f_2 = (
            c.c_4 * np.log(np.sqrt(dist_rup**2 + (c.c_7 * np.exp(c.c_8 * mag)) ** 2))
            + (c.c_5 + c.c_6 * mag) * dist_rup
        )
        # Geometric attenuation, which is zero within r_1 and has an additional
        # term beyond r_2. The distance is limited to r_1 in the logarithm to
        # avoid the logarithm of zero, where the term is not used.
        r_1 = 70.0
        r_2 = 130.0
        ln_dist = np.log(np.maximum(dist_rup, r_1))
        f_3 = np.where(dist_rup <= r_1, 0.0, c.c_9 * (ln_dist - np.log(r_1)))
        f_3 = f_3 + np.where(r_2 < dist_rup, c.c_10 * (ln_dist - np.log(r_2)), 0.0)

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
        c = self._coeff_rows(self.COEFF)
        mag = model.as_column(self._scenario.mag)

        ln_std = np.where(mag < 7.16, c.c_11 + c.c_12 * mag, c.c_13)

        return ln_std
