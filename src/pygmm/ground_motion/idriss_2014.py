"""Idriss (2014, :cite:p:`idriss14`) model."""

import numpy as np

from .. import model
from ..registry import register

__author__ = "Albert Kottke"


@register(tectonic="active_crustal")
class Idriss2014(model.GroundMotionModel):
    """Idriss (2014, :cite:p:`idriss14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``v_s30``, and ``mechanism``) can be a scalar or an array, and the arrays
    are broadcast against each other. For a scalar scenario, the response and
    standard deviation have one value per period, as in other models. For
    arrays of N scenarios, the response and standard deviation have shape
    (N, periods), so, for example, ``pga`` has shape (N,) and
    ``spec_accels`` has shape (N, 22).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods), and/or "psa_all" (all 22 periods). Computing only
        the needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 50.0]),
    ...     v_s30=760.0, mechanism=np.array(["SS", "RS"]))
    >>> pygmm.Idriss2014(s).pga.shape
    (2,)
    >>> pygmm.Idriss2014(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.Idriss2014(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

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
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)
        v_s30 = model.as_column(s.v_s30)

        # Coefficients for small (M <= 6.75) or large magnitudes. Coefficients
        # that are the same for both are used directly, without selecting by
        # magnitude for each scenario.
        is_small = mag <= 6.75
        coeffs = {}
        for name in self.COEFF["small"].dtype.names or ():
            small = self._coeff_rows(self.COEFF["small"][name])
            large = self._coeff_rows(self.COEFF["large"][name])
            if np.array_equal(small, large):
                coeffs[name] = small
            else:
                coeffs[name] = np.where(is_small, small, large)
        c = model.Coefficients(**coeffs)

        # Only reverse earthquakes have a style of faulting term (SS/NS/U are 0)
        flag_mech = np.where(model.as_column(model.equals(s.mechanism, "RS")), 1, 0)

        f_mag = c.alpha_1 + c.alpha_2 * mag + c.alpha_3 * (8.5 - mag) ** 2
        f_dst = (
            -(c.beta_1 + c.beta_2 * mag) * np.log(dist_rup + 10) + c.gamma * dist_rup
        )
        f_ste = c.epsilon * np.log(v_s30)
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
            + 0.035 * np.log(np.clip(self._coeff_rows(self.PERIODS), 0.05, 3.0))
            - 0.06 * np.clip(model.as_column(s.mag), 5.0, 7.5)
        )
        return ln_std
