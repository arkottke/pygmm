"""Macedo, Abrahamson, and Liu (2021, :cite:`macedo21`) CAV model."""

import numpy as np

from .. import model
from ..registry import register

__author__ = "Mahdi Bahrampouri"


@register(provides=("cav",), input="conditional")
class MacedoAbrahamsonLiu2021(model.Model):
    """Macedo, Abrahamson, and Liu (2021, :cite:`macedo21`) CAV model.

    Conditional ground-motion model (CGMM) and scenario-based models for
    cumulative absolute velocity (CAV, m/s) for shallow crustal tectonic
    settings. Calibrated on the NGA-West2 database (~14,000 recordings,
    287 events, M 3.0–7.9).

    Two evaluation modes are supported:

    **Conditional mode** (``pga_model=None``): CAV is estimated from a
    user-supplied PGA value (in g) using Eq. (12). The aleatory variability is
    :math:`\\sigma = \\sqrt{\\tau^2 + \\phi^2} \\approx 0.31`.

    **Scenario mode** (``pga_model`` is one of ``'ASK14'``, ``'BSSA14'``,
    ``'CB14'``, ``'CY14'``, ``'I14'``): The median PGA and its lognormal
    standard deviation are taken from the chosen NGA-West2 GMM, then combined
    with the conditional model via propagation of errors (Eq. 17) to yield a
    full scenario-based CAV estimate.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        Earthquake scenario. When using scenario mode the scenario must also
        contain all parameters required by the chosen backbone PGA GMM (e.g.
        ``dist_x``, ``dip``, ``depth_tor`` for CY14).
    pga_model : str or None, optional
        Backbone NGA-West2 GMM for scenario mode. Accepted values:
        ``'ASK14'``, ``'BSSA14'``, ``'CB14'``, ``'CY14'``, ``'I14'``.
        Default ``None`` selects conditional mode.
    pga : float or None, optional
        Peak ground acceleration in g for conditional mode. Required when
        ``pga_model`` is ``None``.

    References
    ----------
    .. [1] Macedo, J., Abrahamson, N., and Liu, C. (2021). "New Scenario-Based
           Cumulative Absolute Velocity Models for Shallow Crustal Tectonic
           Settings." Bulletin of the Seismological Society of America,
           111(1), 157–172. doi: 10.1785/0120190321.

    """

    NAME = "Macedo, Abrahamson, & Liu (2021)"
    ABBREV = "MAL21"

    # Regression coefficients from Table 1 of Macedo et al. (2021).
    # Eq. (12): ln(CAV[m/s]) = c1 + c2*ln(PGA[g]) + c3*Mw
    #                         + c4*ln(Vs30[m/s]) + c5*ln(Rrup[km])
    #                         + c6 * F_HW * T1(dip) * T2(Mw) * T5(RJB)
    C1 = 1.79
    C2 = 0.67
    C3 = 0.57
    C4 = -0.47
    C5 = -0.0026
    C6 = 0.17

    # Aleatory variability (Table 1).
    # Note: sqrt(0.17² + 0.26²) ≈ 0.31.
    TAU = 0.17  # between-event
    PHI = 0.26  # within-event
    SIGMA_COND = float(np.sqrt(TAU**2 + PHI**2))

    SUPPORTED_PGA_MODELS = ("ASK14", "BSSA14", "CB14", "CY14", "I14")

    PARAMS = [
        model.NumericParameter("mag", True, 3.0, 8.0),
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("dist_jb", True, 0, 300),
        model.NumericParameter("v_s30", True, 150, 1500),
        model.NumericParameter("dip", False, 15, 90, default=90.0),
        model.CategoricalParameter("on_hanging_wall", False, [True, False], False),
    ]

    def __init__(self, scenario, pga_model=None, pga=None):
        """Initialize the model."""
        # Preserve the full scenario so backbone GMMs can access their own params.
        self._full_scenario = scenario
        super().__init__(scenario)
        s = self._scenario

        if pga_model is not None:
            if pga_model not in self.SUPPORTED_PGA_MODELS:
                raise ValueError(
                    f"pga_model '{pga_model}' is not supported. "
                    f"Choose from: {', '.join(self.SUPPORTED_PGA_MODELS)}"
                )
            gmm = self._build_pga_gmm(pga_model)
            self._ln_cav = self._calc_ln_cav(np.log(gmm.pga), s)
            self._ln_std = _calc_scenario_ln_std(gmm.ln_std_pga)
            self._tau = None
            self._phi = None
        else:
            if pga is None:
                raise ValueError(
                    "Either 'pga_model' (scenario mode) or 'pga' "
                    "(conditional mode) must be provided."
                )
            self._ln_cav = self._calc_ln_cav(np.log(float(pga)), s)
            self._ln_std = self.SIGMA_COND
            self._tau = self.TAU
            self._phi = self.PHI

    # ------------------------------------------------------------------ #
    # Public properties                                                    #
    # ------------------------------------------------------------------ #

    @property
    def ln_cav(self):
        """Natural log of the median CAV (m/s)."""
        return self._ln_cav

    @property
    def cav(self):
        """Median CAV (m/s)."""
        return np.exp(self._ln_cav)

    @property
    def ln_std(self):
        """Total lognormal standard deviation of CAV."""
        return self._ln_std

    @property
    def tau(self):
        """Between-event standard deviation (conditional mode only)."""
        return self._tau

    @property
    def phi(self):
        """Within-event standard deviation (conditional mode only)."""
        return self._phi

    @property
    def cav_plus_sigma(self):
        """84th-percentile CAV (m/s): ``exp(ln_cav + ln_std)``."""
        return np.exp(self._ln_cav + self._ln_std)

    @property
    def cav_minus_sigma(self):
        """16th-percentile CAV (m/s): ``exp(ln_cav - ln_std)``."""
        return np.exp(self._ln_cav - self._ln_std)

    # ------------------------------------------------------------------ #
    # Private helpers                                                      #
    # ------------------------------------------------------------------ #

    def _calc_ln_cav(self, ln_pga, s):
        """Evaluate Eq. (12) / Eq. (14) of Macedo et al. (2021)."""
        hw = self._hanging_wall_term(s)
        return (
            self.C1
            + self.C2 * ln_pga
            + self.C3 * s.mag
            + self.C4 * np.log(s.v_s30)
            + self.C5 * np.log(s.dist_rup)
            + self.C6 * hw
        )

    def _hanging_wall_term(self, s):
        """F_HW * T1(dip) * T2(Mw) * T5(RJB) — last term of Eq. (12)."""
        if not s.on_hanging_wall:
            return 0.0
        return _t1_dip(s.dip) * _t2_mag(s.mag) * _t5_rjb(s.dist_jb)

    def _build_pga_gmm(self, pga_model):
        """Lazily import and instantiate the chosen backbone PGA GMM."""
        # Imported from the package namespace rather than by module path, so
        # this does not need updating when a model moves between subpackages.
        from .. import (
            AbrahamsonSilvaKamai2014,
            BooreStewartSeyhanAtkinson2014,
            CampbellBozorgnia2014,
            ChiouYoungs2014,
            Idriss2014,
        )

        gmm_map = {
            "ASK14": AbrahamsonSilvaKamai2014,
            "BSSA14": BooreStewartSeyhanAtkinson2014,
            "CB14": CampbellBozorgnia2014,
            "CY14": ChiouYoungs2014,
            "I14": Idriss2014,
        }
        return gmm_map[pga_model](self._full_scenario)


# ------------------------------------------------------------------ #
# Module-level taper functions (Eqs. 13a–c)                          #
# ------------------------------------------------------------------ #


def _t1_dip(dip):
    """Dip taper T1 — Eq. (13a) of Macedo et al. (2021)."""
    if dip > 30:
        return (90.0 - dip) / 45.0
    return 60.0 / 45.0


def _t2_mag(mag):
    """Magnitude taper T2 — Eq. (13b) of Macedo et al. (2021)."""
    if mag >= 6.5:
        return 1.0
    elif mag <= 5.5:
        return 0.0
    return 1.0 + 0.2 * (mag - 6.5) - 0.8 * (mag - 6.5) ** 2


def _t5_rjb(dist_jb):
    """Distance taper T5 — Eq. (13c) of Macedo et al. (2021)."""
    if dist_jb == 0:
        return 1.0
    elif dist_jb < 15:
        return 1.0 - dist_jb / 15.0
    return 0.0


def _calc_scenario_ln_std(sigma_ln_pga):
    """Total sigma for scenario mode — Eq. (17) of Macedo et al. (2021)."""
    return np.sqrt(0.097 + 0.4541 * sigma_ln_pga**2)
