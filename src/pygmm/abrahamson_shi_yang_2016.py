"""Abrahamson, Shi, and Yang (2016, :cite:`abrahamson16arias`) Arias intensity model."""

import numpy as np

from . import model

__author__ = "Albert Kottke"


class AbrahamsonShiYang2016(model.Model):
    """Abrahamson, Shi, and Yang (2016, :cite:`abrahamson16arias`) Arias intensity model.

    Conditional ground-motion model (CGMM) and scenario-based models for
    Arias intensity (Ia, m/s) for shallow crustal tectonic settings, consistent
    with the NGA-West2 ground-motion models. Calibrated on the NGA-West2
    database (~15,209 recordings, 326 events, M 3.0–7.9).

    Two evaluation modes are supported:

    **Conditional mode** (``backbone=None``): Ia is estimated from user-supplied
    PGA and SA(T=1 s) values (both in g) using Equations (3.3) and (3.4).
    The aleatory variability is
    :math:`\\sigma = \\sqrt{\\tau^2 + \\phi^2} \\approx 0.38`.

    **Scenario mode** (``backbone`` is one of ``'ASK14'``, ``'BSSA14'``,
    ``'CB14'``, ``'CY14'``, ``'I14'``): The median PGA and SA(T=1 s), and
    their lognormal standard deviations, are taken from the chosen NGA-West2
    GMM. These are then combined with the conditional model via propagation of
    errors (Equation 4.4) to yield a full scenario-based Ia estimate.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        Earthquake scenario. When using scenario mode the scenario must also
        contain all parameters required by the chosen backbone GMM.
    backbone : str or None, optional
        Backbone NGA-West2 GMM for scenario mode. Accepted values:
        ``'ASK14'``, ``'BSSA14'``, ``'CB14'``, ``'CY14'``, ``'I14'``.
        Default ``None`` selects conditional mode.
    pga : float or None, optional
        Peak ground acceleration (g) for conditional mode. Required when
        ``backbone`` is ``None``.
    sa_t1 : float or None, optional
        Spectral acceleration at T=1 s (g) for conditional mode. Required
        when ``backbone`` is ``None``.

    References
    ----------
    .. [1] Abrahamson, C., Shi, H.-J. M., and Yang, B. (2016). "Ground-Motion
           Prediction Equations for Arias Intensity Consistent with the
           NGA-West2 Ground-Motion Models." PEER Report No. 2016/05, Pacific
           Earthquake Engineering Research Center, University of California,
           Berkeley.

    """

    NAME = "Abrahamson, Shi, & Yang (2016)"
    ABBREV = "ASY16"

    # Regression coefficients — Table 3.1, Equation 3.3
    C1 = 0.47
    C2 = -0.28
    C3 = 0.50
    C4 = 1.52   # ln(PGA) scaling
    C5 = 0.21   # ln(SA_T1) scaling
    C8 = 0.09   # hanging-wall coefficient (Equation 3.4)

    # Aleatory variability — conditional model (Section 3.2)
    TAU = 0.15   # between-event standard deviation
    PHI = 0.35   # within-event standard deviation
    SIGMA_COND = float(np.sqrt(TAU**2 + PHI**2))

    # Error-propagation coefficients — Table 4.1 / Equation 4.4
    # (partial derivatives rounded as given in paper)
    C4_DERIV = 1.53
    C5_DERIV = 0.20
    # Baker & Jayaram (2008) correlation between ln(PGA) and ln(SA[T=1s])
    RHO_PGA_SAT1 = 0.52

    SUPPORTED_BACKBONE_MODELS = ("ASK14", "BSSA14", "CB14", "CY14", "I14")

    LIMITS = dict(mag=(3.0, 7.9), dist_rup=(0, 400), v_s30=(180, 2000))

    PARAMS = [
        model.NumericParameter("mag", True, 3.0, 7.9),
        model.NumericParameter("v_s30", True, 150, 1500),
        model.NumericParameter("dist_rup", False, 0, 400),
        model.NumericParameter("dist_jb", False, 0, 400),
        model.NumericParameter("dip", False, 15, 90, default=90.0),
        model.CategoricalParameter("on_hanging_wall", False, [True, False], False),
    ]

    def __init__(self, scenario, backbone=None, pga=None, sa_t1=None):
        """Initialize the model."""
        self._full_scenario = scenario
        super().__init__(scenario)
        s = self._scenario

        if backbone is not None:
            if backbone not in self.SUPPORTED_BACKBONE_MODELS:
                raise ValueError(
                    f"backbone '{backbone}' is not supported. "
                    f"Choose from: {', '.join(self.SUPPORTED_BACKBONE_MODELS)}"
                )
            gmm = self._build_backbone(backbone)
            pga_med = gmm.pga
            sa_t1_med = gmm.interp_spec_accels([1.0])[0]
            sigma_pga = gmm.ln_std_pga
            sigma_sa_t1 = gmm.interp_ln_stds([1.0])[0]
            self._ln_ia = self._calc_ln_ia(np.log(pga_med), np.log(sa_t1_med))
            self._ln_std = _calc_scenario_ln_std(sigma_pga, sigma_sa_t1)
            self._tau = None
            self._phi = None
        else:
            if pga is None or sa_t1 is None:
                raise ValueError(
                    "Either 'backbone' (scenario mode) or both 'pga' and "
                    "'sa_t1' (conditional mode) must be provided."
                )
            self._ln_ia = self._calc_ln_ia(np.log(float(pga)), np.log(float(sa_t1)))
            self._ln_std = self.SIGMA_COND
            self._tau = self.TAU
            self._phi = self.PHI

    # ------------------------------------------------------------------ #
    # Public properties                                                    #
    # ------------------------------------------------------------------ #

    @property
    def ln_ia(self):
        """Natural log of the median Arias intensity (ln(m/s))."""
        return self._ln_ia

    @property
    def ia(self):
        """Median Arias intensity (m/s)."""
        return np.exp(self._ln_ia)

    @property
    def ln_std(self):
        """Total lognormal standard deviation of Arias intensity."""
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
    def ia_plus_sigma(self):
        """84th-percentile Arias intensity (m/s): ``exp(ln_ia + ln_std)``."""
        return np.exp(self._ln_ia + self._ln_std)

    @property
    def ia_minus_sigma(self):
        """16th-percentile Arias intensity (m/s): ``exp(ln_ia - ln_std)``."""
        return np.exp(self._ln_ia - self._ln_std)

    # ------------------------------------------------------------------ #
    # Private helpers                                                      #
    # ------------------------------------------------------------------ #

    def _calc_ln_ia(self, ln_pga, ln_sa_t1):
        """Evaluate Equations (3.3) and (3.4) of Abrahamson et al. (2016)."""
        s = self._scenario
        return (
            self.C1
            + self.C2 * np.log(s.v_s30)
            + self.C3 * s.mag
            + self.C4 * ln_pga
            + self.C5 * ln_sa_t1
            + self._hanging_wall_term()
        )

    def _hanging_wall_term(self):
        """C8 · F_HW · T1(dip) · T2(Mw) · T5(RJB) — Equation (3.4)."""
        s = self._scenario
        if not s.on_hanging_wall:
            return 0.0
        rjb = s.dist_jb if s.dist_jb is not None else 999.0
        return self.C8 * _t1_dip(s.dip) * _t2_mag(s.mag) * _t5_rjb(rjb)

    def _build_backbone(self, backbone):
        """Lazily import and instantiate the chosen backbone NGA-West2 GMM."""
        from . import (
            abrahamson_silva_kamai_2014,
            boore_stewart_seyhan_atkinson_2014,
            campbell_bozorgnia_2014,
            chiou_youngs_2014,
            idriss_2014,
        )

        gmm_map = {
            "ASK14": abrahamson_silva_kamai_2014.AbrahamsonSilvaKamai2014,
            "BSSA14": boore_stewart_seyhan_atkinson_2014.BooreStewartSeyhanAtkinson2014,
            "CB14": campbell_bozorgnia_2014.CampbellBozorgnia2014,
            "CY14": chiou_youngs_2014.ChiouYoungs2014,
            "I14": idriss_2014.Idriss2014,
        }
        return gmm_map[backbone](self._full_scenario)


# ------------------------------------------------------------------ #
# Module-level taper functions (Equations 3.5–3.7)                   #
# ------------------------------------------------------------------ #

def _t1_dip(dip):
    """Dip taper T1 — Equation (3.5) of Abrahamson et al. (2016)."""
    if dip > 30:
        return (90.0 - dip) / 45.0
    return 60.0 / 45.0


def _t2_mag(mag):
    """Magnitude taper T2 — Equation (3.6) of Abrahamson et al. (2016)."""
    if mag >= 6.5:
        return 1.0
    if mag <= 5.5:
        return 0.0
    return 1.0 + 0.2 * (mag - 6.5) - 0.8 * (mag - 6.5) ** 2


def _t5_rjb(dist_jb):
    """Distance taper T5 — Equation (3.7) of Abrahamson et al. (2016)."""
    if dist_jb == 0:
        return 1.0
    if dist_jb < 15:
        return 1.0 - dist_jb / 15.0
    return 0.0


def _calc_scenario_ln_std(sigma_pga, sigma_sa_t1):
    """Total sigma for scenario mode — Equation (4.4) of Abrahamson et al. (2016).

    Parameters
    ----------
    sigma_pga : float
        Lognormal standard deviation of PGA from backbone GMM.
    sigma_sa_t1 : float
        Lognormal standard deviation of SA(T=1 s) from backbone GMM.

    Returns
    -------
    float
        Total lognormal standard deviation of Arias intensity.

    """
    C4D = AbrahamsonShiYang2016.C4_DERIV
    C5D = AbrahamsonShiYang2016.C5_DERIV
    rho = AbrahamsonShiYang2016.RHO_PGA_SAT1
    sigma_cond_sq = AbrahamsonShiYang2016.SIGMA_COND ** 2
    return np.sqrt(
        C4D**2 * sigma_pga**2
        + C5D**2 * sigma_sa_t1**2
        + 2 * rho * C4D * C5D * sigma_pga * sigma_sa_t1
        + sigma_cond_sq
    )
