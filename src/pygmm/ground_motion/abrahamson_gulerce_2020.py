"""Abrahamson and Gülerce (2020) ground motion model."""

import os

import numpy as np

from .. import model
from ..registry import register


def _load_ak_adjustment(periods: np.ndarray) -> np.ndarray:
    """Load the USGS Alaska interface adjustment for the model periods."""
    fname = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "nga_sub_ak_interface_adjustment.csv",
    )
    adjustments = {}
    with open(fname) as fp:
        lines = [line.strip() for line in fp if not line.startswith("#")]
    for line in lines[1:]:
        if not line:
            continue
        imt, value = line.split(",")
        if imt == "PGV":
            continue
        period = 0.0 if imt == "PGA" else float(imt)
        adjustments[period] = float(value)
    return np.array([adjustments[float(p)] for p in periods])


@register(tectonic="subduction")
class AbrahamsonGulerce2020(model.GroundMotionModel):
    r"""Abrahamson and Gülerce (2020) ground motion model.

    Regionalized NGA-Subduction ground motion model for subduction interface
    and intraslab events, developed as part of the PEER NGA-Subduction project.
    The model provides the RotD50 average horizontal component.

    The model is described in:

    * Abrahamson, N., and Gülerce, Z. (2020). Regionalized ground-motion models
      for subduction earthquakes based on the NGA-SUB database. PEER Report No.
      2020/25, Pacific Earthquake Engineering Research Center, Berkeley, CA.
    * Abrahamson, N. A., and Gülerce, Z. (2022). Summary of the Abrahamson and
      Gulerce NGA-SUB ground-motion model for subduction earthquakes.
      Earthquake Spectra, 38(4), 2638-2681. doi:10.1177/87552930221114374.

    The implementation follows ``AbrahamsonGulerce_2020`` in the USGS nshmp-lib
    code (commit 44728a7d), including its coefficients (``AG20.csv`` and, for
    Puerto Rico and the U.S. Virgin Islands, ``AG20_PRVI.csv``). As in
    nshmp-lib:

    * The global model and the Alaska, Cascadia, and Puerto Rico and U.S.
      Virgin Islands ("prvi") regional models are provided. The Alaska and
      Cascadia models use the regional constant (``a31`` and ``a32``), linear
      site (``a17`` and ``a18``), and anelastic attenuation (``a24`` and
      ``a25``) terms, and magnitude scaling break points for intraslab events
      of 7.9 and 7.1. The prvi model uses the global model with the constant
      and linear site coefficients from ``AG20_PRVI.csv``.
    * The 0.01 s coefficients are used for PGA.
    * There is no forearc-backarc term and no aftershock term.
    * ``depth_tor`` is only used for intraslab events (the depth term uses
      :math:`\min(Z_{tor}, 200) - 50`).
    * If ``depth_2_5`` is given, the Cascadia basin term
      (:math:`a_{39} \ln Z'`) is used for all regions and variants, but only
      when :math:`Z_{2.5}` is greater than the reference value computed from
      :math:`V_{s30}` (the term only amplifies). With ``basin=True`` (the
      ``*_BASIN`` variants), :math:`\ln Z'` is also multiplied by the USGS
      deep-basin scale factor, which is 0 for periods of 0.5 s or less and for
      :math:`Z_{2.5} \le 1` km, increases linearly from 0 at 1 km to 1 at 3
      km, and is multiplied by 0.585 at 0.75 s.
    * With ``adjusted=True`` (the ``*_ADJUSTED`` variants), the regional
      adjustment to the constant (``cAk`` for Alaska and ``cCasc`` for
      Cascadia) is added to the mean, including the reference PGA used by the
      nonlinear site term. It has no effect for the global and prvi models.
    * With ``ak_adjusted=True`` (``AG_20_GLOBAL_INTERFACE_AK_ADJUSTED``), the
      USGS Alaska data adjustment for the NGA-Subduction models is added to the
      mean (and the reference PGA) of interface events. It has no effect for
      intraslab events.
    * The standard deviation uses the linear :math:`\phi_1` model with
      :math:`\tau = 0.47` and the nonlinear site effects on :math:`\phi` and
      :math:`\tau`.

    The nshmp-lib Gmm ids map to the region and options as follows (the event
    type is "interface" for ``*_INTERFACE*`` and "intraslab" for
    ``*_SLAB*``):

    =============================================== =========== ==============
    nshmp-lib Gmm id                                region      options
    =============================================== =========== ==============
    ``AG_20_GLOBAL_{INTERFACE,SLAB}``               "global"
    ``AG_20_GLOBAL_{INTERFACE,SLAB}_NO_EPI``        "global"    epistemic=False
    ``AG_20_GLOBAL_INTERFACE_AK_ADJUSTED``          "global"    ak_adjusted=True
    ``AG_20_CASCADIA_{INTERFACE,SLAB}``             "cascadia"
    ``AG_20_CASCADIA_{INTERFACE,SLAB}_BASIN``       "cascadia"  basin=True
    ``AG_20_CASCADIA_{INTERFACE,SLAB}_ADJUSTED``    "cascadia"  adjusted=True
    ``AG_20_CASCADIA_{..}_ADJUSTED_BASIN``          "cascadia"  adjusted=True,
                                                                basin=True
    ``AG_20_ALASKA_{INTERFACE,SLAB}``               "alaska"
    ``AG_20_ALASKA_{INTERFACE,SLAB}_ADJUSTED``      "alaska"    adjusted=True
    ``AG_20_PRVI_{INTERFACE,SLAB}``                 "prvi"
    =============================================== =========== ==============

    With ``epistemic=True`` (default), the epistemic uncertainty of the median
    is represented by three branches with means :math:`\mu - 1.645
    \varepsilon`, :math:`\mu`, and :math:`\mu + 1.645 \varepsilon` (the 5th,
    50th, and 95th percentiles) and weights of 0.185, 0.63, and 0.185.
    :math:`\varepsilon` is 0.21 for Alaska, 0.27 for Cascadia, and, for the
    global and prvi models, 0.26 plus :math:`0.15 \ln(T / 3)` for periods
    longer than 3 s (N. Abrahamson, personal communication to the USGS,
    2023). As in nshmp-lib, the branches are collapsed to a single median,
    :math:`\ln \sum_i w_i \exp(\mu_i)`, and the standard deviation, which is
    the same for all branches. The collapsed values are not suitable for
    hazard calculations, which should sum the weighted exceedance
    probabilities of the branches, as nshmp-lib does; the branches are given
    by :meth:`~pygmm.model.GroundMotionModel.ln_branches`. With
    ``epistemic=False``, the central branch
    (:math:`\mu`) is used.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``depth_tor``, ``depth_2_5``, ``v_s30``, ``event_type``, and ``region``) can
    be a scalar or an array, and the arrays are broadcast against each other.
    For a scalar scenario, the response and standard deviation have one value
    per period, as in other models. For arrays of N scenarios, they have shape
    (N, periods), so, for example, ``pga`` has shape (N,). ``depth_tor`` is
    only needed for intraslab events, so it can be *None* if none of the events
    are intraslab events. ``depth_2_5`` can be *None* or *NaN* for no basin
    term. Results for events that are neither interface nor intraslab events
    are *NaN*. Invalid regions are replaced with "global" (with a warning).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    basin : bool, optional
        if *True*, the basin term is scaled by the USGS deep-basin scale factor
        (the ``*_BASIN`` variants). Default is *False*.
    adjusted : bool, optional
        if *True*, the Alaska or Cascadia adjustment to the constant is added
        (the ``*_ADJUSTED`` variants). Default is *False*.
    ak_adjusted : bool, optional
        if *True*, the USGS Alaska adjustment is added for interface events
        (``AG_20_GLOBAL_INTERFACE_AK_ADJUSTED``). Default is *False*.
    epistemic : bool, optional
        if *True* (default), the epistemic branches of the median are included.
        If *False*, only the central branch is used (the ``*_NO_EPI``
        variants).
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods, which are all of the periods), and/or "psa_all" (all 21
        periods). Computing only the needed intensity measures is much faster
        for large vectorized scenarios. If *None* (default), all intensity
        measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> from pygmm.model import Scenario
    >>> from pygmm.ground_motion.abrahamson_gulerce_2020 import AbrahamsonGulerce2020
    >>> s = Scenario(
    ...     mag=np.array([7.0, 9.0]), dist_rup=np.array([80.0, 100.0]),
    ...     depth_tor=np.array([50.0, np.nan]), depth_2_5=np.array([np.nan, 6.0]),
    ...     v_s30=400.0, event_type=np.array(["intraslab", "interface"]),
    ...     region="cascadia")
    >>> m = AbrahamsonGulerce2020(s, basin=True, adjusted=True)
    >>> m.pga.shape, m.spec_accels.shape
    ((2,), (2, 21))
    >>> AbrahamsonGulerce2020(s, epistemic=False, ims=["pga"]).pga.shape
    (2,)
    >>> m = AbrahamsonGulerce2020(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

    """

    NAME = "Abrahamson & Gülerce (2020)"
    ABBREV = "AG20"

    # Load the coefficients for the model
    COEFF = model.load_data_file("abrahamson_gulerce_2020.csv", 1)
    COEFF_PRVI = model.load_data_file("abrahamson_gulerce_2020-prvi.csv", 1)
    PERIODS = COEFF["period"]
    # USGS Alaska interface adjustment (nshmp-lib adj_ak)
    ADJ_AK = _load_ak_adjustment(PERIODS)

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 22)

    REGIONS = ["global", "alaska", "cascadia", "prvi"]

    # Reference shear-wave velocity (m/sec)
    V_REF = 1000.0
    # Maximum V_s30 used by the site term (m/sec)
    VSS_MAX = 1000.0
    # Constants
    C = 1.88
    N = 1.18
    C_4 = 10.0
    # Global magnitude scaling break point for intraslab events
    C_1_SLAB = 7.5
    # Regional magnitude scaling break points for intraslab events
    C_1_SLAB_ALASKA = 7.9
    C_1_SLAB_CASCADIA = 7.1
    # Between-event standard deviation
    TAU = 0.47
    # Variance of the site amplification
    PHI_AMP_SQ = 0.09
    # Epistemic uncertainty of the median
    EPI_GLOBAL = 0.26
    EPI_ALASKA = 0.21
    EPI_CASCADIA = 0.27
    # Weights of the 5th, 50th, and 95th percentile branches
    EPI_WEIGHTS = (0.185, 0.63, 0.185)
    EPI_Z_SCORE = 1.645

    LIMITS = dict(
        mag=(5.0, 9.5),
        dist_rup=(0.0, 1000.0),
        v_s30=(150.0, 1000.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 5.0, 9.5),
        model.NumericParameter("dist_rup", True, 0.0, 1000.0),
        # Only used for intraslab events
        model.NumericParameter("depth_tor", False),
        model.NumericParameter("depth_2_5", False),
        model.NumericParameter("v_s30", True, 150.0, 1000.0),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
        model.CategoricalParameter("region", False, REGIONS, "global"),
    ]

    def __init__(
        self,
        scenario: model.Scenario,
        basin: bool = False,
        adjusted: bool = False,
        ak_adjusted: bool = False,
        epistemic: bool = True,
        ims=None,
    ):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        basin : bool, optional
            if *True*, scale the basin term with the USGS deep-basin scale factor.
        adjusted : bool, optional
            if *True*, add the Alaska or Cascadia adjustment to the constant.
        ak_adjusted : bool, optional
            if *True*, add the USGS Alaska adjustment for interface events.
        epistemic : bool, optional
            if *True* (default), include the epistemic branches of the median.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed.
        """
        super().__init__(scenario, ims)
        self._basin = bool(basin)
        self._adjusted = bool(adjusted)
        self._ak_adjusted = bool(ak_adjusted)
        self._epistemic = bool(epistemic)
        self._epi_delta: np.ndarray | None = None
        self._ln_resp, self._ln_std = self._calc()
        if self._epistemic and self._epi_delta is not None:
            self._branches = model.symmetric_branches(
                self._ln_resp, self._ln_std, self._epi_delta, self.EPI_WEIGHTS
            )

    @property
    def basin(self) -> bool:
        """If the basin term is scaled by the USGS deep-basin scale factor."""
        return self._basin

    @property
    def adjusted(self) -> bool:
        """If the Alaska or Cascadia adjustment to the constant is added."""
        return self._adjusted

    @property
    def ak_adjusted(self) -> bool:
        """If the USGS Alaska adjustment is added for interface events."""
        return self._ak_adjusted

    @property
    def epistemic(self) -> bool:
        """If the epistemic branches of the median are included."""
        return self._epistemic

    def _regional(self, rows, region):
        """Region-dependent coefficients.

        Parameters
        ----------
        rows : slice or array_like
            coefficient rows, or *None* for all rows
        region : dict
            boolean masks of the Alaska, Cascadia, and prvi regions

        Returns
        -------
        c : :class:`numpy.recarray`
            coefficients
        r : dict
            regional coefficients
        """
        if rows is None:
            c = self.COEFF
            c_prvi = self.COEFF_PRVI
            adj_ak = self.ADJ_AK
        else:
            c = self.COEFF[rows]
            c_prvi = self.COEFF_PRVI[rows]
            adj_ak = self.ADJ_AK[rows]
        is_ak = region["alaska"]
        is_casc = region["cascadia"]
        is_prvi = region["prvi"]
        r = dict(
            a1=np.where(is_ak, c.a31, np.where(is_casc, c.a32, c.a1)),
            a1_adj=np.where(is_ak, c.cAk, np.where(is_casc, c.cCasc, 0.0)),
            a6=c.a6 + np.where(is_ak, c.a24, np.where(is_casc, c.a25, 0.0)),
            a12=c.a12 + np.where(is_ak, c.a17, np.where(is_casc, c.a18, 0.0)),
            c1s=np.where(
                is_ak,
                self.C_1_SLAB_ALASKA,
                np.where(is_casc, self.C_1_SLAB_CASCADIA, self.C_1_SLAB),
            ),
            adj_ak=adj_ak,
        )
        if np.any(is_prvi):
            r["a1"] = np.where(is_prvi, c_prvi.a1, r["a1"])
            r["a12"] = np.where(is_prvi, c_prvi.a12, r["a12"])
        return c, r

    def _calc(self):
        """Calculate the natural log of the response and the standard deviation.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        s = self._scenario
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        event_type = model.as_column(s.event_type)
        is_slab = model.equals(event_type, "intraslab")
        is_interface = model.equals(event_type, "interface")
        region_values = model.as_column(s.region)
        region = {
            name: model.equals(region_values, name)
            for name in ["alaska", "cascadia", "prvi"]
        }

        if s.depth_tor is None:
            if np.any(is_slab):
                raise ValueError("depth_tor is required for intraslab events")
            # Not used for interface events
            depth_tor = np.nan
        else:
            depth_tor = model.as_column(s.depth_tor)
        depth_2_5 = None if s.depth_2_5 is None else model.as_column(s.depth_2_5)

        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)
        v_s30 = model.as_column(s.v_s30)

        c, r = self._regional(self._indices, region)
        periods = self._coeff_rows(self.PERIODS)

        # Reference PGA (V_s30 = 1000 m/sec) without the basin term, which only
        # needs the PGA coefficients. It is only used by the nonlinear site
        # term (V_s30 < V_lin).
        c_pga, r_pga = self._regional([self.INDEX_PGA], region)
        ln_pga_rock = self._calc_ln_mean(
            c_pga, r_pga, is_slab, 0.0, mag, dist_rup, depth_tor, self.V_REF
        )
        pga_rock = np.where(v_s30 < c.vlin, np.exp(ln_pga_rock), 0.0)

        ln_resp = self._calc_ln_mean(
            c, r, is_slab, pga_rock, mag, dist_rup, depth_tor, v_s30
        )
        if depth_2_5 is not None:
            ln_resp = ln_resp + self._calc_f_basin(c, periods, v_s30, depth_2_5)

        if self._epistemic:
            # Collapse the branches: ln(sum_i w_i exp(mu + delta_i))
            #   = mu + ln(sum_i w_i exp(delta_i))
            epi = self.EPI_Z_SCORE * self._calc_epi(periods, region)
            self._epi_delta = epi
            w_lo, w_mid, w_hi = self.EPI_WEIGHTS
            ln_resp = ln_resp + np.log(w_lo * np.exp(-epi) + w_mid + w_hi * np.exp(epi))

        # The model is not defined for other event types
        ln_resp = np.where(is_slab | is_interface, ln_resp, np.nan)

        ln_std = self._calc_ln_std(c, c_pga, pga_rock, dist_rup, v_s30)
        ln_std = np.broadcast_to(ln_std, np.shape(ln_resp))
        return ln_resp, ln_std

    def _calc_ln_mean(
        self, c, r, is_slab, pga_rock, mag, dist_rup, depth_tor, v_s30
    ) -> np.ndarray:
        """Calculate the mean without the basin term and epistemic branches.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute
        r : dict
            regional coefficients for the periods to compute
        is_slab : array_like
            if the event is an intraslab event
        pga_rock : float or array_like
            peak ground acceleration [g] at the reference condition
            (V_s30 = 1000 m/sec). Only used if V_s30 < V_lin.
        mag : array_like
            moment magnitude
        dist_rup : array_like
            closest distance to the rupture [km]
        depth_tor : array_like
            depth to the top of the rupture [km], only used for intraslab
            events
        v_s30 : float or array_like
            time-averaged shear-wave velocity over the top 30 m [m/sec]

        Returns
        -------
        ln_mean : class:`np.array`:
            natural log of the response

        """
        # Finite fault term
        ln_dist = np.log(dist_rup + self.C_4 * np.exp((mag - 6.0) * c.a9))

        # Magnitude scaling
        c_1 = np.where(is_slab, r["c1s"], c.c1i)
        c_4s = c.a4 + np.where(is_slab, c.a45, 0.0)
        f_mag = np.where(mag <= c_1, c_4s, c.a5) * (mag - c_1) + c.a13 * (10 - mag) * (
            10 - mag
        )

        # Additional slab scaling and depth scaling (slab only)
        f_slab = np.where(
            is_slab,
            c.a10
            + (c.a4 + c.a45) * (r["c1s"] - self.C_1_SLAB)
            + c.a14 * ln_dist
            + np.where(depth_tor <= 50.0, c.a8, c.a11)
            * (np.minimum(depth_tor, 200.0) - 50.0),
            0.0,
        )

        # Site response
        vs_star = np.minimum(v_s30, self.VSS_MAX)
        vs_ratio = vs_star / c.vlin
        ln_vs = np.log(vs_ratio)
        f_site = r["a12"] * ln_vs + np.where(
            vs_star < c.vlin,
            -c.b * np.log(pga_rock + self.C)
            + c.b * np.log(pga_rock + self.C * vs_ratio**self.N),
            c.b * self.N * ln_vs,
        )

        ln_mean = (
            r["a1"]
            + (c.a2 + c.a3 * (mag - 7.0)) * ln_dist
            + r["a6"] * dist_rup
            + f_mag
            + f_site
            + f_slab
        )
        if self._adjusted:
            ln_mean = ln_mean + r["a1_adj"]
        if self._ak_adjusted:
            ln_mean = ln_mean + np.where(is_slab, 0.0, r["adj_ak"])
        return ln_mean

    def _calc_f_basin(self, c, periods, v_s30, depth_2_5) -> np.ndarray:
        """Calculate the basin term (Cascadia), which only amplifies."""
        ln_z_ref = np.where(
            v_s30 >= 570.0,
            7.6,
            np.where(v_s30 > 200.0, 8.52 - 0.88 * np.log(v_s30 / 200.0), 8.52),
        )
        ln_z_prime = np.log((depth_2_5 * 1000.0 + 50.0) / (np.exp(ln_z_ref) + 50.0))
        if self._basin:
            # USGS deep-basin scale factor for periods longer than 0.5 s
            scale = np.where(
                (periods > 0.5) & (depth_2_5 > 1.0),
                (np.clip(depth_2_5, 1.0, 3.0) - 1.0) / 2.0,
                0.0,
            )
            scale = np.where(periods == 0.75, scale * 0.585, scale)
            ln_z_prime = ln_z_prime * scale
        # NaN depths have no basin term
        return np.where(ln_z_prime > 0, c.a39 * ln_z_prime, 0.0)

    def _calc_epi(self, periods, region) -> np.ndarray:
        """Epistemic uncertainty of the median."""
        epi_global = self.EPI_GLOBAL + np.where(
            periods > 3.0, 0.15 * np.log(np.maximum(periods, 3.0) / 3.0), 0.0
        )
        return np.where(
            region["alaska"],
            self.EPI_ALASKA,
            np.where(region["cascadia"], self.EPI_CASCADIA, epi_global),
        )

    def _calc_ln_std(self, c, c_pga, pga_rock, dist_rup, v_s30) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        The epistemic branches all have the same standard deviation, so it does
        not depend on ``epistemic``.
        """
        tau_sq: float | np.ndarray = self.TAU**2
        phi_sq = self._calc_phi_lin_sq(c, dist_rup)

        # Nonlinear site effects. The partial derivative is 0 if pga_rock is 0
        # (V_s30 >= V_lin).
        phi_b_sq = phi_sq - self.PHI_AMP_SQ
        phi_b_sq_pga = self._calc_phi_lin_sq(c_pga, dist_rup) - self.PHI_AMP_SQ
        d_site = (
            c.b
            * pga_rock
            * (
                1.0 / (pga_rock + self.C * (v_s30 / c.vlin) ** self.N)
                - 1.0 / (pga_rock + self.C)
            )
        )
        nonlinear = np.minimum(v_s30, self.VSS_MAX) < c.vlin
        phi_sq = np.where(
            nonlinear,
            phi_sq
            + d_site**2 * phi_b_sq_pga
            + 2.0 * d_site * np.sqrt(phi_b_sq_pga) * np.sqrt(phi_b_sq) * c.rhoW,
            phi_sq,
        )
        tau_sq = np.where(
            nonlinear,
            tau_sq + d_site**2 * tau_sq + 2.0 * d_site * tau_sq * c.rhoB,
            tau_sq,
        )
        return np.sqrt(tau_sq + phi_sq)

    @staticmethod
    def _calc_phi_lin_sq(c, dist_rup) -> np.ndarray:
        """Linear within-event variance (phi_1^2)."""
        return c.d1 + c.d2 * (np.clip(dist_rup, 150.0, 450.0) - 150.0) / 300.0
