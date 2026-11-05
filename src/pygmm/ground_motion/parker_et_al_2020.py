"""Parker et al. (2020) NGA-Subduction ground motion model."""

import os
from typing import Any

import numpy as np

from .. import model
from ..registry import register


def _erf(x):
    """Approximate error function used by nshmp-lib (``Maths.erf``).

    This is the Abramowitz and Stegun (1964) equation 7.1.26 approximation
    (absolute error less than 1.5e-7), which is used instead of the exact error
    function to reproduce nshmp-lib.
    """
    ax = np.abs(x)
    t = 1 / (1 + 0.3275911 * ax)
    tsq = t * t
    value = 1 - (
        0.254829592 * t
        + -0.284496736 * tsq
        + 1.421413741 * tsq * t
        + -1.453152027 * tsq * tsq
        + 1.061405429 * tsq * tsq * t
    ) * np.exp(-ax * ax)
    return np.where(x < 0.0, -value, value)


def _load_ak_adjustment(periods):
    """Load the shared nshmp-lib Alaska interface adjustment by period."""
    fname = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "nga_sub_ak_interface_adjustment.csv",
    )
    data = np.genfromtxt(
        fname, skip_header=1, delimiter=",", names=True, dtype=None, encoding="utf-8"
    )
    by_period = {}
    for name, value in zip(data["T"], data["adj_ak"]):
        period = {"PGA": 0.0, "PGV": -1.0}.get(name)
        by_period[float(name) if period is None else period] = float(value)
    return np.array([by_period[float(p)] for p in periods])


@register(tectonic="subduction")
class ParkerEtAl2020(model.GroundMotionModel):
    r"""Parker et al. (2020) NGA-Subduction ground motion model.

    Ground motion model for subduction interface and intraslab earthquakes
    developed as part of the NGA-Subduction project, with the global model and
    the Alaska, Cascadia, and Puerto Rico and U.S. Virgin Islands (PRVI)
    regionalizations used by the U.S. Geological Survey. The component is
    RotD50.

    References: Parker, G. A., Stewart, J. P., Boore, D. M., Atkinson, G. M.,
    and Hassani, B. (2020). NGA-Subduction global ground motion models with
    regional adjustment factors. PEER Report No. 2020/03, Pacific Earthquake
    Engineering Research Center, Berkeley, CA. Parker, G. A., Stewart, J. P.,
    Boore, D. M., Atkinson, G. M., and Hassani, B. (2022). NGA-Subduction
    global ground motion models with regional adjustment factors. Earthquake
    Spectra, 38(1), 456-493. https://doi.org/10.1177/87552930211034889

    The implementation follows the ``ParkerEtAl_2020`` class in the USGS
    nshmp-lib code (commit 44728a7d), including its coefficients
    (``PSBAH20_interface.csv``, ``PSBAH20_slab.csv``, and ``PSBAH20_PRVI.csv``)
    and the shared nshmp-lib Alaska interface adjustment. As in nshmp-lib:

    * The PRVI model uses the PRVI constant (:math:`c_0`) and linear site
      coefficient (:math:`s_2`), and the global values of the other regional
      coefficients, corner magnitudes, and epistemic uncertainty.
    * For intraslab events, the hypocentral depth used by the source depth
      term is estimated from the depth to the top of rupture,
      :math:`Z_{hyp} = Z_{tor} + 0.48 \times 6.5` km.
    * The reference PGA (:math:`V_{s30}` = 760 m/sec) used by the nonlinear
      site term excludes the site and basin terms, and includes the Alaska
      adjustment if it is applied.
    * If ``depth_2_5`` is given (not *NaN*), the Cascadia basin term
      (equations 11 to 13 with the Cascadia :math:`Z_{2.5}` model) is applied
      for all regions and event types. Without ``basin``, it is applied to all
      periods without the USGS scaling.
    * The standard deviation is :math:`\sqrt{\tau^2 + \phi^2}`, with the
      distance and :math:`V_{s30}` dependent :math:`\phi` (equations 18 to
      20).
    * With ``epistemic=True`` (default), the epistemic uncertainty (equation
      27 and Table E4) is represented by three branches with means
      :math:`\mu - 1.645 \sigma_\epsilon`, :math:`\mu`, and :math:`\mu +
      1.645 \sigma_\epsilon` (the 5th, 50th, and 95th percentiles) and
      weights of 0.185, 0.63, and 0.185. As in nshmp-lib
      (``GroundMotions.combine``), the branches are collapsed to a single
      median, :math:`\ln \sum_i w_i \exp(\mu_i)`, and the standard
      deviation, which is the same for all branches. The collapsed values are
      not suitable for hazard calculations, which should sum the weighted
      exceedance probabilities of the branches, as nshmp-lib does; the
      branches are given by :meth:`~pygmm.model.GroundMotionModel.ln_branches`.
      With ``epistemic=False``, the central branch is used (as in the
      nshmp-lib ``*_NO_EPI`` variants).

    The nshmp-lib ``Gmm`` ids correspond to the ``event_type`` and ``region``
    scenario values and the model options as follows:

        ===================================== ========== ======== ===================
        nshmp-lib Gmm id                      event_type region   options
        ===================================== ========== ======== ===================
        PSBAH_20_GLOBAL_INTERFACE             interface  global
        PSBAH_20_GLOBAL_INTERFACE_NO_EPI      interface  global   epistemic=False
        PSBAH_20_GLOBAL_INTERFACE_AK_ADJUSTED interface  global   ak_adjusted=True
        PSBAH_20_GLOBAL_SLAB                  intraslab  global
        PSBAH_20_GLOBAL_SLAB_NO_EPI           intraslab  global   epistemic=False
        PSBAH_20_ALASKA_INTERFACE             interface  alaska
        PSBAH_20_ALASKA_SLAB                  intraslab  alaska
        PSBAH_20_CASCADIA_INTERFACE           interface  cascadia
        PSBAH_20_CASCADIA_INTERFACE_BASIN     interface  cascadia basin=True
        PSBAH_20_CASCADIA_INTERFACE_BASIN_M9  interface  cascadia basin=True, m9=True
        PSBAH_20_CASCADIA_SLAB                intraslab  cascadia
        PSBAH_20_CASCADIA_SLAB_BASIN          intraslab  cascadia basin=True
        PSBAH_20_PRVI_INTERFACE               interface  prvi
        PSBAH_20_PRVI_SLAB                    intraslab  prvi
        ===================================== ========== ======== ===================

    nshmp-lib only provides the ``basin``, ``m9``, and ``ak_adjusted`` options
    for some regions and event types, so they only affect those scenarios:
    ``basin`` affects Cascadia events, ``m9`` affects Cascadia interface
    events, and ``ak_adjusted`` affects global interface events.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``depth_tor``, ``depth_2_5``, ``v_s30``, ``event_type``, and ``region``)
    can be a scalar or an array, and the arrays are broadcast against each
    other, so one array can mix event types and regions. For a scalar
    scenario, the response and standard deviation have one value per period,
    as in other models. For arrays of N scenarios, they have shape (N,
    periods), so, for example, ``pga`` has shape (N,). ``depth_tor`` is only
    needed for intraslab events, so it can be *None* if none of the events are
    intraslab events. ``depth_2_5`` can be *None* (or *NaN* for some
    scenarios) for no basin term. Results for events that are neither
    interface nor intraslab events are *NaN*.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    basin : bool, optional
        if *True*, the USGS scaling of the basin term is applied to Cascadia
        events (the ``*_BASIN`` variants): the basin term is only used for
        periods longer than 0.5 s and is scaled from zero at :math:`Z_{2.5}`
        of 1 km to the full term at 3 km (and by 0.585 at 0.75 s). Default is
        *False*.
    ak_adjusted : bool, optional
        if *True*, the USGS Alaska data adjustment is added to global
        interface events (``PSBAH_20_GLOBAL_INTERFACE_AK_ADJUSTED``). Default
        is *False*.
    m9 : bool, optional
        if *True*, the USGS adjustment based on the M9 Cascadia simulations is
        applied to Cascadia interface events at sites with :math:`Z_{2.5}`
        greater than 6 km and periods of 2 s and longer
        (``PSBAH_20_CASCADIA_INTERFACE_BASIN_M9``). It requires
        ``basin=True``. Default is *False*.
    epistemic : bool, optional
        if *True* (default), the epistemic branches are included. If *False*,
        only the central branch is used.
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods, which are all of the periods), and/or "psa_all". Computing
        only the needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> from pygmm.model import Scenario
    >>> from pygmm.ground_motion.parker_et_al_2020 import ParkerEtAl2020
    >>> s = Scenario(
    ...     mag=np.array([7.0, 9.0]), dist_rup=np.array([80.0, 100.0]),
    ...     depth_tor=50.0, depth_2_5=np.array([np.nan, 4.0]), v_s30=400.0,
    ...     event_type=np.array(["intraslab", "interface"]), region="cascadia")
    >>> m = ParkerEtAl2020(s, basin=True)
    >>> m.pga.shape, m.spec_accels.shape
    ((2,), (2, 21))
    >>> ParkerEtAl2020(s, basin=True, ims=["pga"]).pga.shape
    (2,)
    >>> m = ParkerEtAl2020(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

    """

    NAME = "Parker et al. (2020)"
    ABBREV = "PSBAH20"

    # Load the coefficients for the model
    COEFF_INTERFACE = model.load_data_file("parker_et_al_2020-interface.csv", 1)
    COEFF_SLAB = model.load_data_file("parker_et_al_2020-slab.csv", 1)
    COEFF_PRVI = model.load_data_file("parker_et_al_2020-prvi.csv", 1)
    # Period (s), with 0 for PGA and -1 for PGV
    PERIODS = COEFF_INTERFACE["period"]
    # Alaska interface adjustment shared by the NGA-Subduction models
    ADJ_AK = _load_ak_adjustment(PERIODS)

    INDEX_PGV = 0
    INDEX_PGA = 1
    INDICES_PSA = np.arange(2, 23)

    REGIONS = ["global", "alaska", "cascadia", "prvi"]

    # Constants
    B_4 = 0.1
    V_1 = 270.0
    V_REF = 760.0
    # Nonlinear site term
    F_3 = 0.05
    # Source depth term (intraslab events)
    D_B1 = 20.0
    D_B2 = 67.0
    THETA_W = 0.48
    Z_SLAB = 6.5
    # Cascadia basin depth model
    THETA_0 = 3.94
    THETA_1 = -0.42
    NU_MU = 200.0
    NU_SIGMA = 0.2
    # Reference Z_2.5 (m) for the M9 adjustment
    Z_2_5_REF_M9 = 1279.0
    # Phi model
    PHI_R1 = 200.0
    PHI_R2 = 500.0
    PHI_V1 = 200.0
    PHI_V2 = 500.0
    # Weights of the 5th, 50th, and 95th percentile epistemic branches and the
    # z-score of the 5th and 95th percentiles
    EPI_WEIGHTS = (0.185, 0.63, 0.185)
    EPI_Z_SCORE = 1.645

    # Corner magnitudes, (interface, intraslab), Table 4.1
    M_C = dict(
        global_=(7.9, 7.6),
        alaska=(8.6, 7.2),
        cascadia=(7.7, 7.2),
    )
    # Epistemic uncertainty (sigma_eps_1, sigma_eps_2, t_1, t_2) for
    # (interface, intraslab), Table E4
    EPI = dict(
        global_=((0.40, 0.40, 0.20, 0.40), (0.35, 0.22, 0.15, 2.00)),
        alaska=((0.15, 0.10, 1.00, 4.00), (0.15, 0.12, 0.50, 1.00)),
        cascadia=((0.43, 0.33, 0.20, 0.50), (0.35, 0.16, 0.20, 3.00)),
    )

    LIMITS = dict(
        mag=(4.5, 9.5),
        dist_rup=(0.0, 1000.0),
        v_s30=(150.0, 2000.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 4.5, 9.5),
        model.NumericParameter("dist_rup", True, 0.0, 1000.0),
        # Only used for intraslab events
        model.NumericParameter("depth_tor", False),
        model.NumericParameter("depth_2_5", False),
        model.NumericParameter("v_s30", True, 150.0, 2000.0),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
        model.CategoricalParameter("region", False, REGIONS, "global"),
    ]

    def __init__(
        self,
        scenario: model.Scenario,
        basin: bool = False,
        ak_adjusted: bool = False,
        m9: bool = False,
        epistemic: bool = True,
        ims=None,
    ):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        basin : bool, optional
            if *True*, apply the USGS scaling of the basin term for Cascadia events.
        ak_adjusted : bool, optional
            if *True*, add the USGS Alaska adjustment to global interface events.
        m9 : bool, optional
            if *True*, apply the USGS M9 adjustment for Cascadia interface events.
            Requires ``basin=True``.
        epistemic : bool, optional
            if *True* (default), include the epistemic branches.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed.
        """
        if m9 and not basin:
            raise ValueError(
                "m9=True requires basin=True "
                "(nshmp-lib PSBAH_20_CASCADIA_INTERFACE_BASIN_M9)"
            )
        super().__init__(scenario, ims)
        self._basin = bool(basin)
        self._ak_adjusted = bool(ak_adjusted)
        self._m9 = bool(m9)
        self._epistemic = bool(epistemic)
        self._epi_delta: np.ndarray | None = None
        self._ln_resp, self._ln_std = self._calc()
        if self._epistemic and self._epi_delta is not None:
            self._branches = model.symmetric_branches(
                self._ln_resp, self._ln_std, self._epi_delta, self.EPI_WEIGHTS
            )

    @property
    def basin(self) -> bool:
        """If the USGS scaling of the basin term is applied."""
        return self._basin

    @property
    def ak_adjusted(self) -> bool:
        """If the USGS Alaska adjustment is applied to global interface events."""
        return self._ak_adjusted

    @property
    def m9(self) -> bool:
        """If the USGS M9 adjustment is applied to Cascadia interface events."""
        return self._m9

    @property
    def epistemic(self) -> bool:
        """If the epistemic branches are included."""
        return self._epistemic

    def _scenario_values(self) -> dict:
        """Scenario values with a trailing axis for the periods."""
        s = self._scenario
        event_type = model.as_column(s.event_type)
        region = model.as_column(s.region)
        v: dict[str, Any] = dict(
            is_slab=model.equals(event_type, "intraslab"),
            is_interface=model.equals(event_type, "interface"),
            is_alaska=model.equals(region, "alaska"),
            is_cascadia=model.equals(region, "cascadia"),
            is_prvi=model.equals(region, "prvi"),
            mag=model.as_column(s.mag),
            dist_rup=model.as_column(s.dist_rup),
            v_s30=model.as_column(s.v_s30),
        )
        if s.depth_tor is None:
            if np.any(v["is_slab"]):
                raise ValueError("depth_tor is required for intraslab events")
            # Not used for interface events
            v["depth_tor"] = np.nan
        else:
            v["depth_tor"] = model.as_column(s.depth_tor).astype(float)
        if s.depth_2_5 is None:
            v["depth_2_5"] = np.nan
        else:
            v["depth_2_5"] = model.as_column(s.depth_2_5).astype(float)
        return v

    def _by_region(self, v, values) -> np.ndarray:
        """Select (global, alaska, cascadia, prvi) values for each scenario."""
        global_, alaska, cascadia, prvi = values
        return np.select(
            [v["is_alaska"], v["is_cascadia"], v["is_prvi"]],
            [alaska, cascadia, prvi],
            global_,
        )

    def _constants(self, v, table) -> list:
        """Select region and event type dependent constants (PRVI is global)."""
        values = []
        for i in range(len(table["global_"][0])):

            def pick(region):
                interface, slab = table[region]
                return np.where(v["is_slab"], slab[i], interface[i])

            g = pick("global_")
            values.append(self._by_region(v, (g, pick("alaska"), pick("cascadia"), g)))
        return values

    def _calc_ln_mean_ref(self, v, rows) -> tuple:
        """Mean at the reference condition (V_s30 = 760 m/sec, no basin).

        Equation 1 without the site term.

        Returns
        -------
        ln_mean : :class:`np.ndarray`
            natural log of the response at the reference condition
        c : dict
            coefficients for each scenario and period

        """
        is_slab = v["is_slab"]
        ci = self.COEFF_INTERFACE[rows]
        cs = self.COEFF_SLAB[rows]
        cp = self.COEFF_PRVI[rows]

        def by_type(name):
            return np.where(is_slab, cs[name], ci[name])

        def by_region(name, prvi):
            return self._by_region(
                v,
                (
                    by_type("Global_" + name),
                    by_type("Alaska_" + name),
                    by_type("Cascadia_" + name),
                    prvi,
                ),
            )

        c = {
            name: by_type(name)
            for name in [
                "c1",
                "c4",
                "c5",
                "c6",
                "V2",
                "f4",
                "f5",
                "C_e1",
                "C_e2",
                "C_e3",
                "Tau",
                "phi21",
                "phi22",
                "phi2V",
            ]
        }
        c["c0"] = by_region("c0", np.where(is_slab, cp.c0slab, cp.c0))
        # PRVI uses the global anelastic attenuation
        c["a0"] = by_region("a0", by_type("Global_a0"))
        c["s2"] = by_region("s2", cp.s2)
        # s1 = s2 for these regions
        c["s1"] = c["s2"]

        mag = v["mag"]
        (m_c,) = self._constants(v, {k: ((a,), (b,)) for k, (a, b) in self.M_C.items()})

        # Near-source saturation, equation 4
        h = np.where(
            is_slab,
            np.where(
                mag <= m_c,
                10.0 ** ((1.050 / (m_c - 4.0)) * (mag - m_c) + 1.544),
                35.0,
            ),
            10.0 ** (-0.82 + 0.252 * mag),
        )
        # Path term, equations 2 and 3
        r_ref = np.sqrt(1 + h * h)
        r = np.sqrt(v["dist_rup"] * v["dist_rup"] + h * h)
        f_p = c["c1"] * np.log(r) + self.B_4 * mag * np.log(r / r_ref) + c["a0"] * r

        # Magnitude term, equation 5
        d_m = mag - m_c
        f_m = np.where(mag <= m_c, c["c4"] * d_m + c["c5"] * d_m * d_m, c["c6"] * d_m)

        # Source depth term, equation 6, only for intraslab events
        depth_hyp = v["depth_tor"] + self.THETA_W * self.Z_SLAB
        d, m = cs.d, cs.m
        f_d = np.where(
            is_slab,
            np.where(
                depth_hyp < self.D_B1,
                m * (self.D_B1 - self.D_B2) + d,
                np.where(depth_hyp > self.D_B2, d, m * (depth_hyp - self.D_B2) + d),
            ),
            0.0,
        )
        ln_mean = c["c0"] + f_p + f_m + f_d

        # USGS Alaska adjustment of global interface events
        is_ak_adjusted = (
            self._ak_adjusted
            & ~is_slab
            & ~v["is_alaska"]
            & ~v["is_cascadia"]
            & ~v["is_prvi"]
        )
        adj_ak = np.where(is_ak_adjusted, self.ADJ_AK[rows], 0.0)
        return ln_mean, adj_ak, c

    def _calc_basin(self, v, c, periods) -> np.ndarray:
        """Basin term, equations 11 to 13, with the USGS adjustments."""
        v_s30 = v["v_s30"]
        depth_2_5 = v["depth_2_5"]
        z_2_5 = depth_2_5 * 1000.0
        with np.errstate(divide="ignore", invalid="ignore"):
            ln_z_2_5 = np.log(z_2_5)
        ln_mu_z_2_5 = np.log(10) * (
            self.THETA_1
            * (
                1
                + _erf(
                    (np.log10(v_s30) - np.log10(self.NU_MU))
                    / self.NU_SIGMA
                    / np.sqrt(2)
                )
            )
            + self.THETA_0
        )

        def calc_f_b(delta):
            return np.where(
                delta <= c["C_e1"] / c["C_e3"],
                c["C_e1"],
                np.where(delta >= c["C_e2"] / c["C_e3"], c["C_e2"], c["C_e3"] * delta),
            )

        f_b = calc_f_b(ln_z_2_5 - ln_mu_z_2_5)

        if self._basin:
            is_cascadia = v["is_cascadia"]
            if self._m9:
                # Cascadia interface events, Z_2.5 > 6 km, and T > 1.9 s
                apply_m9 = (
                    is_cascadia & ~v["is_slab"] & (depth_2_5 > 6.0) & (periods > 1.9)
                )
                adj_m9 = np.log(2.0) - calc_f_b(ln_z_2_5 - np.log(self.Z_2_5_REF_M9))
                f_b = np.where(apply_m9, f_b + adj_m9, f_b)

            # USGS scaling for periods longer than 0.5 s and Z_2.5 > 1 km
            # (GmmUtils.deltaZ25scale)
            scale = (np.clip(depth_2_5, 1.0, 3.0) - 1.0) / (3.0 - 1.0)
            scale = np.where(np.isclose(periods, 0.75), scale * 0.585, scale)
            scale = np.where((periods > 0.5) & (depth_2_5 > 1.0), scale, 0.0)
            f_b = np.where(is_cascadia, f_b * scale, f_b)

        return np.where(np.isnan(depth_2_5), 0.0, f_b)

    def _calc_ln_std(self, v, c) -> np.ndarray:
        """Standard deviation, equations 18 to 20."""
        dist_rup = v["dist_rup"]
        v_s30 = v["v_s30"]
        ln_r_ratio = np.log(self.PHI_R2 / self.PHI_R1)
        r_prime = np.maximum(self.PHI_R1, np.minimum(self.PHI_R2, dist_rup))
        var_r = np.log(self.PHI_R2 / r_prime) / ln_r_ratio
        d_var = np.where(
            v_s30 <= self.PHI_V1,
            c["phi2V"] * var_r,
            np.where(
                v_s30 >= self.PHI_V2,
                0.0,
                c["phi2V"]
                * (np.log(self.PHI_V2 / v_s30) / np.log(self.PHI_V2 / self.PHI_V1))
                * var_r,
            ),
        )
        phi_2 = np.where(
            dist_rup <= self.PHI_R1,
            c["phi21"],
            np.where(
                dist_rup >= self.PHI_R2,
                c["phi22"],
                (c["phi22"] - c["phi21"]) / ln_r_ratio * np.log(dist_rup / self.PHI_R1)
                + c["phi21"],
            ),
        )
        phi_total = np.sqrt(d_var + phi_2)
        return np.sqrt(c["Tau"] * c["Tau"] + phi_total * phi_total)

    def _calc_epistemic(self, v, periods) -> np.ndarray:
        """Epistemic standard deviation (sigma_eps), equation 27."""
        sigma_1, sigma_2, t_1, t_2 = self._constants(v, self.EPI)
        with np.errstate(divide="ignore", invalid="ignore"):
            interp = sigma_1 - (sigma_1 - sigma_2) * np.log(periods / t_1) / np.log(
                t_2 / t_1
            )
        # PGA and PGV (periods 0 and -1) use sigma_eps_1
        return np.where(
            (periods <= 0) | (periods < t_1),
            sigma_1,
            np.where(periods < t_2, interp, sigma_2),
        )

    def _calc(self) -> tuple:
        """Calculate the natural log of the response and the standard deviation.

        Returns
        -------
        ln_resp : :class:`np.ndarray`
            natural log of the response
        ln_std : :class:`np.ndarray`
            logarithmic standard deviation

        """
        v = self._scenario_values()
        rows = self._coeff_rows(np.arange(len(self.PERIODS)))
        periods = self.PERIODS[rows]

        # Reference PGA (V_s30 = 760 m/sec) for the nonlinear site term
        ln_pga_ref, adj_ak_pga, _ = self._calc_ln_mean_ref(v, [self.INDEX_PGA])
        pga_ref = np.exp(ln_pga_ref + adj_ak_pga)

        ln_mean, adj_ak, c = self._calc_ln_mean_ref(v, rows)
        v_s30 = v["v_s30"]

        # Linear site term, equation 8
        f_lin = np.where(
            v_s30 <= self.V_1,
            c["s1"] * np.log(v_s30 / self.V_1)
            + c["s2"] * np.log(self.V_1 / self.V_REF),
            np.where(
                v_s30 > c["V2"],
                c["s2"] * np.log(c["V2"] / self.V_REF),
                c["s2"] * np.log(v_s30 / self.V_REF),
            ),
        )
        # Nonlinear site term, equations 9 and 10 (f_1 = 0)
        f_2 = c["f4"] * (
            np.exp(c["f5"] * (np.minimum(v_s30, self.V_REF) - 200.0))
            - np.exp(c["f5"] * (self.V_REF - 200.0))
        )
        f_nl = f_2 * np.log((pga_ref + self.F_3) / self.F_3)
        f_b = self._calc_basin(v, c, periods)

        ln_resp = ln_mean + (f_lin + f_nl + f_b) + adj_ak

        if self._epistemic:
            # Collapse the branches: ln(sum_i w_i exp(mu + delta_i))
            #   = mu + ln(sum_i w_i exp(delta_i))
            delta = self.EPI_Z_SCORE * self._calc_epistemic(v, periods)
            self._epi_delta = delta
            w_lo, w_mid, w_hi = self.EPI_WEIGHTS
            ln_resp = ln_resp + np.log(
                w_lo * np.exp(-delta) + w_mid + w_hi * np.exp(delta)
            )

        ln_std = self._calc_ln_std(v, c)
        # The model is not defined for other event types
        valid = v["is_slab"] | v["is_interface"]
        ln_resp = np.where(valid, ln_resp, np.nan)
        ln_std = np.where(valid, ln_std, np.nan)
        # The standard deviation does not depend on all of the scenario values
        shape = np.broadcast_shapes(np.shape(ln_resp), np.shape(ln_std))
        if np.shape(ln_std) != shape:
            ln_std = np.array(np.broadcast_to(ln_std, shape))
        if np.shape(ln_resp) != shape:
            ln_resp = np.array(np.broadcast_to(ln_resp, shape))
        return ln_resp, ln_std
