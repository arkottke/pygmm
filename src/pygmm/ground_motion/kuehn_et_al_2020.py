"""Kuehn et al. (2020) ground motion model."""

import csv
import gzip
import os

import numpy as np

from .. import model
from ..registry import register


def _load_epistemic():
    """Load the epistemic uncertainty tables used by nshmp-lib.

    Returns
    -------
    mags : :class:`np.ndarray`
        magnitudes of the tables
    tables : :class:`np.ndarray`
        epistemic standard deviation of the median (natural log units) at a
        rupture distance of 10 km, with shape (regions, event types,
        magnitudes, periods) for the regions "global", "alaska", and
        "cascadia", and the event types "interface" and "intraslab". The
        periods are in the order of the coefficient file (PGV, PGA, and the
        spectral periods).

    """
    fname = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "kuehn_et_al_2020-epistemic.csv.gz",
    )
    with gzip.open(fname, "rt") as fp:
        rows = [row for row in csv.reader(fp) if not row[0].startswith("#")]
    header, rows = rows[0], rows[1:]
    periods = np.array(header[4:], dtype=float)
    np.testing.assert_array_equal(periods, KuehnEtAl2020.COEFF.period)

    # As in nshmp-lib, only the values at the shortest distance (10 km) are
    # used. See the class documentation.
    rows = [row for row in rows if float(row[3]) == 10.0]
    mags = np.unique(np.array([row[2] for row in rows], dtype=float))
    tables = np.full((3, 2, len(mags), len(periods)), np.nan)
    for row in rows:
        i = KuehnEtAl2020.EPI_REGIONS.index(row[0])
        j = ["interface", "intraslab"].index(row[1])
        k = np.flatnonzero(mags == float(row[2]))[0]
        tables[i, j, k] = np.array(row[4:], dtype=float)
    assert not np.any(np.isnan(tables))
    return mags, tables


@register(tectonic="subduction")
class KuehnEtAl2020(model.GroundMotionModel):
    r"""Kuehn et al. (2020) ground motion model.

    Partially non-ergodic NGA-Subduction ground motion model for subduction
    interface and intraslab events (KBCG20), with the global model and the
    regional models for Alaska, Cascadia, and Puerto Rico and the U.S. Virgin
    Islands (PRVI).

    Reference: Kuehn, N., Bozorgnia, Y., Campbell, K. W., and Gregor, N.
    (2020). Partially non-ergodic ground-motion model for subduction regions
    using the NGA-Subduction database. PEER Report No. 2020/04, Pacific
    Earthquake Engineering Research Center, Berkeley, CA.

    The model follows the ``KuehnEtAl_2020`` implementation in the USGS
    nshmp-lib code (commit 44728a7d), including its coefficients
    (``KBCG20.csv``, which agrees with ``COEFFS_SEP21/coefficients_KBCG20.csv``
    of https://github.com/nikuehn/KBCG20, and ``KBCG20_PRVI.csv``), its
    epistemic uncertainty tables, and its USGS additions:

    * The magnitude break, :math:`M_b`, is 7.9 (interface) and 7.6
      (intraslab) for the global and PRVI models, 8.0 and 7.2 for Cascadia, and
      8.6 and 7.2 for Alaska. The reference depths of the depth term are 15 km
      (interface) and 50 km (intraslab), and :math:`\theta_{10} = 0`.
    * The regional models use the regional constant (:math:`\theta_1`), the
      regional :math:`V_{s30}` scaling (:math:`\theta_7`), and the regional
      anelastic attenuation coefficient ``theta_6_2_reg`` (forearc and backarc
      sites are not distinguished). The PRVI model uses :math:`\theta_1` and
      :math:`\theta_7` from ``KBCG20_PRVI.csv`` with the global
      :math:`\theta_6`.
    * Spectral accelerations at periods of 0.1 s or less are not less than PGA.
    * If the basin depth (``depth_2_5``) is given, the Cascadia basin term,
      :math:`\min(\theta_{11} + \theta_{12} \Delta \ln Z_{2.5},
      \theta_{11,S})`, is used for all regions (the term is limited by the
      Seattle basin mean residual, :math:`\theta_{11,S}`). Without
      ``depth_2_5`` (*None* or *NaN*), there is no basin term.
    * The standard deviation is :math:`\sqrt{\phi^2 + \tau^2}` with the
      Campbell and Bozorgnia (2014) adjustment for nonlinear site
      amplification, which uses the coefficient :math:`\rho` and
      :math:`\phi_{\ln AF} = 0.3`.

    The options select the nshmp-lib variants (``Gmm`` ids). Each option only
    applies to the regions and event types for which nshmp-lib has the
    variant, and has no effect on the other scenarios:

    * ``basin=True``: the USGS basin depth and period scaling of the basin term
      (zero for periods of 0.5 s or less and for :math:`Z_{2.5} \le` 1 km,
      linear to one at 3 km, times 0.585 at 0.75 s). Cascadia only.
    * ``seattle_basin=True``: the basin term is the Seattle basin mean
      residual, :math:`\theta_{11,S}`, for sites with a basin depth. Cascadia
      only.
    * ``m9=True``: with ``seattle_basin=True``, the basin term is
      :math:`\ln 2` for periods of 2 s and longer and :math:`Z_{2.5}` > 6 km
      (USGS adjustment based on M9 simulations). Cascadia interface events
      only.
    * ``ak_adjusted=True``: the USGS Alaska data adjustment of the global
      interface model is added to the mean (including the reference PGA). Global
      interface events only.
    * ``epistemic=False``: only the central branch of the epistemic
      uncertainty (the nshmp-lib ``*_NO_EPI`` variants).

    The nshmp-lib ``Gmm`` ids correspond to the options as follows (``region``
    and ``event_type`` are scenario values; ``KBCG_20_<REGION>_<TYPE>`` is
    ``region`` "global", "cascadia", "alaska", or "prvi" and ``event_type``
    "interface" or "intraslab" with the default options):

    ================================================ ==========================
    ``Gmm``                                          Options
    ================================================ ==========================
    ``KBCG_20_<REGION>_<TYPE>``                      defaults
    ``KBCG_20_GLOBAL_<TYPE>_NO_EPI``                 ``epistemic=False``
    ``KBCG_20_GLOBAL_INTERFACE_AK_ADJUSTED``         ``ak_adjusted=True``
    ``KBCG_20_CASCADIA_<TYPE>_BASIN``                ``basin=True``
    ``KBCG_20_CASCADIA_INTERFACE_SEATTLE_BASIN``     ``basin=True``,
                                                     ``seattle_basin=True``
    ``KBCG_20_CASCADIA_INTERFACE_SEATTLE_BASIN_M9``  ``basin=True``,
                                                     ``seattle_basin=True``,
                                                     ``m9=True``
    ``KBCG_20_CASCADIA_SLAB_SEATTLE_BASIN``          ``seattle_basin=True``
    ================================================ ==========================

    The epistemic uncertainty of the median is represented as in nshmp-lib by
    three branches with means :math:`\mu - 1.645 \psi`, :math:`\mu`, and
    :math:`\mu + 1.645 \psi` and weights of 0.185, 0.63, and 0.185, where
    :math:`\psi` is the standard deviation of the median from the tables of
    Kuehn et al. (global, Cascadia, or Alaska; PRVI uses the global tables).
    The branches are collapsed to a single median, :math:`\ln \sum_i w_i
    \exp(\mu_i)`, and the standard deviation, which is the same for all
    branches. The collapsed values are not suitable for hazard calculations,
    which should sum the weighted exceedance probabilities of the branches,
    as nshmp-lib does; the branches are given by
    :meth:`~pygmm.model.GroundMotionModel.ln_branches`. nshmp-lib interpolates
    the tables linearly in magnitude (with
    the magnitude limited to 4 to 9.5) and is meant to interpolate in the
    logarithm of the distance, but its distance lookup compares the base-10
    logarithm of the distance with the distances of the tables (10 to 1000 km)
    and therefore always uses the values at 10 km. That is reproduced here.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``depth_tor``, ``v_s30``, ``depth_2_5``, ``event_type``, and ``region``)
    can be a scalar or an array, and the arrays are broadcast against each
    other, so events of different types and regions can be combined. For a
    scalar scenario, the response and standard deviation have one value per
    period, as in other models. For arrays of N scenarios, they have shape
    (N, periods), so, for example, ``pga`` has shape (N,). ``depth_tor`` is
    required. ``depth_2_5`` can be *None* or contain *NaN* for sites without a
    basin term. ``region`` defaults to "global". Results for events that are
    neither interface nor intraslab events are *NaN*.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    basin : bool, optional
        if *True*, apply the USGS basin depth and period scaling to the basin
        term of Cascadia events. Default is *False*.
    seattle_basin : bool, optional
        if *True*, use the Seattle basin term for Cascadia events. Default is
        *False*.
    m9 : bool, optional
        if *True*, use the USGS M9 adjustment of the Seattle basin term for
        Cascadia interface events. Requires ``seattle_basin=True``. Default
        is *False*.
    ak_adjusted : bool, optional
        if *True*, add the USGS Alaska adjustment to global interface events.
        Default is *False*.
    epistemic : bool, optional
        if *True* (default), include the epistemic uncertainty branches of
        the median. If *False*, only the central branch is used.
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods), and/or "psa_all" (all 21 periods). Computing only the
        needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> from pygmm.model import Scenario
    >>> from pygmm.ground_motion.kuehn_et_al_2020 import KuehnEtAl2020
    >>> s = Scenario(
    ...     mag=np.array([7.0, 9.0]), dist_rup=np.array([80.0, 100.0]),
    ...     depth_tor=np.array([50.0, 15.0]), v_s30=400.0,
    ...     depth_2_5=np.array([np.nan, 4.0]),
    ...     event_type=np.array(["intraslab", "interface"]),
    ...     region="cascadia")
    >>> m = KuehnEtAl2020(s, basin=True)
    >>> m.pga.shape, m.spec_accels.shape
    ((2,), (2, 21))
    >>> KuehnEtAl2020(s, epistemic=False, ims=["pga"]).pga.shape
    (2,)
    >>> m = KuehnEtAl2020(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

    """

    NAME = "Kuehn et al. (2020)"
    ABBREV = "KBCG20"

    # Load the coefficients for the model
    #: USGS Alaska interface adjustment, and the epistemic uncertainty tables
    #: (set below the class)
    _AK_ADJ: np.ndarray
    _EPI_MAGS: np.ndarray
    _EPI_TABLES: np.ndarray

    COEFF = model.load_data_file("kuehn_et_al_2020.csv", 1)
    COEFF_PRVI = model.load_data_file("kuehn_et_al_2020-prvi.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGV = 0
    INDEX_PGA = 1
    INDICES_PSA = np.arange(2, 23)

    # PGV is in cm/sec
    PGV_SCALE = 1.0

    REGIONS = ["global", "alaska", "cascadia", "prvi"]
    # Regions of the epistemic uncertainty tables (PRVI uses the global tables)
    EPI_REGIONS = ["global", "alaska", "cascadia"]

    # Magnitude breaks (interface, intraslab) by region
    MAG_BREAKS = {
        "global": (7.9, 7.6),
        "prvi": (7.9, 7.6),
        "cascadia": (8.0, 7.2),
        "alaska": (8.6, 7.2),
    }

    # Constants
    DELTA_M = 0.1
    MAG_REF = 6.0
    DELTA_Z = 1.0
    # Depth break offsets and reference depths (interface, intraslab)
    DEPTH_BREAKS = (30.0, 80.0)
    DEPTH_REFS = (15.0, 50.0)
    # Nonlinear site model of Campbell and Bozorgnia (2014)
    C = 1.88
    N = 1.18
    V_S30_ROCK = 1100.0
    PHI_LN_AF_SQ = 0.09
    # Basin reference depth model
    THETA_Z = (8.29404964010203, 2.30258509299405, 6.39692965521615, 0.27081459)
    # USGS basin scaling (km)
    BASIN_Z2P5_UPPER = 1.0
    BASIN_Z2P5_LOWER = 3.0
    # Epistemic branches
    EPI_Z_SCORE = 1.645
    EPI_WEIGHTS = (0.185, 0.63, 0.185)

    LIMITS = dict(
        mag=(5.0, 9.5),
        dist_rup=(10.0, 1000.0),
        depth_tor=(0.0, 200.0),
        v_s30=(150.0, 1500.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 5.0, 9.5),
        model.NumericParameter("dist_rup", True, 10.0, 1000.0),
        # Required, but checked in the model for a clearer error message
        model.NumericParameter("depth_tor", False, 0.0, 200.0),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.NumericParameter("depth_2_5", False),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
        model.CategoricalParameter("region", False, REGIONS, "global"),
    ]

    def __init__(
        self,
        scenario: model.Scenario,
        basin: bool = False,
        seattle_basin: bool = False,
        m9: bool = False,
        ak_adjusted: bool = False,
        epistemic: bool = True,
        ims=None,
    ):
        """Initialize the model."""
        super().__init__(scenario, ims)
        self._basin = bool(basin)
        self._seattle_basin = bool(seattle_basin)
        self._m9 = bool(m9)
        self._ak_adjusted = bool(ak_adjusted)
        self._epistemic = bool(epistemic)
        self._epi_delta: np.ndarray | None = None
        if self._m9 and not self._seattle_basin:
            raise ValueError("m9=True requires seattle_basin=True")
        self._ln_resp, self._ln_std = self._calc()
        if self._epistemic and self._epi_delta is not None:
            self._branches = model.symmetric_branches(
                self._ln_resp, self._ln_std, self._epi_delta, self.EPI_WEIGHTS
            )

    @property
    def basin(self) -> bool:
        """If the USGS basin scaling is applied to Cascadia events."""
        return self._basin

    @property
    def seattle_basin(self) -> bool:
        """If the Seattle basin term is used for Cascadia events."""
        return self._seattle_basin

    @property
    def m9(self) -> bool:
        """If the M9 adjustment is used for Cascadia interface events."""
        return self._m9

    @property
    def ak_adjusted(self) -> bool:
        """If the Alaska adjustment is added to global interface events."""
        return self._ak_adjusted

    @property
    def epistemic(self) -> bool:
        """If the epistemic uncertainty branches of the median are included."""
        return self._epistemic

    def _calc(self):
        """Calculate the natural logarithm of the response and its standard deviation.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        s = self._scenario
        if s.depth_tor is None:
            raise ValueError("depth_tor is required")

        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        event_type = model.as_column(s.event_type)
        region = model.as_column(s.region)
        is_slab = model.equals(event_type, "intraslab")
        is_interface = model.equals(event_type, "interface")
        regions = {r: model.equals(region, r) for r in self.REGIONS}
        is_cascadia = regions["cascadia"]

        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)
        depth_tor = model.as_column(s.depth_tor)
        v_s30 = model.as_column(s.v_s30)
        depth_2_5 = model.as_column(np.nan if s.depth_2_5 is None else s.depth_2_5)

        # Magnitude break
        mag_break = np.select(
            [regions[r] for r in ["cascadia", "alaska"]],
            [
                np.where(is_slab, *self.MAG_BREAKS[r][::-1])
                for r in ["cascadia", "alaska"]
            ],
            np.where(is_slab, *self.MAG_BREAKS["global"][::-1]),
        )

        # Options
        use_basin = self._basin & is_cascadia
        use_seattle = self._seattle_basin & is_cascadia
        use_m9 = self._m9 & is_cascadia & is_interface
        use_ak = self._ak_adjusted & regions["global"] & is_interface

        args = (
            is_slab,
            regions,
            mag_break,
            use_basin,
            use_seattle,
            use_m9,
            use_ak,
            mag,
            dist_rup,
            depth_tor,
            depth_2_5,
        )

        # Reference PGA (V_s30 = 1100 m/sec)
        idx_pga = [self.INDEX_PGA]
        # The reference PGA does not include the basin term
        pga_ref = np.exp(
            self._calc_ln_mean(idx_pga, *args, self.V_S30_ROCK, 0.0, with_basin=False)
        )

        indices = (
            np.arange(len(self.PERIODS)) if self._indices is None else self._indices
        )
        ln_resp = self._calc_ln_mean(indices, *args, v_s30, pga_ref)

        # Spectral accelerations at short periods are not less than PGA
        period = self.PERIODS[indices]
        is_short = (0 < period) & (period <= 0.1)
        if np.any(is_short):
            ln_pga = self._calc_ln_mean(idx_pga, *args, v_s30, pga_ref)
            ln_resp = np.where(is_short, np.maximum(ln_resp, ln_pga), ln_resp)

        if self._epistemic:
            psi = self._calc_epistemic(indices, is_slab, regions, mag)
            # Collapse the branches: ln(sum_i w_i exp(mu + delta_i))
            #   = mu + ln(sum_i w_i exp(delta_i))
            w_lo, w_mid, w_hi = self.EPI_WEIGHTS
            delta = self.EPI_Z_SCORE * psi
            self._epi_delta = delta
            ln_resp = ln_resp + np.log(
                w_lo * np.exp(-delta) + w_mid + w_hi * np.exp(delta)
            )

        ln_std = self._calc_ln_std(indices, v_s30, pga_ref)

        # The model is not defined for other event types
        valid = is_slab | is_interface
        ln_resp = np.where(valid, ln_resp, np.nan)
        # The standard deviation does not depend on all of the scenario values
        ln_std = np.where(valid, np.broadcast_to(ln_std, ln_resp.shape), np.nan)
        return ln_resp, ln_std

    def _calc_ln_mean(
        self,
        indices,
        is_slab,
        regions,
        mag_break,
        use_basin,
        use_seattle,
        use_m9,
        use_ak,
        mag,
        dist_rup,
        depth_tor,
        depth_2_5,
        v_s30,
        pga_ref,
        with_basin=True,
    ) -> np.ndarray:
        """Calculate the mean (natural log units) without the epistemic branches.

        Parameters
        ----------
        indices : array_like
            coefficient rows (periods) to compute
        is_slab : array_like
            if the event is an intraslab event
        regions : dict
            if the event is in each region
        mag_break : array_like
            magnitude break
        use_basin, use_seattle, use_m9, use_ak : array_like
            if the USGS basin scaling, the Seattle basin term, the M9
            adjustment, and the Alaska adjustment are used
        mag : array_like
            moment magnitude
        dist_rup : array_like
            closest distance to the rupture [km]
        depth_tor : array_like
            depth to the top of the rupture [km]
        depth_2_5 : array_like
            depth to a shear-wave velocity of 2.5 km/sec [km], *NaN* for no
            basin term
        v_s30 : float or array_like
            time-averaged shear-wave velocity over the top 30 m [m/sec]
        pga_ref : float or array_like
            peak ground acceleration [g] at the reference condition
            (V_s30 = 1100 m/sec)
        with_basin : bool, optional
            if *False*, the basin term is not included (for the reference PGA)

        Returns
        -------
        ln_mean : class:`np.array`:
            natural log of the response

        """
        c = self.COEFF[indices]
        c_prvi = self.COEFF_PRVI[indices]
        period = c.period

        # Regional coefficients
        cascadia, alaska, prvi = (regions[r] for r in ["cascadia", "alaska", "prvi"])
        theta_1 = np.select(
            [cascadia, alaska, prvi],
            [
                np.where(is_slab, c.theta_1_slab_reg_Ca, c.theta_1_if_reg_Ca),
                np.where(is_slab, c.theta_1_slab_reg_Al, c.theta_1_if_reg_Al),
                np.where(is_slab, c_prvi.mu_c_1_slab, c_prvi.mu_c_1_if),
            ],
            np.where(is_slab, c.mu_theta_1_slab, c.mu_theta_1_if),
        )
        theta_6 = np.select(
            [cascadia, alaska], [c.theta_6_2_reg_Ca, c.theta_6_2_reg_Al], c.mu_theta_6
        )
        theta_7 = np.select(
            [cascadia, alaska, prvi],
            [c.theta_7_reg_Ca, c.theta_7_reg_Al, c_prvi.mu_c_7],
            c.mu_theta_7,
        )

        # Magnitude term
        theta_4 = np.where(is_slab, c.theta_4_slab, c.theta_4_if)
        f_mag = self._log_hinge(
            mag,
            mag_break,
            theta_4 * (mag_break - self.MAG_REF),
            theta_4,
            c.theta_5,
            self.DELTA_M,
        )

        # Geometrical spreading with the finite-fault term
        theta_2 = np.where(is_slab, c.theta_2_slab, c.theta_2_if)
        f_geom = (theta_2 + c.theta_3 * mag) * np.log(
            dist_rup + 10 ** (c.nft_1 + c.nft_2 * (mag - 6.0))
        )

        # Depth term
        depth_break = np.where(
            is_slab, self.DEPTH_BREAKS[1] + c.dzb_slab, self.DEPTH_BREAKS[0] + c.dzb_if
        )
        depth_ref = np.where(is_slab, *self.DEPTH_REFS[::-1])
        theta_9 = np.where(is_slab, c.theta_9_slab, c.theta_9_if)
        f_depth = self._log_hinge(
            depth_tor,
            depth_break,
            theta_9 * (depth_break - depth_ref),
            theta_9,
            0.0,
            self.DELTA_Z,
        )

        # Anelastic attenuation
        f_atten = theta_6 * dist_rup

        # Site term
        vs_ratio = v_s30 / c.k1
        ln_vs = np.log(vs_ratio)
        f_site = np.where(
            v_s30 <= c.k1,
            theta_7 * ln_vs
            + c.k2
            * (np.log(pga_ref + self.C * vs_ratio**self.N) - np.log(pga_ref + self.C)),
            (theta_7 + c.k2 * self.N) * ln_vs,
        )

        ln_mean = theta_1 + f_mag + f_geom + f_depth + f_atten + f_site

        if with_basin:
            ln_mean = ln_mean + self._calc_f_basin(
                c, period, is_slab, use_basin, use_seattle, use_m9, depth_2_5, v_s30
            )

        # USGS Alaska adjustment of global interface events
        if np.any(use_ak):
            ln_mean = ln_mean + np.where(use_ak, self._AK_ADJ[indices], 0.0)

        return ln_mean

    @staticmethod
    def _log_hinge(x, x_0, a, b_0, b_1, delta):
        """Smoothed hinge function (Eq. 4.3)."""
        return (
            a
            + b_0 * (x - x_0)
            + (b_1 - b_0) * delta * np.log(1 + np.exp((x - x_0) / delta))
        )

    def _calc_f_basin(
        self, c, period, is_slab, use_basin, use_seattle, use_m9, depth_2_5, v_s30
    ) -> np.ndarray:
        """Calculate the basin term."""
        theta_z1, theta_z2, theta_z3, theta_z4 = self.THETA_Z
        vs_term = np.exp((np.log(v_s30) - theta_z3) / theta_z4)
        ln_z_ref = theta_z1 + (theta_z2 - theta_z1) * vs_term / (1.0 + vs_term)
        with np.errstate(invalid="ignore", divide="ignore"):
            delta_ln_z = np.log(depth_2_5 * 1000.0) - ln_z_ref

        f_basin = np.where(
            use_seattle,
            np.where(
                use_m9 & ~is_slab & (depth_2_5 > 6.0) & (period > 1.9),
                np.log(2.0),
                c.mean_residual_Seattle_basin,
            ),
            np.minimum(
                c.intercept_Ca + c.slope_Ca * delta_ln_z,
                c.mean_residual_Seattle_basin,
            ),
        )

        if np.any(use_basin):
            # USGS basin depth and period scaling (GmmUtils.deltaZ25scale)
            scale = (
                np.clip(depth_2_5, self.BASIN_Z2P5_UPPER, self.BASIN_Z2P5_LOWER)
                - self.BASIN_Z2P5_UPPER
            ) / (self.BASIN_Z2P5_LOWER - self.BASIN_Z2P5_UPPER)
            scale = np.where(period == 0.75, scale * 0.585, scale)
            scale = np.where(
                (depth_2_5 > self.BASIN_Z2P5_UPPER) & (period > 0.5), scale, 0.0
            )
            f_basin = np.where(use_basin, f_basin * scale, f_basin)

        # No basin term without the basin depth
        return np.where(np.isnan(depth_2_5), 0.0, f_basin)

    def _calc_epistemic(self, indices, is_slab, regions, mag) -> np.ndarray:
        """Standard deviation of the median (natural log units)."""
        mags, tables = self._EPI_MAGS, self._EPI_TABLES[..., indices]
        i_region = np.select([regions["alaska"], regions["cascadia"]], [1, 2], 0)[
            ..., 0
        ]
        i_type = np.where(is_slab, 1, 0)[..., 0]
        m = np.asarray(mag)[..., 0]
        # Clamping linear interpolation in magnitude (GroundMotionTables)
        i_mag = np.clip(np.searchsorted(mags, m, side="right") - 1, 0, len(mags) - 2)
        m_lo = mags[i_mag]
        m_hi = mags[i_mag + 1]
        frac = np.where(
            m < m_lo, 0.0, np.where(m > m_hi, 1.0, (m - m_lo) / (m_hi - m_lo))
        )[..., np.newaxis]
        lo = tables[i_region, i_type, i_mag]
        hi = tables[i_region, i_type, i_mag + 1]
        return lo + frac * (hi - lo)

    def _calc_ln_std(self, indices, v_s30, pga_ref) -> np.ndarray:
        """Calculate the standard deviation with the nonlinear site adjustment."""
        c = self.COEFF[indices]
        c_pga = self.COEFF[self.INDEX_PGA]

        alpha = np.where(
            v_s30 < c.k1,
            c.k2
            * pga_ref
            * (
                1 / (pga_ref + self.C * (v_s30 / c.k1) ** self.N)
                - 1 / (pga_ref + self.C)
            ),
            0.0,
        )
        alpha_tau = alpha * c_pga.tau
        tau_sq = (
            c.tau * c.tau
            + alpha_tau * alpha_tau
            + 2.0 * alpha * c.rho * c.tau * c_pga.tau
        )
        phi_b = np.sqrt(c.phi * c.phi - self.PHI_LN_AF_SQ)
        alpha_phi_pga_b = alpha * np.sqrt(c_pga.phi * c_pga.phi - self.PHI_LN_AF_SQ)
        phi_sq = (
            c.phi * c.phi
            + alpha_phi_pga_b * alpha_phi_pga_b
            + 2.0 * c.rho * phi_b * alpha_phi_pga_b
        )
        return np.sqrt(phi_sq + tau_sq)


def _load_ak_adjustment():
    """Load the USGS Alaska adjustment in the order of the coefficients."""
    fname = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        "data",
        "nga_sub_ak_interface_adjustment.csv",
    )
    with open(fname) as fp:
        rows = [row for row in csv.reader(fp) if not row[0].startswith("#")]
    names = {"PGV": -1.0, "PGA": 0.0}
    values = {names[t] if t in names else float(t): float(v) for t, v in rows[1:]}
    return np.array([values[p] for p in KuehnEtAl2020.PERIODS])


KuehnEtAl2020._AK_ADJ = _load_ak_adjustment()
KuehnEtAl2020._EPI_MAGS, KuehnEtAl2020._EPI_TABLES = _load_epistemic()
