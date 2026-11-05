"""Model for the Campbell and Bozorgnia (2014) ground motion model."""

import logging

import numpy as np

from .. import model
from .._types import ArrayLike
from ..registry import register
from . import nshmp as usgs

__author__ = "Albert Kottke"


def _min(a: ArrayLike, b: ArrayLike) -> ArrayLike:
    """Element-wise ``min(a, b)``.

    Same as the builtin ``min``, which returns `a` unless `b` is less, so a NaN
    value of `b` is ignored, unlike :func:`np.minimum`.
    """
    return np.where(b < a, b, a)[()]


def _pow(x: ArrayLike, y: float) -> ArrayLike:
    """Element-wise ``x ** y``, computed with the C ``pow``.

    For arrays, ``x ** 2`` is computed as ``x * x``, which can differ by one unit
    in the last place from ``pow(x, 2)`` used by ``**`` for scalars on some
    platforms. This gives arrays of scenario values the same results as scalars.
    """
    return np.float_power(x, y)


@register(tectonic="active_crustal")
class CampbellBozorgnia2014(usgs.NshmpOptions, model.GroundMotionModel):
    r"""Campbell and Bozorgnia (2014, :cite:p:`campbell14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``dist_jb``, ``dist_x``, ``dip``, ``v_s30``, ``depth_1_0``, ``depth_2_5``,
    ``depth_tor``, ``depth_bor``, ``depth_bot``, ``depth_hyp``, ``width``,
    ``mechanism``, and ``region``) can be a scalar or an array, and the arrays
    are broadcast against each other. Values that are not provided (e.g.,
    ``depth_tor``, ``width``, and ``depth_hyp``) are estimated for each
    scenario. For a scalar scenario, the response and standard deviation have
    one value per period, as in other models. For arrays of N scenarios, they
    have shape (N, periods), so, for example, ``pga`` has shape (N,) and
    ``spec_accels`` has shape (N, 21).

    The options select the variants of the model in the USGS nshmp-lib library
    (``CampbellBozorgnia_2014``, commit 44728a7d; see
    :mod:`pygmm.ground_motion.nshmp`), and are off by default, which gives the
    published model. With ``nshmp=True``, the model is computed as in
    nshmp-lib, which:

    - computes ``width`` from the depth to the bottom of the rupture,
      :math:`(Z_{bor} - Z_{tor}) / \sin \delta`, if ``width`` is *None* and
      ``depth_bor`` is given, and ``depth_hyp`` as the middle of the rupture,
      :math:`Z_{tor} + W \sin(\delta) / 2`, if ``depth_hyp`` is *None*
      (``Faults.hypocentralDepth``),
    - uses the :math:`V_{S30}` based reference :math:`Z_{2.5}` (equation 33)
      if ``depth_2_5`` is *None* or NaN,
    - computes the reference rock PGA (:math:`V_{S30}` = 1100 m/sec) with
      :math:`Z_{2.5}` = 0.398 km (or, with ``basin=True``, the reference depth
      for 1100 m/sec, as the published model), and only uses it if
      :math:`V_{S30} < k_1` of the period (otherwise it is 0, also for the
      PGA limit of the short periods),
    - for periods from 0.01 to 0.25 s (the published model: shorter than 0.25
      s), uses the PGA if the spectral acceleration is smaller,
    - uses the :math:`f_2(R_x)` hanging-wall taper for :math:`R_x \ge R_1`,
    - and treats an unspecified style of faulting (``mechanism`` of "U") as
      strike-slip.

    nshmp-lib only has the California model, so ``region`` should not be given
    with ``nshmp=True``. The nshmp-lib ``Gmm`` ids are (see :attr:`GMM_IDS`):

    ========================  =========================================
    nshmp-lib ``Gmm``         Options
    ========================  =========================================
    ``CB_14``                 ``nshmp=True, epistemic=True``
    ``CB_14_BASE``            ``nshmp=True``
    ``CB_14_BASIN``           ``nshmp=True, epistemic=True, basin=True``
    ``CB_14_CYBERSHAKE``      ``nshmp=True, epistemic=True, basin=True,
                              cybershake=True``
    ``CB_14_PRVI``            ``nshmp=True, prvi=True``
    ========================  =========================================

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods), and/or "psa_all" (all 21 periods). Computing only
        the needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.
    nshmp : bool, optional
        compute the model as in nshmp-lib (default *False*)
    epistemic : bool, optional
        include the USGS epistemic branches (default *False*)
    basin : bool, optional
        use the USGS deep basin model (default *False*)
    cybershake : bool, optional
        use the CyberShake basin terms (default *False*, requires
        ``basin=True``)
    prvi : bool, optional
        use the PRVI coefficients (default *False*)

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 50.0]),
    ...     dist_jb=np.array([8.0, 49.0]), dist_x=np.array([-5.0, 20.0]),
    ...     dip=np.array([90.0, 45.0]), v_s30=400.0,
    ...     mechanism=np.array(["SS", "RS"]))
    >>> pygmm.CampbellBozorgnia2014(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.CampbellBozorgnia2014(s, ims=["pga", "psa_ngawest2_21"])
    >>> m.spec_accels.shape
    (2, 21)

    The nshmp-lib ``CB_14_BASIN`` variant, with the width and hypocentral
    depth from the depth to the bottom of the rupture:

    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 50.0]),
    ...     dist_jb=np.array([8.0, 48.0]), dist_x=np.array([8.0, -48.0]),
    ...     dip=60.0, depth_tor=2.0, depth_bor=15.0, mechanism="RS",
    ...     v_s30=400.0, depth_2_5=np.array([np.nan, 4.0]))
    >>> options = pygmm.CampbellBozorgnia2014.GMM_IDS["CB_14_BASIN"]
    >>> pygmm.CampbellBozorgnia2014(s, **options).spec_accels.shape
    (2, 21)
    """

    NAME = "Campbell & Bozorgnia (2014)"
    ABBREV = "CB14"

    # Reference velocity (m/sec)
    V_REF = 1100.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("campbell_bozorgnia_2014.csv", 2)

    PERIODS = COEFF["period"]

    #: nshmp-lib PRVI coefficients
    COEFF_PRVI = usgs.load_coeffs("campbell_bozorgnia_2014-prvi.csv", PERIODS)
    #: nshmp-lib CyberShake basin coefficients (``slope`` replaces
    #: :math:`c_{16} k_3`)
    COEFF_CYBERSHAKE = usgs.load_coeffs(
        "campbell_bozorgnia_2014-cybershake.csv", PERIODS
    )

    #: Model options of the nshmp-lib ``Gmm`` ids
    GMM_IDS = {
        "CB_14": dict(nshmp=True, epistemic=True),
        "CB_14_BASE": dict(nshmp=True),
        "CB_14_BASIN": dict(nshmp=True, epistemic=True, basin=True),
        "CB_14_CYBERSHAKE": dict(
            nshmp=True, epistemic=True, basin=True, cybershake=True
        ),
        "CB_14_PRVI": dict(nshmp=True, prvi=True),
    }

    #: :math:`Z_{2.5}` (km) of the reference rock PGA in nshmp-lib
    DEPTH_2_5_REF_NSHMP = 0.398

    # Period independent model coefficients
    COEFF_C = 1.88
    COEFF_N = 1.18
    COEEF_H_4 = 1

    INDICES_PSA = np.arange(21)
    INDEX_PGA = -2
    INDEX_PGV = -1

    PARAMS = [
        model.NumericParameter("depth_1_0", False),
        model.NumericParameter("depth_2_5", False, 0, 10),
        model.NumericParameter("depth_bor", False),
        model.NumericParameter("depth_bot", False, default=15.0),
        model.NumericParameter("depth_hyp", False, 0, 20),
        model.NumericParameter("depth_tor", False, 0, 20),
        model.NumericParameter("dip", True, 15, 90),
        model.NumericParameter("dist_jb", True),
        model.NumericParameter("dist_rup", True, None, 300),
        model.NumericParameter("dist_x", True),
        model.NumericParameter("mag", True, 3.3, 8.5),
        model.NumericParameter("v_s30", True, 150, 1500),
        model.NumericParameter("width", False),
        model.CategoricalParameter(
            "region",
            False,
            ["global", "california", "japan", "italy", "china"],
            "global",
        ),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def _check_inputs(self) -> None:
        """Check the inputs."""
        s = self._scenario
        if self._nshmp:
            # nshmp-lib defaults
            if s["mechanism"] is not None:
                # No style of faulting term for an unspecified mechanism
                mechanism = np.asarray(s["mechanism"])
                s["mechanism"] = np.where(
                    model.equals(mechanism, "U"), "SS", mechanism
                )[()]
            if s["depth_2_5"] is None:
                # Reference depth from V_s30
                s["depth_2_5"] = np.nan
        super()._check_inputs()

        if np.ndim(s.mechanism) > 0 or np.ndim(s.mag) > 0:
            self._check_mag_by_mechanism()
        else:
            for mech, limit in [("SS", 8.5), ("RS", 8.0), ("NS", 7.5)]:
                if mech == s.mechanism and s.mag > limit:
                    logging.warning(
                        "Magnitude of %g is greater than the recommended limit of"
                        "%g for %s style faults",
                        s.mag,
                        limit,
                        mech,
                    )

        # Estimate the values that are not provided. The estimates accept
        # arrays, so each scenario of a vectorized scenario gets its own value.
        if s.depth_2_5 is None:
            s.depth_2_5 = self.calc_depth_2_5(s.v_s30, s.region, s.depth_1_0)

        if s.depth_tor is None:
            s.depth_tor = self._calc_depth_tor(s.mag, s.mechanism)

        if s.width is None:
            if self._nshmp and s.depth_bor is not None:
                s.width = (np.asarray(s.depth_bor) - np.asarray(s.depth_tor)) / np.sin(
                    np.radians(s.dip)
                )
            else:
                s.width = CampbellBozorgnia2014.calc_width(
                    s.mag, s.dip, s.depth_tor, s.depth_bot
                )

        if s.depth_bor is None:
            s.depth_bor = self.calc_depth_bor(s.depth_tor, s.dip, s.width)

        if s.depth_hyp is None:
            if self._nshmp:
                # Middle of the rupture (nshmp-lib Faults.hypocentralDepth)
                s.depth_hyp = (
                    np.asarray(s.depth_tor)
                    + np.sin(np.radians(s.dip)) * np.asarray(s.width) / 2.0
                )
            else:
                s.depth_hyp = CampbellBozorgnia2014.calc_depth_hyp(
                    s.mag, s.dip, s.depth_tor, s.depth_bor
                )

    def _check_mag_by_mechanism(self) -> None:
        """Check mechanism specific magnitude limits for vectorized scenarios."""
        s = self._scenario
        mag, mechanism = np.broadcast_arrays(np.asarray(s.mag), np.asarray(s.mechanism))
        for mech, limit in [("SS", 8.5), ("RS", 8.0), ("NS", 7.5)]:
            is_mech = model.equals(mechanism, mech)
            above = is_mech & (mag > limit)
            if np.any(above):
                logging.warning(
                    "Magnitude (%d of %d earthquakes, maximum of %g) is greater "
                    "than the recommended limit of %g for %s style faults",
                    np.count_nonzero(above),
                    np.count_nonzero(is_mech),
                    np.max(mag[above]),
                    limit,
                    mech,
                )

    def __init__(
        self,
        scenario: model.Scenario,
        ims=None,
        nshmp: bool = False,
        epistemic: bool = False,
        basin: bool = False,
        cybershake: bool = False,
        prvi: bool = False,
    ):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed. nshmp, epistemic, basin, cybershake, prvi (bool, optional):
            the nshmp-lib options (see :mod:`pygmm.ground_motion.nshmp`).
        """
        self._set_options(nshmp, epistemic, basin, cybershake, prvi)
        super().__init__(scenario, ims)
        self._coeff = self._options_coeff(self.COEFF)

        # The nonlinear site term and the standard deviation depend on PGA at
        # the reference condition, which only needs the PGA coefficients
        c_ref = self._coeff[[self.INDEX_PGA]]
        pga_ref = np.exp(self._calc_ln_resp(np.nan, self.V_REF, c_ref)[..., 0])
        pga_rock = None
        if self._nshmp:
            # nshmp-lib only uses the reference PGA for V_s30 < k_1 of the period
            pga_rock = np.where(
                model.as_column(self._scenario.v_s30)
                < self._coeff_rows(self._coeff).k_1,
                model.as_column(pga_ref),
                0.0,
            )
        ln_resp = self._calc_ln_resp(pga_ref, self._scenario.v_s30, pga_rock=pga_rock)

        # PSA at periods shorter than 0.25 s (from 0.01 to 0.25 s in nshmp-lib)
        # is not less than PGA
        periods = self._coeff_rows(self.PERIODS)
        if self._nshmp:
            is_short = (0.01 <= periods) & (periods <= 0.25)
        else:
            is_short = (0 < periods) & (periods < 0.25)
        if np.any(is_short):
            c_pga = self._coeff[[self.INDEX_PGA]]
            ln_pga = self._calc_ln_resp(
                pga_ref, self._scenario.v_s30, c_pga, pga_rock=pga_rock
            )
            ln_resp = np.where(is_short, np.maximum(ln_resp, ln_pga), ln_resp)
        self._ln_resp = ln_resp
        self._ln_std, self._tau, self._phi = self._calc_ln_std(pga_ref)
        if self._epistemic:
            self._ln_resp, self._branches = usgs.epistemic_branches(
                self._ln_resp,
                self._ln_std,
                model.as_column(self._scenario.mag),
                model.as_column(self._scenario.dist_jb),
            )
        if self._nshmp or self._epistemic:
            self._ln_resp, self._ln_std = usgs.broadcast(self._ln_resp, self._ln_std)

    def _calc_ln_resp(
        self, pga_ref: ArrayLike, v_s30: ArrayLike, c=None, pga_rock=None
    ) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        pga_ref : float or array_like
            peak ground acceleration [g] at the reference
            condition. If :class:`np.nan`, then the reference condition is
            computed, which uses the basin depth estimated from `v_s30`.
        v_s30 : float or array_like
            time-averaged shear-wave velocity over the top 30 m
            of the site (:math:`V_{s30}`, m/s).
        c : :class:`numpy.recarray`, optional
            coefficients for the periods to compute. If *None*, the
            coefficients for the requested intensity measures are used.
        pga_rock : array_like, optional
            reference PGA [g] of each computed period (with a trailing axis
            for the periods), which replaces `pga_ref` (nshmp-lib).

        Returns
        -------
            class:`np.array`: Natural log of the response.

        """
        s = self._scenario
        if c is None:
            c = self._coeff_rows(self._coeff)
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        mag = model.as_column(s.mag)
        dip = model.as_column(s.dip)
        dist_rup = model.as_column(s.dist_rup)
        dist_jb = model.as_column(s.dist_jb)
        dist_x = model.as_column(s.dist_x)
        depth_tor = model.as_column(s.depth_tor)
        mechanism = model.as_column(s.mechanism)
        region = model.as_column(s.region)

        # Magnitude term. Each slope applies above its magnitude.
        f_mag = c.c_0 + c.c_1 * mag
        for min_mag, slope in ([4.5, c.c_2], [5.5, c.c_3], [6.5, c.c_4]):
            f_mag = f_mag + np.where(min_mag < mag, slope * (mag - min_mag), 0.0)

        # Geometric attenuation term
        f_dis = (c.c_5 + c.c_6 * mag) * np.log(np.sqrt(_pow(dist_rup, 2) + c.c_7**2))

        # Style of faulting term
        taper = np.clip(mag - 4.5, 0, 1)
        f_flt = np.select(
            [model.equals(mechanism, "RS"), model.equals(mechanism, "NS")],
            [c.c_8 * taper, c.c_9 * taper],
            # Strike-slip
            default=0.0,
        )

        # Hanging-wall term
        R_1 = model.as_column(s.width) * np.cos(np.radians(dip))
        R_2 = 62 * mag - 350
        # Both branches are computed for every scenario, and the ratio of
        # the branch that is not used may divide by zero
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio_1 = dist_x / R_1
            ratio_2 = (dist_x - R_1) / (R_2 - R_1)
            f_hngRrup = np.where(dist_rup == 0, 1.0, (dist_rup - dist_jb) / dist_rup)
        f_hngRx = np.select(
            # nshmp-lib uses f_2 for R_x = R_1
            [dist_x < 0, dist_x < R_1 if self._nshmp else dist_x <= R_1],
            [0.0, c.h_1 + c.h_2 * ratio_1 + c.h_3 * _pow(ratio_1, 2)],
            default=np.maximum(0, c.h_4 + c.h_5 * ratio_2 + c.h_6 * _pow(ratio_2, 2)),
        )

        f_hngM = np.where(
            mag <= 5.5, 0.0, np.minimum(mag - 5.5, 1) * (1 + c.a_2 * (mag - 6.5))
        )

        f_hngZ = np.where(depth_tor > 16.66, 0.0, 1 - 0.06 * depth_tor)
        f_hngDip = (90 - dip) / 45

        f_hng = c.c_10 * f_hngRx * f_hngRrup * f_hngM * f_hngZ * f_hngDip

        if not (self._nshmp or self._basin or self._cybershake):
            site_term = self._calc_site_term(
                c,
                model.as_column(pga_ref),
                model.as_column(v_s30),
                model.as_column(s.depth_2_5),
                region,
            )
        else:
            site_term = self._calc_site_term(
                c,
                model.as_column(pga_ref) if pga_rock is None else pga_rock,
                model.as_column(v_s30),
                None,
                region,
                basin=False,
            ) + self._calc_basin_term(
                c, model.as_column(v_s30), model.as_column(s.depth_2_5), region, pga_ref
            )

        # Hypocentral depth term
        f_hypH = np.clip(model.as_column(s.depth_hyp) - 7, 0, 13)
        f_hypM = c.c_17 + (c.c_18 - c.c_17) * np.clip(mag - 5.5, 0, 1)
        f_hyp = f_hypH * f_hypM

        # Fault dip term
        f_dip = c.c_19 * dip * np.clip(5.5 - mag, 0, 1)

        # Anaelastic attenuation term
        dc_20 = np.select(
            [
                model.equals(region, "japan") | model.equals(region, "italy"),
                model.equals(region, "china"),
            ],
            [c.dc_20jp, c.dc_20ch],
            # 'global', 'california'
            default=c.dc_20ca,
        )

        f_atn = (c.c_20 + dc_20) * np.maximum(dist_rup - 80, 0)

        ln_resp = f_mag + f_dis + f_flt + f_hng + site_term + f_hyp + f_dip + f_atn
        return ln_resp

    @classmethod
    def calc_site_term(
        cls,
        pga_ref: float,
        v_s30: float,
        depth_2_5: float,
        region: str = "california",
    ) -> ArrayLike:
        """Calculate the site term, which includes site and basin effects.

        Parameters
        ----------
        pga_ref : float or array_like
            peak ground acceleration [g] at the reference condition. If
            :class:`np.nan`, then the basin depth is estimated from `v_s30`.
        v_s30 : float or array_like
            site condition. Set `v_s30` to the reference
            velocity (e.g., 1180 m/s) for the reference response.
        depth_2_5 : float or array_like
            depth to the 2.5 km∕s shear-wave velocity horizon beneath the site,
            :math:`Z_{2.5}` in [km].
        region : str or array_like, optional
            region of basin model. Valid options:
            "global", "california", "japan", "italy", "china". If
            *None*, then "global" is used as the default value.

        Returns
        -------
        site_term: :class:`np.ndarray`
            site term that is applied to the natural log response. Arrays of
            inputs need a trailing axis (e.g., shape (N, 1)) to broadcast
            against the periods.
        """
        return cls._calc_site_term(cls.COEFF, pga_ref, v_s30, depth_2_5, region)

    @classmethod
    def _calc_site_term(
        cls, c, pga_ref, v_s30, depth_2_5, region, basin: bool = True
    ) -> ArrayLike:
        """Calculate the site term for the periods in the coefficients `c`.

        If `basin` is *False*, the basin response term is not included.
        """
        # Site term
        vs_ratio = v_s30 / c.k_1
        with np.errstate(invalid="ignore"):
            # The nonlinear branch is NaN for the reference condition, but the
            # reference velocity is greater than k_1 for PGA
            f_site = np.where(
                v_s30 <= c.k_1,
                c.c_11 * np.log(vs_ratio)
                + c.k_2
                * (
                    np.log(pga_ref + cls.COEFF_C * vs_ratio**cls.COEFF_N)
                    - np.log(pga_ref + cls.COEFF_C)
                ),
                (c.c_11 + c.k_2 * cls.COEFF_N) * np.log(vs_ratio),
            )

        is_japan = model.equals(region, "japan")
        if np.any(is_japan):
            # Apply regional correction for Japan
            f_site_jp = np.where(
                v_s30 <= 200,
                (c.c_12 + c.k_2 * cls.COEFF_N)
                * (np.log(vs_ratio) - np.log(200 / c.k_1)),
                (c.c_13 + c.k_2 * cls.COEFF_N) * np.log(vs_ratio),
            )
            f_site = f_site + np.where(is_japan, f_site_jp, 0.0)
        if not basin:
            return f_site

        # Basin response term
        is_ref = np.isnan(pga_ref)
        if np.any(is_ref):
            # Use model to compute depth_2_5 for the reference velocity case
            depth_2_5_ref = cls.calc_depth_2_5(v_s30, region)
            if depth_2_5 is None or np.all(is_ref):
                depth_2_5 = depth_2_5_ref
            else:
                depth_2_5 = np.where(is_ref, depth_2_5_ref, depth_2_5)

        return f_site + cls._calc_f_sed(c, depth_2_5, is_japan)

    @staticmethod
    def _calc_f_sed(c, depth_2_5, is_japan, slope=None) -> np.ndarray:
        """Calculate the basin response term (equation 20).

        `slope` replaces :math:`c_{16} k_3` (the CyberShake basin term).
        """
        f_sed_shallow = c.c_14 * (depth_2_5 - 1)
        if np.any(is_japan):
            f_sed_shallow = f_sed_shallow + np.where(
                is_japan, c.c_15 * (depth_2_5 - 1), 0.0
            )
        if slope is None:
            slope = c.c_16 * c.k_3
        return np.select(
            [depth_2_5 <= 1, depth_2_5 <= 3],
            [f_sed_shallow, 0.0],
            default=slope * np.exp(-0.75) * (1 - np.exp(-0.25 * (depth_2_5 - 3))),
        )

    def _calc_basin_term(self, c, v_s30, depth_2_5, region, pga_ref) -> np.ndarray:
        """Basin response term with the nshmp-lib options.

        The reference condition (`pga_ref` of NaN) uses the reference depth of
        the published model or, with ``nshmp`` (and without ``basin``), 0.398 km.
        """
        periods = c.period
        is_japan = model.equals(region, "japan")
        # Reference depth from V_s30 -- equation 33
        depth_ref = self.calc_depth_2_5(v_s30, region)
        term_ref = self._calc_f_sed(c, depth_ref, is_japan)
        if np.all(np.isnan(pga_ref)):
            # Reference rock condition
            if self._nshmp and not self._basin:
                return self._calc_f_sed(c, self.DEPTH_2_5_REF_NSHMP, is_japan)
            return term_ref
        cybershake = self._cybershake_rows(periods)
        slope = None
        if np.any(cybershake):
            # CyberShake basin coefficients
            c_cy = self._coeff_rows(self.COEFF_CYBERSHAKE)
            slope = np.where(cybershake, c_cy.slope, c.c_16 * c.k_3)
        with np.errstate(invalid="ignore"):
            term = self._calc_f_sed(c, depth_2_5, is_japan, slope)
        if self._basin:
            # USGS deep basin model: from the reference term at 1 km to the
            # full term at 3 km
            upper, lower = usgs.BASIN_Z2P5
            in_basin = (periods > 0.5) & (depth_2_5 > upper)
            scale = (np.clip(depth_2_5, upper, lower) - upper) / (lower - upper)
            scaled = term_ref * (1.0 - scale) + term * scale
            scaled = np.where(cybershake, scaled + usgs.CYBERSHAKE_CONSTANT, scaled)
            scaled = np.where(periods == 0.75, scaled * 0.585, scaled)
            term = np.where(in_basin, scaled, term_ref)
        if self._nshmp:
            # Reference term for NaN depths
            term = np.where(np.isnan(depth_2_5), term_ref, term)
        return term

    def _calc_ln_std(
        self, pga_ref: ArrayLike
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate the logarithmic standard deviation.

        Parameters
        ----------
        pga_ref : float or array_like
            peak ground acceleration [g] at the reference
            condition.

        Returns
        -------
            class:`np.array`: Logarithmic standard deviation.

        """
        c = self._coeff_rows(self._coeff)
        # The standard deviation of PGA is needed for all periods
        c_pga = self._coeff[[self.INDEX_PGA]]
        s = self._scenario
        mag = model.as_column(s.mag)
        v_s30 = model.as_column(s.v_s30)
        pga_ref = model.as_column(pga_ref)

        def interp_mag(value_1, value_2):
            return value_2 + (value_1 - value_2) * np.clip(5.5 - mag, 0, 1)

        tau_lnY = interp_mag(c.tau_1, c.tau_2)
        phi_lnY = interp_mag(c.phi_1, c.phi_2)

        vs_ratio = v_s30 / c.k_1
        alpha = np.where(
            v_s30 < c.k_1,
            c.k_2
            * pga_ref
            * (
                (pga_ref + self.COEFF_C * vs_ratio**self.COEFF_N) ** (-1)
                - _pow(pga_ref + self.COEFF_C, -1)
            ),
            0.0,
        )

        tau_lnPGA = interp_mag(c_pga.tau_1, c_pga.tau_2)
        tau = np.sqrt(
            tau_lnY**2
            + alpha**2 * _pow(tau_lnPGA, 2)
            + 2 * alpha * c.rho_lnPGAlnY * tau_lnY * tau_lnPGA
        )

        phi_lnPGA = interp_mag(c_pga.phi_1, c_pga.phi_2)
        phi_lnAF_PGA = c_pga.phi_lnAF
        phi_lnPGA_B = np.sqrt(_pow(phi_lnPGA, 2) - _pow(phi_lnAF_PGA, 2))
        phi_lnY_B = np.sqrt(phi_lnY**2 - c.phi_lnAF**2)

        phi = np.sqrt(
            phi_lnY_B**2
            + c.phi_lnAF**2
            + alpha**2 * (_pow(phi_lnPGA, 2) - _pow(phi_lnAF_PGA, 2))
            + 2 * alpha * c.rho_lnPGAlnY * phi_lnY_B * phi_lnPGA_B
        )

        ln_std = np.sqrt(phi**2 + tau**2)

        return ln_std, tau, phi

    @staticmethod
    def calc_depth_2_5(
        v_s30: float, region: str = "global", depth_1_0: float | None = None
    ) -> float:
        """Calculate the depth to a shear-wave velocity of 2.5 km/sec.

        This is :math:`Z_{2.5}`.

        Provide either `v_s30` or `depth_1_0`.

        Parameters
        ----------
        v_s30 : Optional[float or array_like]
            time-averaged shear-wave velocity over
            the top 30 m of the site (:math:`V_{s30}`, m/s).
            Keyword Args:
        region : Optional[str or array_like]
            region of the basin model. Valid values:
            "california", "japan". (Default value = 'global')
        depth_1_0 : Optional[float or array_like]
            depth to the 1.0 km∕s shear-wave
            velocity horizon beneath the site, :math:`Z_{1.0}` in [km].
            (Default value = None)

        Returns
        -------
        float or :class:`np.ndarray`
            estimated depth to a shear-wave velocity of 2.5 km/sec
            [km]. Arrays of `v_s30` (or `depth_1_0`) and `region` are
            broadcast against each other.

        """
        is_japan = model.equals(region, "japan")
        if v_s30 is not None and (np.ndim(v_s30) > 0 or v_s30):
            param = v_s30
            # From equations 6.10 (Japan) and 6.9 on page 63
            intercept = np.where(is_japan, 5.359, 7.089)
            slope = np.where(is_japan, 1.102, 1.144)

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 6.510
            # bar = 1.181
        elif depth_1_0 is not None and (np.ndim(depth_1_0) > 0 or depth_1_0):
            param = depth_1_0
            # From equations 6.13 (Japan) and 6.12 on page 64
            intercept = np.where(is_japan, 0.408, 1.392)
            slope = np.where(is_japan, 1.745, 1.798)

            # Global model
            # Not supported by NGA-West2 spreadsheet, and therefore removed.
            # foo = 0.748
            # bar = 2.128
        else:
            raise NotImplementedError

        return np.exp(intercept - slope * np.log(param))

    @staticmethod
    def _calc_depth_tor(mag: ArrayLike, mechanism: ArrayLike) -> ArrayLike:
        """Estimate the depth to the top of rupture [km].

        Same as :meth:`pygmm.ChiouYoungs2014.calc_depth_tor`, but accepts
        arrays of `mag` and `mechanism`.
        """
        depth_tor = np.where(
            model.equals(mechanism, "RS"),
            # Reverse and reverse-oblique faulting
            2.704 - 1.226 * np.maximum(mag - 5.849, 0),
            # Combined strike-slip and normal faulting
            2.673 - 1.136 * np.maximum(mag - 4.970, 0),
        )
        return _pow(np.maximum(depth_tor, 0), 2)

    @staticmethod
    def calc_depth_hyp(
        mag: float, dip: float, depth_tor: float, depth_bor: float
    ) -> float:
        r"""Estimate the depth to hypocenter.

        Parameters
        ----------
        mag : float or array_like
            moment magnitude of the event (:math:`M_w`)
        dip : float or array_like
            fault dip angle (:math:`\\phi`, deg).
        depth_tor : float or array_like
            depth to the top of the rupture
            plane (:math:`Z_{tor}`, km).
        depth_bor : float or array_like
            depth to the bottom of the rupture
            plane (:math:`Z_{bor}`, km).

        Returns
        -------
        float or :class:`np.ndarray`
            estimated hypocenter depth [km]

        """
        # Equations 35, 36, and 37 of journal article
        ln_dZ = _min(
            _min(-4.317 + 0.984 * mag, 2.325) + _min(0.0445 * (dip - 40), 0),
            np.log(0.9 * (depth_bor - depth_tor)),
        )

        depth_hyp = depth_tor + np.exp(ln_dZ)

        return depth_hyp

    @staticmethod
    def calc_width(
        mag: float, dip: float, depth_tor: float, depth_bot: float = 15.0
    ) -> float:
        r"""Estimate the fault width using Equation (39) of CB14.

        Parameters
        ----------
        mag : float or array_like
            moment magnitude of the event (:math:`M_w`)
        dip : float or array_like
            fault dip angle (:math:`\\phi`, deg).
        depth_tor : float or array_like
            depth to the top of the rupture
            plane (:math:`Z_{tor}`, km).
            Keyword Args:
        depth_bot : Optional[float or array_like]
            depth to bottom of seismogenic crust
            [km]. Used to calculate fault width if none is specified. If
            *None*, then a value of 15 km is used. (Default value = 15.0)

        Returns
        -------
        float or :class:`np.ndarray`
            estimated fault width [km]
        """
        return _min(
            np.sqrt(10 ** ((mag - 4.07) / 0.98)),
            (depth_bot - depth_tor) / np.sin(np.radians(dip)),
        )

    @staticmethod
    def calc_depth_bor(depth_tor: float, dip: float, width: float) -> float:
        r"""Compute the depth to bottom of the rupture [km].

        Parameters
        ----------
        dip : float or array_like
            fault dip angle (:math:`\\phi`, deg).
        depth_tor : float or array_like
            depth to the top of the rupture
            plane (:math:`Z_{tor}`, km).
        width : float or array_like
            Down-dip width of the fault.

        Returns
        -------
        float or :class:`np.ndarray`
            depth to bottom of the fault rupture [km]
        """
        return depth_tor + width * np.sin(np.radians(dip))
