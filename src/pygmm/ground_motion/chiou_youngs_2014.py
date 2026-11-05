"""Chiou and Youngs (2014, :cite:p:`chiou14`) model."""

import logging

import numpy as np

from .. import model
from .._types import ArrayLike
from ..registry import register
from . import nshmp as usgs

__author__ = "Albert Kottke"


@register(tectonic="active_crustal")
class ChiouYoungs2014(usgs.NshmpOptions, model.GroundMotionModel):
    r"""Chiou and Youngs (2014, :cite:p:`chiou14`) model.

    This model was developed for active tectonic regions as part of the
    NGA-West2 effort.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``dist_jb``, ``dist_x``, ``dip``, ``v_s30``, ``depth_tor``, ``depth_1_0``,
    ``dpp_centered``, ``mechanism``, ``on_hanging_wall``, ``region``, and
    ``vs_source``) can be a scalar or an array, and the arrays are broadcast
    against each other. If ``depth_tor`` or ``depth_1_0`` is not provided, it
    is estimated for each scenario from the magnitude and mechanism, or from
    ``v_s30`` and the region. For a scalar scenario, the response and standard
    deviation have one value per period, as in other models. For arrays of N
    scenarios, they have shape (N, periods), so, for example, ``pga`` has
    shape (N,) and ``spec_accels`` has shape (N, 24).

    The options select the variants of the model in the USGS nshmp-lib library
    (``ChiouYoungs_2014``, commit 44728a7d; see
    :mod:`pygmm.ground_motion.nshmp`), and are off by default, which gives the
    published model. With ``nshmp=True``, the model is computed as in
    nshmp-lib, which:

    - applies the hanging-wall term for sites with :math:`R_x \ge 0` if
      ``on_hanging_wall`` is not given,
    - has no basin term if ``depth_1_0`` is *None* or NaN, and always uses the
      California reference :math:`Z_{1.0}`,
    - has no directivity term (``dpp_centered`` is not used),
    - and uses the aleatory variability for an inferred :math:`V_{S30}` if
      ``vs_source`` is not given.

    nshmp-lib only has the California model, so ``region`` should not be given
    with ``nshmp=True``. The nshmp-lib ``Gmm`` ids are (see :attr:`GMM_IDS`,
    whose values include the ``vs_source`` scenario value):

    ========================  =========================================
    nshmp-lib ``Gmm``         Options
    ========================  =========================================
    ``CY_14``                 ``nshmp=True, epistemic=True``
    ``CY_14_BASE``            ``nshmp=True``
    ``CY_14_BASIN``           ``nshmp=True, epistemic=True, basin=True``
    ``CY_14_CYBERSHAKE``      ``nshmp=True, epistemic=True, basin=True,
                              cybershake=True``
    ``CY_14_VS30_MEASURED``   ``nshmp=True`` (``vs_source="measured"``)
    ``CY_14_PRVI``            ``nshmp=True, prvi=True``
    ========================  =========================================

    The other ids have ``vs_source="inferred"``. The PRVI coefficients are NaN
    at 0.04, 0.12, and 0.17 s, which nshmp-lib does not provide.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods), and/or "psa_all" (all 24 periods). Computing only
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
    ...     dist_jb=np.array([8.0, 48.0]), dist_x=np.array([5.0, -20.0]),
    ...     dip=np.array([45.0, 90.0]), v_s30=300.0,
    ...     mechanism=np.array(["RS", "SS"]),
    ...     on_hanging_wall=np.array([True, False]))
    >>> pygmm.ChiouYoungs2014(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.ChiouYoungs2014(s, ims=["pga", "psa_ngawest2_21"])
    >>> m.spec_accels.shape
    (2, 21)

    The nshmp-lib ``CY_14_BASIN`` variant (with the hanging wall given by
    ``dist_x``), with the ``vs_source`` of the id in the scenario:

    >>> options = dict(pygmm.ChiouYoungs2014.GMM_IDS["CY_14_BASIN"])
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 50.0]),
    ...     dist_jb=np.array([8.0, 48.0]), dist_x=np.array([8.0, -48.0]),
    ...     dip=60.0, depth_tor=2.0, mechanism="RS", v_s30=400.0,
    ...     depth_1_0=np.array([np.nan, 0.6]), vs_source=options.pop("vs_source"))
    >>> pygmm.ChiouYoungs2014(s, **options).spec_accels.shape
    (2, 24)
    """

    NAME = "Chiou and Youngs (2014)"
    ABBREV = "CY14"

    # Reference velocity (m/s)
    V_REF = 1130.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("chiou_youngs_2014.csv", 2)

    PERIODS = COEFF["period"]

    #: nshmp-lib PRVI coefficients (NaN at 0.04, 0.12, and 0.17 s)
    COEFF_PRVI = usgs.load_coeffs("chiou_youngs_2014-prvi.csv", PERIODS)
    #: nshmp-lib CyberShake basin coefficients
    COEFF_CYBERSHAKE = usgs.load_coeffs("chiou_youngs_2014-cybershake.csv", PERIODS)

    #: Model options and scenario values of the nshmp-lib ``Gmm`` ids
    GMM_IDS = {
        "CY_14": dict(nshmp=True, epistemic=True, vs_source="inferred"),
        "CY_14_BASE": dict(nshmp=True, vs_source="inferred"),
        "CY_14_BASIN": dict(
            nshmp=True, epistemic=True, basin=True, vs_source="inferred"
        ),
        "CY_14_CYBERSHAKE": dict(
            nshmp=True,
            epistemic=True,
            basin=True,
            cybershake=True,
            vs_source="inferred",
        ),
        "CY_14_VS30_MEASURED": dict(nshmp=True, vs_source="measured"),
        "CY_14_PRVI": dict(nshmp=True, prvi=True, vs_source="inferred"),
    }

    #: Reference velocity (m/sec) of the reference Z_1.0 in nshmp-lib
    V_REF_Z1_NSHMP = 571.0

    INDICES_PSA = np.arange(24)
    INDEX_PGA = 24
    INDEX_PGV = 25

    PARAMS = [
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("dist_x", True),
        model.NumericParameter("dist_jb", True),
        model.NumericParameter("mag", True, 3.5, 8.5),
        model.NumericParameter("v_s30", True, 180, 1500),
        model.NumericParameter("depth_tor", False, 0, 20),
        model.NumericParameter("depth_1_0", False),
        model.NumericParameter("dpp_centered", False, default=0),
        model.NumericParameter("dip", True),
        model.CategoricalParameter("mechanism", False, ["U", "SS", "NS", "RS"], "U"),
        model.CategoricalParameter("on_hanging_wall", False, [True, False], False),
        model.CategoricalParameter(
            "region", False, ["california", "china", "italy", "japan"], "california"
        ),
        model.CategoricalParameter(
            "vs_source", False, ["measured", "inferred"], "measured"
        ),
    ]

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
        """Initialize the model."""
        self._set_options(nshmp, epistemic, basin, cybershake, prvi)
        super().__init__(scenario, ims)
        # The nonlinear site term of each period depends on the reference
        # response at the same period, so only the requested periods are needed
        c = self._coeff_rows(self._options_coeff(self.COEFF))
        ln_resp_ref = self._calc_ln_resp_ref(c)
        self._ln_resp = self._calc_ln_resp_site(ln_resp_ref, c)
        self._ln_std, self._tau, self._phi = self._calc_ln_std(np.exp(ln_resp_ref), c)
        if self._epistemic:
            self._ln_resp, self._branches = usgs.epistemic_branches(
                self._ln_resp,
                self._ln_std,
                model.as_column(self._scenario.mag),
                model.as_column(self._scenario.dist_jb),
            )
        if self._nshmp or self._epistemic:
            self._ln_resp, self._ln_std = usgs.broadcast(self._ln_resp, self._ln_std)

    def _calc_ln_resp_ref(self, c=None) -> np.ndarray:
        """Calculate the natural logarithm of the reference response.

        Parameters
        ----------
        c : :class:`numpy.recarray`, optional
            coefficients for the periods to compute. If *None*, the
            coefficients for the requested intensity measures are used.

        Returns
        -------
        ln_resp_ref : class:`np.array`:
            natural log of the response

        """
        if c is None:
            c = self._coeff_rows(self.COEFF)
        s = self._scenario
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)
        mechanism = model.as_column(s.mechanism)
        region = model.as_column(s.region)
        cos_dip = np.cos(np.radians(model.as_column(s.dip)))

        cosh_mag = np.cosh(2 * np.maximum(mag - 4.5, 0))
        ln_resp = np.array(c.c_1)

        # Reverse and normal fault terms. Other mechanisms have no term.
        ln_resp = ln_resp + np.select(
            [model.equals(mechanism, "RS"), model.equals(mechanism, "NS")],
            [c.c_1a + c.c_1c / cosh_mag, c.c_1b + c.c_1d / cosh_mag],
            default=0.0,
        )

        # Magnitude scaling
        ln_resp = ln_resp + c.c_2 * (mag - 6)
        ln_resp = ln_resp + (c.c_2 - c.c_3) / c.c_n * np.log(
            1 + np.exp(c.c_n * (c.c_m - mag))
        )

        # Top of rupture term relative to model average
        diff_depth_tor = model.as_column(s.depth_tor) - model.as_column(
            self.calc_depth_tor(s.mag, s.mechanism)
        )
        ln_resp = ln_resp + (c.c_7 + c.c_7b / cosh_mag) * diff_depth_tor

        # Dip angle term
        ln_resp = ln_resp + (c.c_11 + c.c_11b / cosh_mag) * cos_dip**2

        # Distance terms
        ln_resp = ln_resp + c.c_4 * np.log(
            dist_rup + c.c_5 * np.cosh(c.c_6 * np.maximum(mag - c.c_hm, 0))
        )
        ln_resp = ln_resp + (c.c_4a - c.c_4) * np.log(np.sqrt(dist_rup**2 + c.c_rb**2))

        # Regional adjustment
        scale = np.select(
            [
                (model.equals(region, "japan") | model.equals(region, "italy"))
                & (6 < mag)
                & (mag < 6.9),
                model.equals(region, "china"),
            ],
            [c.gamma_ji, c.gamma_c],
            default=1.0,
        )
        ln_resp = ln_resp + (
            scale
            * (c.c_gamma1 + c.c_gamma2 / np.cosh(np.maximum(mag - c.c_gamma3, 0)))
            * dist_rup
        )

        # Directivity term, which is not in nshmp-lib
        if not self._nshmp:
            ln_resp = ln_resp + (
                c.c_8
                * np.maximum(1 - np.maximum(dist_rup - 40, 0) / 30, 0)
                * np.minimum(np.maximum(mag - 5.5, 0) / 0.8, 1)
                * np.exp(-c.c_8a * (mag - c.c_8b) ** 2)
                * model.as_column(s.dpp_centered)
            )

        # Hanging wall term, which is only applied to sites on the hanging wall
        on_hanging_wall = model.as_column(model.equals(s.on_hanging_wall, True))
        if np.any(on_hanging_wall):
            ln_resp = ln_resp + np.where(
                on_hanging_wall,
                c.c_9
                * cos_dip
                * (c.c_9a + (1 - c.c_9a) * np.tanh(model.as_column(s.dist_x) / c.c_9b))
                * (
                    1
                    - np.sqrt(
                        model.as_column(s.dist_jb) ** 2
                        + model.as_column(s.depth_tor) ** 2
                    )
                    / (dist_rup + 1)
                ),
                0.0,
            )

        return ln_resp

    def _calc_ln_resp_site(self, ln_resp_ref: np.ndarray, c=None) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        ln_resp_ref : :class:`np.array`
            natural logarithm of the response at the reference site condition
            at each of the periods specified by the model coefficients.
        c : :class:`numpy.recarray`, optional
            coefficients for the periods to compute. If *None*, the
            coefficients for the requested intensity measures are used.

        Returns
        -------
        ln_resp_rsite : class:`np.ndarray`:
            natural log of the response including the site effects

        """
        if c is None:
            c = self._coeff_rows(self.COEFF)
        s = self._scenario

        v_s30 = model.as_column(s.v_s30)
        depth_1_0 = model.as_column(s.depth_1_0)
        region = s.region if np.ndim(s.region) == 0 else model.as_column(s.region)
        if not (self._nshmp or self._basin or self._cybershake):
            site_term = self._calc_site_term(
                c, np.exp(ln_resp_ref), v_s30, depth_1_0, region
            )
        else:
            site_term = self._calc_site_term(
                c, np.exp(ln_resp_ref), v_s30, None, region
            ) + self._calc_basin_term(c, v_s30, depth_1_0, region)
        return ln_resp_ref + site_term

    def _calc_basin_term(self, c, v_s30, depth_1_0, region) -> np.ndarray:
        """Basin (sediment thickness) term with the nshmp-lib options."""
        periods = c.period
        _, phi_5, phi_6 = self._region_coeffs(c, region)
        cybershake = self._cybershake_rows(periods)
        if np.any(cybershake):
            # CyberShake basin coefficients
            c_cy = self._coeff_rows(self.COEFF_CYBERSHAKE)
            phi_5 = np.where(cybershake, c_cy.phi_5, phi_5)
            phi_6 = np.where(cybershake, c_cy.phi_6, phi_6)
        if self._nshmp:
            # Reference Z_1.0 (m) with 571 m/sec, and no basin term for NaN
            # depths
            v_s30_pow4 = v_s30 * v_s30 * v_s30 * v_s30
            v_ref_pow4 = self.V_REF_Z1_NSHMP**4
            depth_ref = np.exp(
                -7.15 / 4 * np.log((v_s30_pow4 + v_ref_pow4) / (1360.0**4 + v_ref_pow4))
            )
            diff_depth_1_0 = np.where(
                np.isnan(depth_1_0), 0.0, depth_1_0 * 1000.0 - depth_ref
            )
        else:
            diff_depth_1_0 = 1000 * (depth_1_0 - self.calc_depth_1_0(v_s30, region))
        term = phi_5 * (1.0 - np.exp(-diff_depth_1_0 / phi_6))
        term = np.where(cybershake, term + usgs.CYBERSHAKE_CONSTANT, term)
        if self._basin:
            term = term * usgs.basin_scale(periods, depth_1_0, *usgs.BASIN_Z1P0)
        return term

    @classmethod
    def calc_site_term(
        cls, resp_ref: ArrayLike, v_s30: float, depth_1_0: float, region: str | None
    ) -> ArrayLike:
        """Calculate the site term, which includes site and basin effects.

        Parameters
        ----------
        resp_ref : :class:`np.array`
            response (not its natural logarithm) at the reference site
            condition at each of the periods specified by the model
            coefficients.
        v_s30 : float or array_like
            site condition. Set `v_s30` to the reference
            velocity (e.g., 1180 m/s) for the reference response.
        depth_1_0 : float or array_like
            depth to the 1.0 km∕s shear-wave velocity horizon beneath the site,
            :math:`Z_{1.0}` in [km]. If *None*, then no basin term is applied.
        region : str or array_like, optional
            region of basin model. Valid options: 'california', 'japan'. If
            *None*, then 'california' is used as the default value.

        Returns
        -------
        site_term: :class:`np.ndarray`
            site term that is applied to the natural log response.

        Notes
        -----
        Arrays of scenario values broadcast against the coefficients, which
        have one value per period, so N scenarios need values with shape
        (N, 1) (see :func:`pygmm.model.as_column`).
        """
        return cls._calc_site_term(cls.COEFF, resp_ref, v_s30, depth_1_0, region)

    @classmethod
    def _calc_site_term(cls, c, resp_ref, v_s30, depth_1_0, region) -> ArrayLike:
        """Calculate the site term for the periods in the coefficients `c`."""
        phi_1, phi_5, phi_6 = cls._region_coeffs(c, region)

        ln_resp_ref = np.log(resp_ref)

        # Linear scaling
        site_term = phi_1 * np.minimum(np.log(v_s30 / 1130.0), 0)
        # Nonlinear scaling
        site_term = site_term + (
            c.phi_2
            * (
                np.exp(c.phi_3 * (np.minimum(v_s30, 1130.0) - 360.0))
                - np.exp(c.phi_3 * (1130.0 - 360.0))
            )
            * np.log((np.exp(ln_resp_ref) + c.phi_4) / c.phi_4)
        )
        if depth_1_0 is None:
            return site_term
        # Basin model
        diff_depth_1_0 = 1000 * (depth_1_0 - cls.calc_depth_1_0(v_s30, region))
        site_term = site_term + phi_5 * (1 - np.exp(-diff_depth_1_0 / phi_6))

        return site_term

    @staticmethod
    def _region_coeffs(c, region) -> tuple:
        """Coefficients phi_1, phi_5, and phi_6 of the basin region."""
        if np.ndim(region) > 0:
            is_japan = model.equals(region, "japan")
            return (
                np.where(is_japan, c.phi_1jp, c.phi_1),
                np.where(is_japan, c.phi_5jp, c.phi_5),
                np.where(is_japan, c.phi_6jp, c.phi_6),
            )
        if region in ["japan"]:
            return c.phi_1jp, c.phi_5jp, c.phi_6jp
        return c.phi_1, c.phi_5, c.phi_6

    def _calc_ln_std(
        self, resp_ref: np.ndarray, c=None
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Calculate the logarithmic standard deviation.

        Parameters
        ----------
        resp_ref : :class:`np.array`
            response at the reference site condition at each of the computed
            periods.
        c : :class:`numpy.recarray`, optional
            coefficients for the periods to compute. If *None*, the
            coefficients for the requested intensity measures are used.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation
        tau : class:`np.array`:
            between-event standard deviation, including the nonlinear site term
        phi : class:`np.array`:
            within-event standard deviation

        """
        if c is None:
            c = self._coeff_rows(self.COEFF)
        s = self._scenario

        sigma_2 = np.where(
            model.as_column(model.equals(s.region, "japan")), c.sigma_2jp, c.sigma_2
        )

        clipped_mag = np.clip(model.as_column(s.mag), 5.0, 6.5) - 5.0
        tau = c.tau_1 + (c.tau_2 - c.tau_1) / 1.5 * clipped_mag

        nl_0 = (
            c.phi_2
            * (
                np.exp(c.phi_3 * (np.minimum(model.as_column(s.v_s30), 1130.0) - 360.0))
                - np.exp(c.phi_3 * (1130.0 - 360.0))
            )
            * (resp_ref / (resp_ref + c.phi_4))
        )

        flag_meas = model.as_column(
            np.where(model.equals(s.vs_source, "measured"), 1, 0)
        )
        phi_nl = (c.sigma_1 + (sigma_2 - c.sigma_1) / 1.5 * clipped_mag) * np.sqrt(
            c.sigma_3 * (1 - flag_meas) + 0.7 * flag_meas + (1 + nl_0) ** 2
        )

        ln_std = np.sqrt((1 + nl_0) ** 2 * tau**2 + phi_nl**2)
        return ln_std, np.sqrt((1 + nl_0) ** 2 * tau**2), phi_nl

    def _check_inputs(self) -> None:
        """Check the inputs."""
        s = self._scenario
        if self._nshmp:
            # nshmp-lib defaults
            if s["on_hanging_wall"] is None and s["dist_x"] is not None:
                s["on_hanging_wall"] = (np.asarray(s["dist_x"]) >= 0)[()]
            if s["vs_source"] is None:
                s["vs_source"] = "inferred"
            if s["depth_1_0"] is None:
                # No basin term
                s["depth_1_0"] = np.nan
        super()._check_inputs()

        if np.ndim(s.mechanism) > 0 or np.ndim(s.mag) > 0:
            self._check_mag_by_mechanism()
        else:
            if s.mechanism in ["RS", "NS"]:
                _min, _max = 3.5, 8.0
            else:
                _min, _max = 3.5, 8.5

            if not (_min <= s.mag <= _max):
                logging.warning(
                    "Magnitude (%g) exceeds recommended bounds (%g to %g)"
                    " for a %s earthquake!",
                    s.mag,
                    _min,
                    _max,
                    s.mechanism,
                )

        if s.get("depth_tor", None) is None:
            s["depth_tor"] = self.calc_depth_tor(s.mag, s.mechanism)

        if s.get("depth_1_0", None) is None:
            # Calculate depth (m) and convert to (km)
            s["depth_1_0"] = self.calc_depth_1_0(s.v_s30, s.region)

    def _check_mag_by_mechanism(self) -> None:
        """Check mechanism specific magnitude limits for vectorized scenarios."""
        s = self._scenario
        mag, mechanism = np.broadcast_arrays(np.asarray(s.mag), np.asarray(s.mechanism))
        is_rs_ns = model.equals(mechanism, "RS") | model.equals(mechanism, "NS")
        for is_option, name, _min, _max in [
            (is_rs_ns, "RS or NS", 3.5, 8.0),
            (~is_rs_ns, "U or SS", 3.5, 8.5),
        ]:
            outside = is_option & ~((_min <= mag) & (mag <= _max))
            if np.any(outside):
                logging.warning(
                    "Magnitude (%d of %d %s earthquakes, %g to %g) exceeds"
                    " recommended bounds (%g to %g)!",
                    np.count_nonzero(outside),
                    np.count_nonzero(is_option),
                    name,
                    np.min(mag[outside]),
                    np.max(mag[outside]),
                    _min,
                    _max,
                )

    @staticmethod
    def calc_depth_1_0(v_s30: float, region: str) -> float:
        """Calculate the depth to 1 km/sec (:math:`Z_{1.0}`).

        Equations 1 (California) and 2 (Japan) of Chiou and Youngs (2014).

        Parameters
        ----------
        v_s30 : float or array_like
            time-averaged shear-wave velocity over the top 30 m
            of the site (:math:`V_{s30}`, m/s).
        region : str or array_like
            basin region. Valid options: "california", "japan". An array of
            regions is broadcast against `v_s30`.

        Returns
        -------
        depth_1_0 : float
            estimated depth to a shear-wave velocity of 1 km/sec [km]

        """
        if np.ndim(region) > 0:
            return np.where(
                model.equals(region, "japan"),
                ChiouYoungs2014.calc_depth_1_0(v_s30, "japan"),
                ChiouYoungs2014.calc_depth_1_0(v_s30, "california"),
            )
        if region in ["japan"]:
            # Japan
            power = 2
            v_ref = 412.0
            slope = -5.23 / power
        else:
            # Global
            power = 4
            v_ref = 571.0
            slope = -7.15 / power

        return (
            np.exp(
                slope
                * np.log((v_s30**power + v_ref**power) / (1360.0**power + v_ref**power))
            )
            / 1000
        )

    @staticmethod
    def calc_depth_tor(mag: float, mechanism: str) -> float:
        """Calculate an estimate of the depth to top of rupture [km].

        Parameters
        ----------
        mag : float or array_like
            moment magnitude of the event (:math:`M_w`)
        mechanism : str or array_like
            fault mechanism. Valid options: "U", "SS", "NS",
            "RS". An array of mechanisms is broadcast against `mag`.

        Returns
        -------
        depth_tor : float or :class:`np.ndarray`
            estimated depth to top of rupture [km]

        """
        if np.ndim(mag) > 0 or np.ndim(mechanism) > 0:
            depth_tor = np.where(
                model.equals(mechanism, "RS"),
                # Reverse and reverse-oblique faulting
                2.704 - 1.226 * np.maximum(np.asarray(mag) - 5.849, 0),
                # Combined strike-slip and normal faulting
                2.673 - 1.136 * np.maximum(np.asarray(mag) - 4.970, 0),
            )
            return np.maximum(depth_tor, 0) ** 2

        if mechanism == "RS":
            # Reverse and reverse-oblique faulting
            depth_tor = 2.704 - 1.226 * max(mag - 5.849, 0)
        else:
            # Combined strike-slip and normal faulting
            depth_tor = 2.673 - 1.136 * max(mag - 4.970, 0)

        return max(depth_tor, 0) ** 2
