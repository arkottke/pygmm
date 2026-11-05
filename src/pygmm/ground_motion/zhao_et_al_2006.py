r"""Zhao et al. (2006) subduction model as in nshmp-lib.

The model follows the ``ZhaoEtAl_2006`` class of the USGS nshmp-lib library
(commit 44728a7d), which is used for the current U.S. Geological Survey
national seismic hazard models (NSHMs). nshmp-lib extrapolates the 7.5 and
10 s results with the BC Hydro (2012) model (``BcHydro_2012``), which is
included here only for that purpose.

References
----------
Zhao, J. X., Zhang, J., Asano, A., Ohno, Y., Oouchi, T., Takahashi, T.,
Ogawa, H., Irikura, K., Thio, H. K., Somerville, P. G., Fukushima, Y., and
Fukushima, Y. (2006). Attenuation relations of strong ground motion in Japan
using site classification based on predominant period. Bulletin of the
Seismological Society of America, 96(3), 898-913.
https://doi.org/10.1785/0120050122

Abrahamson, N., Gregor, N., and Addo, K. (2016). BC Hydro ground motion
prediction equations for subduction earthquakes. Earthquake Spectra, 32(1),
23-44. https://doi.org/10.1193/051712EQS188MR

USGS nshmp-lib, ``gov.usgs.earthquake.nshmp.gmm.ZhaoEtAl_2006``,
``BcHydro_2012``, ``ExtrapolatedGmm``, ``InterpolatedGmm``, and
``UsgsPgvSupport``, commit 44728a7d
(https://code.usgs.gov/ghsc/nshmp/nshmp-lib).
"""

import numpy as np

from .. import model
from ..registry import register
from .atkinson_macias_2009 import (
    LN_G,
    AtkinsonMacias2009,
    _NshmpSubductionBase,
    coeff_rows,
    deep_basin_term,
    interpolate_periods,
)

#: nshmp-lib BC Hydro (2012) coefficients (``BCHydro12.csv``)
_COEFF_BCHYDRO = model.load_data_file("bchydro_2012-nshmp.csv", 1)


def bchydro_2012(periods, v, basin=False) -> tuple:
    r"""BC Hydro (2012) model as in nshmp-lib (``BcHydro_2012``).

    The forearc model with the mean :math:`\Delta C_1` for interface events
    and :math:`\Delta C_1 = -0.3` for intraslab events, a reference PGA at
    :math:`V_{S30}` of 1000 m/sec, a standard deviation of 0.74, and, with
    `basin`, the USGS deep basin term (:func:`deep_basin_term`, without the M9
    adjustment, which nshmp-lib does not apply to this model). Used by
    nshmp-lib to extrapolate Zhao et al. (2006) to 7.5 and 10 s.

    Parameters
    ----------
    periods : array_like
        periods [s] with coefficients, with 0 for PGA (0.03 s is not
        included, as nshmp-lib interpolates it)
    v : dict
        ``mag``, ``dist_rup``, ``depth_tor``, ``v_s30``, ``depth_2_5``, and
        ``is_slab`` with a trailing axis for the periods
    basin : bool, optional
        add the USGS deep basin term

    Returns
    -------
    ln_resp : :class:`np.ndarray`
        natural log of the response [g]
    ln_std : :class:`np.ndarray`
        logarithmic standard deviation
    """
    periods = np.asarray(periods, dtype=float)
    c = coeff_rows(_COEFF_BCHYDRO, periods)
    c_pga = coeff_rows(_COEFF_BCHYDRO, [0.0])
    pga_rock = np.exp(_bchydro_ln_mean(c_pga, v, 0.0, 1000.0))
    ln_resp = _bchydro_ln_mean(c, v, pga_rock, v["v_s30"])
    if basin:
        ln_resp = ln_resp + deep_basin_term(periods, v["depth_2_5"])
    return ln_resp, np.full(periods.shape, 0.74)


def _bchydro_ln_mean(c, v, pga_rock, v_s30) -> np.ndarray:
    """BcHydro_2012.calcMean (forearc)."""
    t3, t4, t5, t9 = 0.1, 0.9, 0.0, 0.4
    c_1, c_4, c_, n_ = 7.8, 10.0, 1.88, 1.18
    is_slab = v["is_slab"]
    mag = v["mag"]
    dist_rup = v["dist_rup"]
    delta_c1 = np.where(is_slab, -0.3, c.dC1mid)
    m_cut = c_1 + delta_c1
    t13m = c.t13 * (10 - mag) * (10 - mag)
    f_mag = np.where(mag <= m_cut, t4, t5) * (mag - m_cut) + t13m
    # The depth term is only used for intraslab events
    with np.errstate(invalid="ignore"):
        f_depth = np.where(
            is_slab, c.t11 * (np.minimum(v["depth_tor"], 120.0) - 60.0), 0.0
        )
    v_s = np.minimum(v_s30, 1000.0)
    f_site = c.t12 * np.log(v_s / c.vlin)
    f_site = f_site + np.where(
        v_s30 < c.vlin,
        -c.b * np.log(pga_rock + c_)
        + c.b * np.log(pga_rock + c_ * (v_s / c.vlin) ** n_),
        c.b * n_ * np.log(v_s / c.vlin),
    )
    return (
        c.t1
        + t4 * delta_c1
        + (c.t2 + np.where(is_slab, c.t14, 0.0) + t3 * (mag - 7.8))
        * np.log(dist_rup + c_4 * np.exp((mag - 6.0) * t9))
        + c.t6 * dist_rup
        + np.where(is_slab, c.t10, 0.0)
        + f_mag
        + f_depth
        + f_site
    )


@register(tectonic="subduction")
class ZhaoEtAl2006(_NshmpSubductionBase):
    r"""Zhao et al. (2006) subduction model (nshmp-lib).

    Ground motion model for subduction interface and intraslab earthquakes
    (and crustal earthquakes, which are not implemented) in Japan, as
    implemented in the ``ZhaoEtAl_2006`` class of USGS nshmp-lib (commit
    44728a7d), with the coefficients of ``Zhao06.csv``. nshmp-lib differs from
    the published model as follows:

    * The site class is selected from :math:`V_{S30}`: the hard rock class
      (:math:`C_H`) is not used, and :math:`C_1` (SC I, rock), :math:`C_2` (SC
      II, hard soil), and :math:`C_3` (SC III, medium soil) are used for
      :math:`V_{S30} \geq 600`, :math:`300 \leq V_{S30} < 600`, and
      :math:`V_{S30} < 300` m/sec. The soft soil class (:math:`C_4`) is not
      used.
    * The focal depth is fixed at 20 km for interface events, and is
      :math:`\min(Z_{tor}, 125)` km (the depth to the top of rupture,
      ``depth_tor``) for intraslab events. The depth term is :math:`e (h -
      15)` for depths :math:`h` of 15 km or more, and zero otherwise.
    * The magnitude-squared corrections of Zhao et al. (2006) (Table 6,
      :math:`Q_I` and :math:`W_I` for interface events, and :math:`P_S`,
      :math:`Q_S`, and :math:`W_S` for intraslab events, with reference
      magnitudes of 6.3 and 6.5) are applied, but the standard deviation
      uses the inter-event standard deviation of Table 5 (:math:`\tau` for
      interface events and :math:`\tau_S` for intraslab events) instead of
      the reduced values for the corrected model:
      :math:`\sqrt{\sigma^2 + \tau^2}`.
    * The interface event term is :math:`S_I` (there is no :math:`S_R` term
      for crustal reverse events), and the intraslab terms are :math:`S_S +
      S_{SL} \ln(r)`.
    * The rupture distance is at least 1 km.
    * PGA coefficients are used for 0.01 s. The 0.02, 0.03, 0.075, and 0.75
      s results (including the basin term) are interpolated (linearly in
      period) between 0.01 and 0.05 s, 0.05 and 0.1 s, and 0.5 and 1.0 s
      (``InterpolatedGmm``). So the basin term at 0.75 s is half of that at
      1.0 s (instead of the 0.585 factor of the USGS basin model).
    * The 7.5 and 10 s results are extrapolated (``ExtrapolatedGmm``) from
      the 5 s result with reference models: the BC Hydro (2012) model
      (:func:`bchydro_2012`) for intraslab events, and the average of BC
      Hydro (2012) and
      :class:`~pygmm.ground_motion.atkinson_macias_2009.AtkinsonMacias2009`
      (without the site fix, with weights of 0.5) for interface events, with
      the same basin options. The natural log of the median and the standard
      deviation at 5 s are multiplied by the ratio of those of the reference
      model at 7.5 or 10 s to those at 5 s, :math:`\mu(T) = \mu_{ref}(T)
      \mu(5) / \mu_{ref}(5)`. (Note that this scales the natural log of the
      median, not the median.)
    * PGV [cm/sec] is computed with the Abrahamson and Bhasin (2020)
      conditional model from the response spectrum of the model
      (:func:`~pygmm.ground_motion.atkinson_macias_2009.conditional_pgv`).
    * The scenario values are not clipped to the model limits.

    The nshmp-lib ``Gmm`` ids correspond to the ``event_type`` scenario value
    and the options as follows (see :attr:`GMM_IDS`):

        ========================== ========== =========================
        nshmp-lib Gmm id           event_type options
        ========================== ========== =========================
        ZHAO_06_INTERFACE          interface
        ZHAO_06_INTERFACE_BASIN    interface  basin=True
        ZHAO_06_INTERFACE_BASIN_M9 interface  basin=True, m9=True
        ZHAO_06_SLAB               intraslab
        ZHAO_06_SLAB_BASIN         intraslab  basin=True
        ========================== ========== =========================

    ``ZHAO_06_INTERFACE_BASIN`` (weight of 0.125) and ``ZHAO_06_SLAB_BASIN``
    (weight of 0.25) are used by the Cascadia interface and intraslab logic
    trees of the current conterminous U.S. NSHM (nshm-conus 6.2.0). nshmp-lib
    has no intraslab M9 variant, so ``m9`` only affects interface events.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``depth_tor``, ``v_s30``, ``depth_2_5``, and ``event_type``) can be a
    scalar or an array, and the arrays are broadcast against each other, so
    one array can mix interface and intraslab events. For a scalar scenario,
    the response and standard deviation have one value per period, as in
    other models. For arrays of N scenarios, they have shape (N, periods), so,
    for example, ``pga`` has shape (N,). ``depth_tor`` is only needed for
    intraslab events, so it can be *None* if none of the events are intraslab
    events. ``depth_2_5`` [km] is only used with ``basin=True``, and can be
    *None* (or *NaN* for some scenarios) for no basin term. Results for
    events that are neither interface nor intraslab events are *NaN*.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    basin : bool, optional
        if *True*, add the USGS deep basin term (see
        :func:`~pygmm.ground_motion.atkinson_macias_2009.deep_basin_term`; the
        ``_BASIN`` variants). Default is *False*.
    m9 : bool, optional
        if *True*, apply the Seattle M9 adjustment to interface events, which
        replaces the basin term by :math:`\ln 2` for :math:`Z_{2.5}` greater
        than 6 km and periods of 2 s and longer
        (``ZHAO_06_INTERFACE_BASIN_M9``). It requires ``basin=True``. Default
        is *False*.
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods, which are all of the periods), and/or "psa_all". If *None*
        (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> from pygmm.model import Scenario
    >>> from pygmm.ground_motion.zhao_et_al_2006 import ZhaoEtAl2006
    >>> s = Scenario(
    ...     mag=np.array([7.0, 9.0]), dist_rup=np.array([80.0, 100.0]),
    ...     depth_tor=50.0, v_s30=400.0, depth_2_5=np.array([np.nan, 4.0]),
    ...     event_type=np.array(["intraslab", "interface"]))
    >>> m = ZhaoEtAl2006(s, basin=True)
    >>> m.pga.shape, m.spec_accels.shape
    ((2,), (2, 21))
    >>> ZhaoEtAl2006(s, ims=["pga", "psa_1p000"]).spec_accels.shape
    (2, 1)

    The ``GMM_IDS`` values include the ``event_type`` scenario value:

    >>> options = dict(ZhaoEtAl2006.GMM_IDS["ZHAO_06_SLAB_BASIN"])
    >>> event_type = options.pop("event_type")
    >>> m = ZhaoEtAl2006(s.copy_with(event_type=event_type), **options)
    >>> m.basin, m.m9
    (True, False)
    """

    NAME = "Zhao et al. (2006)"
    ABBREV = "Zea06"

    COEFF = model.load_data_file("zhao_et_al_2006-nshmp.csv", 1)

    INTERPOLATED = {
        0.02: (0.01, 0.05),
        0.03: (0.01, 0.05),
        0.075: (0.05, 0.1),
        0.75: (0.5, 1.0),
    }
    #: Periods (s) that nshmp-lib extrapolates from 5 s
    EXTRAPOLATED = (7.5, 10.0)

    #: Model options and scenario values of the nshmp-lib ``Gmm`` ids
    GMM_IDS = {
        "ZHAO_06_INTERFACE": dict(event_type="interface"),
        "ZHAO_06_INTERFACE_BASIN": dict(event_type="interface", basin=True),
        "ZHAO_06_INTERFACE_BASIN_M9": dict(event_type="interface", basin=True, m9=True),
        "ZHAO_06_SLAB": dict(event_type="intraslab"),
        "ZHAO_06_SLAB_BASIN": dict(event_type="intraslab", basin=True),
    }

    LIMITS = dict(mag=(5.0, 9.5), dist_rup=(0.0, 1000.0), v_s30=(150.0, 1000.0))

    PARAMS = [
        model.NumericParameter("mag", True, 5.0, 9.5),
        model.NumericParameter("dist_rup", True, 0.0, 1000.0),
        # Only used for intraslab events
        model.NumericParameter("depth_tor", False),
        model.NumericParameter("v_s30", True, 150.0, 1000.0),
        model.NumericParameter("depth_2_5", False),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
    ]

    # Constants
    _HC = 15.0
    _MC_S = 6.5
    _MC_I = 6.3
    _MAX_SLAB_DEPTH = 125.0
    _INTERFACE_DEPTH = 20.0

    def __init__(
        self,
        scenario: model.Scenario,
        basin: bool = False,
        m9: bool = False,
        ims=None,
    ):
        """Initialize the model."""
        if m9 and not basin:
            raise ValueError(
                "m9=True requires basin=True (nshmp-lib ZHAO_06_INTERFACE_BASIN_M9)"
            )
        super().__init__(scenario, ims)
        self._basin = bool(basin)
        self._m9 = bool(m9)
        self._ln_resp, self._ln_std = self._calc()

    @property
    def basin(self) -> bool:
        """If the USGS deep basin term is added."""
        return self._basin

    @property
    def m9(self) -> bool:
        """If the Seattle M9 adjustment is applied to interface events."""
        return self._m9

    def _scenario_values(self) -> dict:
        s = self._scenario
        event_type = model.as_column(s.event_type)
        is_slab = model.equals(event_type, "intraslab")
        is_interface = model.equals(event_type, "interface")
        if s.depth_tor is None:
            if np.any(is_slab):
                raise ValueError("depth_tor is required for intraslab events")
            # Not used for interface events
            depth_tor = np.array([np.nan])
        else:
            depth_tor = model.as_column(s.depth_tor).astype(float)
        return dict(
            mag=model.as_column(s.mag).astype(float),
            dist_rup=model.as_column(s.dist_rup).astype(float),
            v_s30=model.as_column(s.v_s30).astype(float),
            depth_tor=depth_tor,
            depth_2_5=self._depth_2_5(),
            is_slab=is_slab,
            is_interface=is_interface,
            valid=is_slab | is_interface,
        )

    def _calc_ln_base(self, periods, v) -> tuple:
        """Ground motions at the periods with coefficients (``calcMean``)."""
        periods = np.asarray(periods, dtype=float)
        c = coeff_rows(self.COEFF, periods)
        is_slab = v["is_slab"]
        mag = v["mag"]
        # Avoid ln(0)
        dist_rup = np.maximum(v["dist_rup"], 1.0)
        with np.errstate(invalid="ignore"):
            depth = np.where(
                is_slab,
                np.minimum(v["depth_tor"], self._MAX_SLAB_DEPTH),
                self._INTERFACE_DEPTH,
            )
            h_fac = np.where(depth < self._HC, 0.0, depth - self._HC)
        m2 = mag - np.where(is_slab, self._MC_S, self._MC_I)
        a_fac = np.where(is_slab, c.Ssl * np.log(dist_rup) + c.Ss, c.Si)
        xm_cor = np.where(
            is_slab, c.Ps * m2 + c.Qs * m2 * m2 + c.Ws, c.Qi * m2 * m2 + c.Wi
        )
        # Site term from V_s30 (siteTermStep)
        v_s30 = v["v_s30"]
        f_site = np.where(v_s30 >= 600.0, c.C1, np.where(v_s30 >= 300.0, c.C2, c.C3))
        r = dist_rup + c.c * np.exp(c.d * mag)
        ln_resp = (
            c.a * mag
            + c.b * dist_rup
            - np.log(r)
            + c.e * h_fac
            + a_fac
            + f_site
            + xm_cor
            - LN_G
        )
        if self._basin:
            f_basin = deep_basin_term(periods, v["depth_2_5"])
            if self._m9:
                f_basin = np.where(
                    is_slab, f_basin, deep_basin_term(periods, v["depth_2_5"], True)
                )
            ln_resp = ln_resp + f_basin
        # Frankel (2007): use the total standard deviation of Table 5
        ln_std = np.sqrt(
            c.sigma * c.sigma + np.where(is_slab, c.tauS * c.tauS, c.tau * c.tau)
        )
        return ln_resp, ln_std

    def _reference(self, periods, v) -> tuple:
        """Reference model of the extrapolation (``extrapolationRefs``)."""
        ln_resp_bc, ln_std_bc = bchydro_2012(periods, v, self._basin)
        is_interface = v["is_interface"]
        if not np.any(is_interface):
            return ln_resp_bc, np.broadcast_to(ln_std_bc, np.shape(ln_resp_bc))
        # AM_09_INTERFACE, AM_09_INTERFACE_BASIN, or AM_09_INTERFACE_BASIN_M9
        ln_resp_am, ln_std_am = AtkinsonMacias2009.calc_ln_base(
            periods, v, self._basin, self._m9, False
        )
        ln_resp = np.where(
            is_interface, 0.0 + ln_resp_bc * 0.5 + ln_resp_am * 0.5, ln_resp_bc
        )
        ln_std = np.where(
            is_interface, 0.0 + ln_std_bc * 0.5 + ln_std_am * 0.5, ln_std_bc
        )
        return ln_resp, ln_std

    def _spectrum(self, periods, v) -> tuple:
        periods = [float(p) for p in periods]
        extrapolated = [p for p in periods if p in self.EXTRAPOLATED]
        others = [p for p in periods if p not in self.EXTRAPOLATED]
        common = 5.0
        targets = sorted(set(others) | ({common} if extrapolated else set()))
        ln_resp_t, ln_std_t = interpolate_periods(
            targets, self.INTERPOLATED, lambda base: self._calc_ln_base(base, v)
        )
        cols = {p: i for i, p in enumerate(targets)}
        resp = {p: ln_resp_t[..., cols[p]] for p in others}
        std = {p: ln_std_t[..., cols[p]] for p in others}
        if extrapolated:
            ref_periods = [common] + extrapolated
            ln_resp_ref, ln_std_ref = self._reference(ref_periods, v)
            shape = np.broadcast_shapes(np.shape(ln_resp_ref), np.shape(ln_std_ref))
            ln_resp_ref = np.broadcast_to(ln_resp_ref, shape)
            ln_std_ref = np.broadcast_to(ln_std_ref, shape)
            # ExtrapolatedGmm
            resp_scale = ln_resp_t[..., cols[common]] / ln_resp_ref[..., 0]
            std_scale = ln_std_t[..., cols[common]] / ln_std_ref[..., 0]
            for i, p in enumerate(extrapolated, start=1):
                resp[p] = ln_resp_ref[..., i] * resp_scale
                std[p] = ln_std_ref[..., i] * std_scale
        shape = np.broadcast_shapes(
            *[np.shape(x) for x in list(resp.values()) + list(std.values())]
        )
        return (
            np.stack([np.broadcast_to(resp[p], shape) for p in periods], axis=-1),
            np.stack([np.broadcast_to(std[p], shape) for p in periods], axis=-1),
        )
