r"""Atkinson and Macias (2009) subduction interface model as in nshmp-lib.

The model follows the ``AtkinsonMacias_2009`` class of the USGS nshmp-lib
library (commit 44728a7d), which is used for the current U.S. Geological
Survey national seismic hazard models (NSHMs). The module also provides the
nshmp-lib helpers that are shared with :mod:`pygmm.ground_motion.zhao_et_al_2006`: the
Campbell and Bozorgnia (2014) deep basin term used by the USGS basin variants
(:func:`deep_basin_term`) and the Abrahamson and Bhasin (2020) conditional PGV
model that nshmp-lib uses for PGV (:func:`conditional_pgv`).

References
----------
Atkinson, G. M., and Macias, D. M. (2009). Predicted ground motions for great
interface earthquakes in the Cascadia subduction zone. Bulletin of the
Seismological Society of America, 99(3), 1552-1578.
https://doi.org/10.1785/0120080147

USGS nshmp-lib, ``gov.usgs.earthquake.nshmp.gmm.AtkinsonMacias_2009``,
``BooreAtkinson_2008``, ``CampbellBozorgnia_2014``, ``GmmUtils``,
``InterpolatedGmm``, and ``UsgsPgvSupport``, commit 44728a7d
(https://code.usgs.gov/ghsc/nshmp/nshmp-lib).
"""

from typing import Any

import numpy as np

from .. import model
from .._types import ArrayLike
from ..registry import register
from .campbell_bozorgnia_2014 import CampbellBozorgnia2014

#: Natural log of 10 (``GmmUtils.BASE_10_TO_E``)
LN_10 = np.log(10.0)
#: Natural log of the acceleration of gravity in cm/sec² (``GmmUtils.LN_G_CM_TO_M``)
LN_G = np.log(980.0)

#: Spectral periods (s) provided by the nshmp-lib subduction models
PERIODS_SA = model.GroundMotionModel.PERIODS_NGAWEST2_21

#: Campbell and Bozorgnia (2014) coefficients sorted by period (for the basin term),
#: which are the same as the nshmp-lib coefficients
_COEFF_CB14 = np.sort(CampbellBozorgnia2014.COEFF, order="period").view(np.recarray)

# Abrahamson and Bhasin (2020) conditional PGV model (UsgsPgvSupport)
_AB20_B1 = -4.09
_AB20_B2 = 0.66
_AB20_A = (5.39, 0.799, 0.654, 0.479, -0.062, -0.359, -0.134, 0.023)
_AB20_SIGMA = 0.33


def coeff_rows(coeff: np.recarray, periods: ArrayLike) -> Any:
    """Rows of a coefficient table (with a ``period`` column) for periods."""
    periods = np.asarray(periods, dtype=float)
    idx = np.searchsorted(coeff.period, periods)
    idx = np.clip(idx, 0, len(coeff) - 1)
    if not np.all(np.isclose(coeff.period[idx], periods, rtol=1e-9, atol=0)):
        raise ValueError(f"periods {periods} are not all in the coefficient table")
    return coeff[idx]


def deep_basin_term(
    periods: ArrayLike, depth_2_5: ArrayLike, m9: bool = False
) -> np.ndarray:
    r"""USGS deep basin term of the subduction models in nshmp-lib.

    nshmp-lib adds the depth-tapered deep basin term of Campbell and
    Bozorgnia (2014) (``CampbellBozorgnia_2014.deepBasinScaling``) to the
    ``_BASIN`` variants of the subduction models (Petersen et al., 2020). The
    term is zero for periods of 0.5 s and shorter (and PGA), and for
    :math:`Z_{2.5}` of 1 km or less or NaN. Otherwise, the CB14 basin term
    (which is zero from 1 to 3 km and :math:`c_{16} k_3 e^{-0.75} (1 -
    e^{-0.25 (Z_{2.5} - 3)})` deeper than 3 km) is multiplied by the taper
    :math:`(\min(Z_{2.5}, 3) - 1) / 2` (which is 1 where the term is not
    zero) and, at 0.75 s, by 0.585. With `m9`, the Seattle M9 adjustment
    replaces the term by :math:`\ln 2` for :math:`Z_{2.5}` greater than 6 km
    and periods longer than 1.9 s.

    Parameters
    ----------
    periods : array_like
        spectral periods [s], with 0 for PGA
    depth_2_5 : float or array_like
        depth to the 2.5 km/sec shear-wave velocity horizon [km], with a
        trailing axis for the periods
    m9 : bool, optional
        apply the Seattle M9 adjustment (default *False*)

    Returns
    -------
    f_basin : :class:`np.ndarray`
        basin term in natural log units

    Examples
    --------
    >>> from pygmm.ground_motion.atkinson_macias_2009 import deep_basin_term
    >>> deep_basin_term([0.0, 0.5, 1.0, 3.0], 2.0).tolist()
    [0.0, 0.0, 0.0, 0.0]
    >>> deep_basin_term([1.0, 3.0], 7.0, m9=True).round(4).tolist()
    [0.4441, 0.6931]
    """
    periods = np.asarray(periods, dtype=float)
    c = coeff_rows(_COEFF_CB14, periods)
    depth_2_5 = np.asarray(depth_2_5, dtype=float)
    with np.errstate(invalid="ignore"):
        scale = (np.clip(depth_2_5, 1.0, 3.0) - 1.0) / (3.0 - 1.0)
        # CampbellBozorgnia_2014.calcBasinTerm for Z2.5 > 1 km
        term = np.where(
            depth_2_5 > 3.0,
            c.c_16 * c.k_3 * np.exp(-0.75) * (1.0 - np.exp(-0.25 * (depth_2_5 - 3.0))),
            0.0,
        )
        f_basin = term * scale
        f_basin = np.where(periods == 0.75, f_basin * 0.585, f_basin)
        # GmmUtils.checkBasin (NaN depths compare False)
        f_basin = np.where((periods > 0.5) & (depth_2_5 > 1.0), f_basin, 0.0)
        if m9:
            f_basin = np.where(
                (depth_2_5 > 6.0) & (periods > 1.9), np.log(2.0), f_basin
            )
    return f_basin


def interpolate_periods(periods, interpolated: dict, calc) -> tuple:
    """Ground motions at periods, interpolating periods without coefficients.

    nshmp-lib ``InterpolatedGmm``: the natural log of the median and the
    standard deviation are interpolated linearly in period between the
    bounding periods.

    Parameters
    ----------
    periods : sequence of float
        target periods [s]
    interpolated : dict
        bounding periods ``(lo, hi)`` of the interpolated periods
    calc : callable
        ``calc(base_periods)`` returns the natural log of the median and the
        standard deviation at periods with coefficients, with a trailing axis
        for the periods

    Returns
    -------
    ln_resp : :class:`np.ndarray`
        natural log of the response at `periods`
    ln_std : :class:`np.ndarray`
        logarithmic standard deviation at `periods`
    """
    periods = [float(p) for p in periods]
    needed = set()
    for p in periods:
        needed.update(interpolated.get(p, (p,)))
    base = sorted(needed)
    ln_resp, ln_std = calc(np.array(base))
    shape = np.broadcast_shapes(np.shape(ln_resp), np.shape(ln_std))
    ln_resp = np.broadcast_to(ln_resp, shape)
    ln_std = np.broadcast_to(ln_std, shape)
    col = {p: i for i, p in enumerate(base)}
    resp, std = [], []
    for p in periods:
        if p in interpolated:
            lo, hi = interpolated[p]
            i, j = col[lo], col[hi]
            # Interpolator.findY
            resp.append(
                ln_resp[..., i]
                + (p - lo) * (ln_resp[..., j] - ln_resp[..., i]) / (hi - lo)
            )
            std.append(
                ln_std[..., i]
                + (p - lo) * (ln_std[..., j] - ln_std[..., i]) / (hi - lo)
            )
        else:
            resp.append(ln_resp[..., col[p]])
            std.append(ln_std[..., col[p]])
    return np.stack(resp, axis=-1), np.stack(std, axis=-1)


def conditional_pgv(
    periods: ArrayLike,
    ln_resp: ArrayLike,
    ln_std: ArrayLike,
    mag: ArrayLike,
    dist_rup: ArrayLike,
    v_s30: ArrayLike,
) -> tuple:
    r"""Abrahamson and Bhasin (2020) conditional PGV as in nshmp-lib.

    nshmp-lib (``UsgsPgvSupport.calcAB20Pgv``) computes the PGV of models
    without PGV coefficients from the spectral acceleration at the period
    :math:`T_{PGV} = \exp(-4.09 + 0.66 M)`, which is interpolated linearly in
    :math:`\ln T` from the natural log of the median and the standard
    deviation at the model's spectral periods (extrapolated below the
    shortest period, and NaN for :math:`T_{PGV}` of the longest period or
    longer), with the horizontal coefficients of Table 3.2 and a
    conditional standard deviation of 0.33.

    Parameters
    ----------
    periods : array_like
        spectral periods [s] of the model, in increasing order
    ln_resp : array_like
        natural log of the spectral accelerations [g], with the periods along
        the last axis
    ln_std : array_like
        logarithmic standard deviations of the spectral accelerations
    mag : float or array_like
        moment magnitude, with a trailing axis of length one
    dist_rup : float or array_like
        rupture distance [km], with a trailing axis of length one
    v_s30 : float or array_like
        time-averaged shear-wave velocity in the top 30 m [m/sec], with a
        trailing axis of length one

    Returns
    -------
    ln_pgv : :class:`np.ndarray`
        natural log of the PGV [cm/sec], with a trailing axis of length one
    ln_std_pgv : :class:`np.ndarray`
        logarithmic standard deviation of the PGV
    """
    a1, a2, a3, a4, a5, a6, a7, a8 = _AB20_A
    periods = np.asarray(periods, dtype=float)
    mag = np.asarray(mag, dtype=float)
    dist_rup = np.asarray(dist_rup, dtype=float)
    v_s30 = np.asarray(v_s30, dtype=float)
    shape = np.broadcast_shapes(np.shape(ln_resp), np.shape(ln_std), mag.shape)
    shape = shape[:-1] + (len(periods),)
    ln_resp = np.broadcast_to(ln_resp, shape)
    ln_std = np.broadcast_to(ln_std, shape)

    period_pgv = np.exp(_AB20_B1 + _AB20_B2 * mag)
    period_pgv = np.broadcast_to(period_pgv, shape[:-1] + (1,))
    # UsgsPgvSupport.calcSaGroundMotion: the lower period is the longest
    # period that is not longer than T_PGV (or the shortest period), and the
    # upper period is the next period (or the longest period)
    n = len(periods)
    k = np.searchsorted(periods, period_pgv, side="right")
    lo = np.maximum(k - 1, 0)
    hi = np.minimum(np.maximum(k, 1), n - 1)
    ln_periods = np.log(periods)
    x1 = ln_periods[lo]
    x2 = ln_periods[hi]
    x = np.log(period_pgv)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu_lo = np.take_along_axis(ln_resp, lo, axis=-1)
        mu_hi = np.take_along_axis(ln_resp, hi, axis=-1)
        sd_lo = np.take_along_axis(ln_std, lo, axis=-1)
        sd_hi = np.take_along_axis(ln_std, hi, axis=-1)
        mu = mu_lo + (x - x1) * (mu_hi - mu_lo) / (x2 - x1)
        sd = sd_lo + (x - x1) * (sd_hi - sd_lo) / (x2 - x1)

    # Equation 3.8
    f1 = np.where(
        mag < 5.0, a2, np.where(mag <= 7.5, a2 + (a3 - a2) * (mag - 5.0) / 2.5, a3)
    )
    # Equation 3.7
    ln_pgv = (
        a1
        + f1 * mu
        + a4 * (mag - 6.0)
        + a5 * (8.5 - mag) ** 2
        + a6 * np.log(dist_rup + 5.0 * np.exp(0.4 * (mag - 6.0)))
        + (a7 + a8 * (mag - 5.0)) * np.log(v_s30 / 425)
    )
    ln_std_pgv = np.sqrt(f1 * f1 * sd * sd + _AB20_SIGMA * _AB20_SIGMA)
    return ln_pgv, ln_std_pgv


class _NshmpSubductionBase(model.GroundMotionModel):
    """Shared code of the nshmp-lib subduction models with conditional PGV."""

    # Period (s), with 0 for PGA and -1 for PGV
    PERIODS = np.r_[-1.0, 0.0, PERIODS_SA]
    INDEX_PGV = 0
    INDEX_PGA = 1
    INDICES_PSA = np.arange(2, 23)
    # PGV is in cm/sec
    PGV_SCALE = 1.0

    #: Periods (s) that nshmp-lib interpolates, with the bounding periods
    INTERPOLATED: dict[float, tuple[float, float]] = {}

    def _scenario_values(self) -> dict:
        raise NotImplementedError

    def _spectrum(self, periods, v) -> tuple:
        """Ground motions at PGA and spectral periods."""
        raise NotImplementedError

    def _calc(self) -> tuple:
        """Natural log of the response and the standard deviation."""
        v = self._scenario_values()
        rows = self._coeff_rows(np.arange(len(self.PERIODS)))
        periods = self.PERIODS[rows]
        has_pgv = np.any(periods < 0)
        needed = set(periods[periods >= 0].tolist())
        if has_pgv:
            # The conditional PGV model needs all of the spectral periods
            needed.update(PERIODS_SA.tolist())
        targets = sorted(needed)
        ln_resp_t, ln_std_t = self._spectrum(targets, v)
        cols = {p: i for i, p in enumerate(targets)}
        resp, std = [], []
        for p in periods:
            if p < 0:
                sa = [cols[float(t)] for t in PERIODS_SA]
                ln_pgv, ln_std_pgv = conditional_pgv(
                    PERIODS_SA,
                    ln_resp_t[..., sa],
                    ln_std_t[..., sa],
                    v["mag"],
                    v["dist_rup"],
                    v["v_s30"],
                )
                resp.append(ln_pgv[..., 0])
                std.append(ln_std_pgv[..., 0])
            else:
                resp.append(ln_resp_t[..., cols[float(p)]])
                std.append(ln_std_t[..., cols[float(p)]])
        ln_resp = np.stack(resp, axis=-1)
        ln_std = np.stack(std, axis=-1)
        valid = v["valid"]
        ln_resp = np.where(valid, ln_resp, np.nan)
        ln_std = np.where(valid, ln_std, np.nan)
        # Some terms do not depend on all of the scenario values (e.g.,
        # depth_2_5 without the basin term)
        shape = np.broadcast_shapes(
            np.shape(ln_resp),
            np.shape(ln_std),
            *[np.shape(value)[:-1] + (1,) for value in v.values() if value is not None],
        )
        return (
            np.array(np.broadcast_to(ln_resp, shape), dtype=float),
            np.array(np.broadcast_to(ln_std, shape), dtype=float),
        )

    def _depth_2_5(self) -> np.ndarray:
        """Depth to the 2.5 km/sec horizon, with NaN for no basin term."""
        if self._scenario["depth_2_5"] is None:
            return np.array([np.nan])
        return model.as_column(self._scenario["depth_2_5"]).astype(float)


@register(tectonic="subduction")
class AtkinsonMacias2009(_NshmpSubductionBase):
    r"""Atkinson and Macias (2009) subduction interface model (nshmp-lib).

    Stochastic finite-fault model for great Cascadia subduction interface
    earthquakes, as implemented in the ``AtkinsonMacias_2009`` class of USGS
    nshmp-lib (commit 44728a7d), with the coefficients of ``AM09.csv``.
    nshmp-lib uses the model with the 2014 NSHM implementation:

    * The model is for a reference site (NEHRP B/C). The site term is the
      Boore and Atkinson (2008) site amplification (``BooreAtkinson_2008``,
      ``BA08.csv``) at the same period, with the PGA of the model at the
      reference site as the reference PGA of the nonlinear term.
    * nshmp-lib passes the natural log of the reference PGA, instead of the
      PGA [g], to the BA08 site term (``lnPgaRock``), so that the nonlinear
      term is always that of a weak-motion reference PGA (0.03 g or less),
      :math:`b_{nl} \ln(0.06 / 0.1)`. ``site_fix=True`` (the ``_SITE_FIX``
      variants, added in nshmp-lib in 2025) uses the PGA.
    * The near-source term uses :math:`h = M^2 - 3.1 M - 14.55` and the
      rupture distance (there is no fixed hypocentral depth).
    * PGA coefficients are used for 0.01 s. The 0.02, 0.03, 0.075, 0.15,
      0.25, and 1.5 s results (including the basin term) are interpolated
      (linearly in period) between 0.01 and 0.05 s, 0.05 and 0.1 s, 0.1 and
      0.2 s, 0.2 and 0.3 s, and 1.0 and 2.0 s (``InterpolatedGmm``). 7.5 s is
      used for the 0.13 Hz values (7.7 s in the NSHM Fortran code).
    * PGV [cm/sec] is computed with the Abrahamson and Bhasin (2020)
      conditional model from the response spectrum of the model
      (:func:`conditional_pgv`).
    * The standard deviation is the total standard deviation of the model
      (``sig``, converted to natural log units).
    * The scenario values are not clipped to the model limits.

    The nshmp-lib ``Gmm`` ids correspond to the options as follows (see
    :attr:`GMM_IDS`):

        ================================= =====================================
        nshmp-lib Gmm id                  options
        ================================= =====================================
        AM_09_INTERFACE
        AM_09_INTERFACE_BASIN             basin=True
        AM_09_INTERFACE_BASIN_M9          basin=True, m9=True
        AM_09_INTERFACE_BASIN_SITE_FIX    basin=True, site_fix=True
        AM_09_INTERFACE_BASIN_M9_SITE_FIX basin=True, m9=True, site_fix=True
        ================================= =====================================

    ``AM_09_INTERFACE_BASIN`` is used, with a weight of 0.125, by the
    Cascadia interface logic tree of the current conterminous U.S. NSHM
    (nshm-conus 6.2.0).

    The model is for interface events. The ``event_type`` scenario value is
    optional ("interface" by default); results for other event types (e.g.,
    "intraslab") are *NaN*.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``v_s30``, ``depth_2_5``, and ``event_type``) can be a scalar or an
    array, and the arrays are broadcast against each other. For a scalar
    scenario, the response and standard deviation have one value per period,
    as in other models. For arrays of N scenarios, they have shape (N,
    periods), so, for example, ``pga`` has shape (N,). ``depth_2_5`` [km] is
    only used with ``basin=True``, and can be *None* (or *NaN* for some
    scenarios) for no basin term.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    basin : bool, optional
        if *True*, add the USGS deep basin term (see
        :func:`deep_basin_term`; the ``_BASIN`` variants). Default is
        *False*.
    m9 : bool, optional
        if *True*, apply the Seattle M9 adjustment, which replaces the basin
        term by :math:`\ln 2` for :math:`Z_{2.5}` greater than 6 km and
        periods of 2 s and longer (and so, by the interpolation, affects 1.5
        s; the ``_M9`` variants). It requires ``basin=True``. Default is
        *False*.
    site_fix : bool, optional
        if *True*, pass the reference PGA [g], instead of its natural log, to
        the BA08 nonlinear site term (the ``_SITE_FIX`` variants). Default is
        *False*, as in ``AM_09_INTERFACE`` and ``AM_09_INTERFACE_BASIN``.
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), "psa_ngawest2_21" (the 21 NGA-West2 comparison
        periods, which are all of the periods), and/or "psa_all". If *None*
        (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> from pygmm.model import Scenario
    >>> from pygmm.ground_motion.atkinson_macias_2009 import AtkinsonMacias2009
    >>> s = Scenario(
    ...     mag=np.array([8.0, 9.0]), dist_rup=np.array([50.0, 100.0]),
    ...     v_s30=400.0, depth_2_5=np.array([np.nan, 7.0]))
    >>> m = AtkinsonMacias2009(s, **AtkinsonMacias2009.GMM_IDS["AM_09_INTERFACE_BASIN"])
    >>> m.pga.shape, m.spec_accels.shape
    ((2,), (2, 21))
    >>> AtkinsonMacias2009(s, ims=["pga"]).pga.shape
    (2,)
    """

    NAME = "Atkinson and Macias (2009)"
    ABBREV = "AM09"

    COEFF = model.load_data_file("atkinson_macias_2009-nshmp.csv", 1)
    #: Boore and Atkinson (2008) coefficients of nshmp-lib (site term)
    COEFF_BA08 = model.load_data_file("boore_atkinson_2008-nshmp.csv", 1)

    INTERPOLATED = {
        0.02: (0.01, 0.05),
        0.03: (0.01, 0.05),
        0.075: (0.05, 0.1),
        0.15: (0.1, 0.2),
        0.25: (0.2, 0.3),
        1.5: (1.0, 2.0),
    }

    #: Model options of the nshmp-lib ``Gmm`` ids
    GMM_IDS = {
        "AM_09_INTERFACE": dict(),
        "AM_09_INTERFACE_BASIN": dict(basin=True),
        "AM_09_INTERFACE_BASIN_M9": dict(basin=True, m9=True),
        "AM_09_INTERFACE_BASIN_SITE_FIX": dict(basin=True, site_fix=True),
        "AM_09_INTERFACE_BASIN_M9_SITE_FIX": dict(basin=True, m9=True, site_fix=True),
    }

    LIMITS = dict(mag=(5.0, 9.5), dist_rup=(0.0, 1000.0), v_s30=(150.0, 1500.0))

    PARAMS = [
        model.NumericParameter("mag", True, 5.0, 9.5),
        model.NumericParameter("dist_rup", True, 0.0, 1000.0),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.NumericParameter("depth_2_5", False),
        model.CategoricalParameter(
            "event_type", False, ["interface", "intraslab"], "interface"
        ),
    ]

    # BA08 site term constants
    _BA08_PGA_LO = 0.06
    _BA08_A1 = 0.03
    _BA08_A2 = 0.09
    _BA08_V1 = 180.0
    _BA08_V2 = 300.0
    _BA08_V_REF = 760.0

    def __init__(
        self,
        scenario: model.Scenario,
        basin: bool = False,
        m9: bool = False,
        site_fix: bool = False,
        ims=None,
    ):
        """Initialize the model."""
        if m9 and not basin:
            raise ValueError(
                "m9=True requires basin=True (nshmp-lib AM_09_INTERFACE_BASIN_M9)"
            )
        super().__init__(scenario, ims)
        self._basin = bool(basin)
        self._m9 = bool(m9)
        self._site_fix = bool(site_fix)
        self._ln_resp, self._ln_std = self._calc()

    @property
    def basin(self) -> bool:
        """If the USGS deep basin term is added."""
        return self._basin

    @property
    def m9(self) -> bool:
        """If the Seattle M9 adjustment is applied."""
        return self._m9

    @property
    def site_fix(self) -> bool:
        """If the reference PGA [g] is used by the nonlinear site term."""
        return self._site_fix

    def _scenario_values(self) -> dict:
        s = self._scenario
        return dict(
            mag=model.as_column(s.mag).astype(float),
            dist_rup=model.as_column(s.dist_rup).astype(float),
            v_s30=model.as_column(s.v_s30).astype(float),
            depth_2_5=self._depth_2_5(),
            valid=model.equals(model.as_column(s.event_type), "interface"),
        )

    def _spectrum(self, periods, v) -> tuple:
        return interpolate_periods(
            periods,
            self.INTERPOLATED,
            lambda base: self.calc_ln_base(
                base, v, self._basin, self._m9, self._site_fix
            ),
        )

    @classmethod
    def _ln_mean_ref(cls, c, mag, dist_rup) -> np.ndarray:
        """Natural log of the reference site response [g]."""
        h = (mag * mag) - (3.1 * mag) - 14.55
        d_m = mag - 8.0
        r = np.sqrt(dist_rup * dist_rup + h * h)
        gnd = c.c0 + (c.c3 * d_m) + (c.c4 * d_m * d_m)
        gnd = gnd + (c.c1 * np.log10(r) + c.c2 * r)
        return gnd * LN_10 - LN_G

    @classmethod
    def ba08_site_term(cls, periods, pga_ref, v_s30) -> np.ndarray:
        """Boore and Atkinson (2008) site term (``BooreAtkinson_2008.siteAmp``).

        Parameters
        ----------
        periods : array_like
            spectral periods [s], with 0 for PGA
        pga_ref : float or array_like
            reference PGA [g] of the nonlinear term
        v_s30 : float or array_like
            time-averaged shear-wave velocity in the top 30 m [m/sec]

        Returns
        -------
        f_site : :class:`np.ndarray`
            site term in natural log units
        """
        c = coeff_rows(cls.COEFF_BA08, periods)
        v_ref = cls._BA08_V_REF
        v_1, v_2 = cls._BA08_V1, cls._BA08_V2
        a_1, a_2 = cls._BA08_A1, cls._BA08_A2
        pga_lo = cls._BA08_PGA_LO
        f_lin = c.b_lin * np.log(v_s30 / v_ref)
        b_nl = np.where(
            v_s30 < v_ref,
            np.where(
                v_s30 > v_2,
                c.b2 * np.log(v_s30 / v_ref) / np.log(v_2 / v_ref),
                np.where(
                    v_s30 > v_1,
                    (c.b1 - c.b2) * np.log(v_s30 / v_2) / np.log(v_1 / v_2) + c.b2,
                    c.b1,
                ),
            ),
            0.0,
        )
        d_x = np.log(a_2 / a_1)
        d_y = b_nl * np.log(a_2 / pga_lo)
        c_ = (3.0 * d_y - b_nl * d_x) / (d_x * d_x)
        d = -(2.0 * d_y - b_nl * d_x) / (d_x * d_x * d_x)
        with np.errstate(divide="ignore", invalid="ignore"):
            p = np.log(pga_ref / a_1)
            ln_pga_ref = np.log(pga_ref / 0.1)
        f_nl = np.where(
            pga_ref <= a_1,
            b_nl * np.log(pga_lo / 0.1),
            np.where(
                pga_ref <= a_2,
                b_nl * np.log(pga_lo / 0.1) + (c_ * p * p) + (d * p * p * p),
                b_nl * ln_pga_ref,
            ),
        )
        return f_lin + f_nl

    @classmethod
    def calc_ln_base(cls, periods, v, basin=False, m9=False, site_fix=False):
        """Ground motions at the periods with coefficients.

        Parameters
        ----------
        periods : array_like
            periods [s] with coefficients, with 0 for PGA
        v : dict
            ``mag``, ``dist_rup``, ``v_s30``, and ``depth_2_5`` with a
            trailing axis for the periods
        basin, m9, site_fix : bool
            model options

        Returns
        -------
        ln_resp : :class:`np.ndarray`
            natural log of the response [g]
        ln_std : :class:`np.ndarray`
            logarithmic standard deviation
        """
        periods = np.asarray(periods, dtype=float)
        c = coeff_rows(cls.COEFF, periods)
        c_pga = coeff_rows(cls.COEFF, [0.0])
        ln_resp_ref = cls._ln_mean_ref(c, v["mag"], v["dist_rup"])
        ln_pga_ref = cls._ln_mean_ref(c_pga, v["mag"], v["dist_rup"])
        pga_ref = np.exp(ln_pga_ref) if site_fix else ln_pga_ref
        ln_resp = ln_resp_ref + cls.ba08_site_term(periods, pga_ref, v["v_s30"])
        if basin:
            ln_resp = ln_resp + deep_basin_term(periods, v["depth_2_5"], m9)
        ln_std = c.sig * LN_10
        return ln_resp, ln_std
