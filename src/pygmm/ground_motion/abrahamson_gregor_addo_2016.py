"""Abrahamson, Gregor, and Addo (2016) ground motion model."""

from typing import Any

import numpy as np

from .. import model
from .._types import ArrayLike
from ..registry import register

__author__ = "Albert Kottke"


@register(tectonic="subduction")
class AbrahamsonGregorAddo2016(model.GroundMotionModel):
    r"""Abrahamson, Gregor, and Addo (2016) ground motion model.

    This model was developed for subduction regions and is commonly referred to
    as the BCHydro model.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``dist_hyp``, ``depth_hyp``, ``v_s30``, ``event_type``, and
    ``tectonic_region``) can be a scalar or an array, and the arrays are
    broadcast against each other. For a scalar scenario, the response and
    standard deviation have one value per period, as in other models. For
    arrays of N scenarios, they have shape (N, periods), so, for example,
    ``pga`` has shape (N,). Intraslab events use ``dist_hyp`` and
    ``depth_hyp``, and interface events use ``dist_rup``, so values that are
    not needed by any of the events can be *None*.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    adjust_c1 : float or array_like, optional
        adjustment to the magnitude scaling break point (:math:`\Delta C_1`)
        as a single value or one value per period. If *None* (default), -0.3
        is used for intraslab events and a period dependent value from 0.2 to
        -0.2 is used for interface events.
    adjust_c4 : float, optional
        adjustment to the finite-fault distance term (:math:`\Delta C_4`,
        km), which is added to :math:`C_4` (10 km). Default is 0.
    scale_atten : float, optional
        scale factor applied to the anelastic attenuation term. Default is 1.
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 22 periods). The model does
        not provide all of the "psa_ngawest2_21" periods. Computing only the
        needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([7.0, 8.5]), dist_rup=np.array([50.0, 100.0]),
    ...     dist_hyp=np.array([80.0, 120.0]), depth_hyp=60.0, v_s30=400.0,
    ...     event_type=np.array(["intraslab", "interface"]),
    ...     tectonic_region="forearc")
    >>> pygmm.AbrahamsonGregorAddo2016(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.AbrahamsonGregorAddo2016(s, ims=["pga", "psa_0p200", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p200', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

    """

    NAME = "Abrahamson, Gregor, & Addo (2016)"
    ABBREV = "AGA16"

    # Reference shear-wave velocity in m/sec
    V_REF = 1000.0

    # Load the coefficients for the model
    COEFF = model.load_data_file("abrahamson_gregor_addo_2016.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDICES_PSA = np.arange(1, 23)

    # FIXME
    LIMITS = dict(
        mag=(3.0, 8.5),
        dist_jb=(0.0, 300.0),
        v_s30=(150.0, 1500.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 3, 8.5),
        model.NumericParameter("dist_rup", False),
        model.NumericParameter("dist_hyp", False),
        model.NumericParameter("depth_hyp", False),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
        model.CategoricalParameter(
            "tectonic_region", False, ["forearc", "backarc", "unknown"], "unknown"
        ),
    ]

    def __init__(
        self, scenario, adjust_c1=None, adjust_c4=0, scale_atten=1.0, ims=None
    ):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        adjust_c1 : float or array_like, optional
            adjustment to the magnitude scaling break point. If *None* (default), the
            adjustment depends on the event type.
        adjust_c4 : float, optional
            adjustment added to the finite-fault distance term, C_4 [km].
        scale_atten : float, optional
            scale factor applied to the anelastic attenuation term.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed.
        """
        super().__init__(scenario, ims)

        n = len(self.COEFF)
        if adjust_c1 is None:
            # Default adjustments depend on the event type of each scenario, so
            # they are computed with the response. See _calc_adjust_c1.
            self._adjust_c1 = None
        else:
            if isinstance(adjust_c1, float):
                self._adjust_c1 = adjust_c1 * np.ones(n)
            else:
                self._adjust_c1 = np.asarray(adjust_c1)

        self._adjust_c4 = adjust_c4
        self._scale_atten = scale_atten

        # The nonlinear site term depends on PGA at the reference condition,
        # which only needs the PGA coefficients
        pga_ref = np.exp(self._calc_ln_resp(np.nan, [self.INDEX_PGA])[..., 0])
        self._ln_resp = self._calc_ln_resp(pga_ref)
        self._ln_std = self._calc_ln_std()

    @property
    def adjust_c1(self):
        r"""Adjustment to the magnitude scaling break point (:math:`\Delta C_1`).

        The default adjustment has one value per computed period, with shape
        (N, periods) for N scenarios with different event types.
        """
        if self._adjust_c1 is None:
            return self._calc_adjust_c1()
        return self._adjust_c1

    @property
    def adjust_c4(self):
        r"""Adjustment to the finite-fault distance term (:math:`\Delta C_4`, km)."""
        return self._adjust_c4

    @property
    def scale_atten(self):
        """Scale factor applied to the attenuation."""
        return self._scale_atten

    def _select_rows(self, values: Any, rows=None) -> Any:
        """Select coefficient rows (periods) from per-period values.

        If `rows` is *None*, the rows for the requested intensity measures are
        selected.
        """
        if rows is None:
            return self._coeff_rows(values)
        if not isinstance(values, np.ndarray):
            values = np.asarray(values)
        # Indexing keeps the array type, so coefficient record arrays keep
        # attribute access (e.g., ``c.t_1``)
        return values[rows]

    def _calc_adjust_c1(self, rows=None) -> np.ndarray:
        """Calculate the adjustment to the magnitude scaling break point.

        Parameters
        ----------
        rows : array_like, optional
            coefficient rows (periods) to compute. If *None*, the rows for the
            requested intensity measures are used.

        Returns
        -------
        adjust_c1 : class:`np.array`:
            adjustment for each scenario and period

        """
        if self._adjust_c1 is not None:
            # Specified adjustment with one value per period
            return self._select_rows(
                np.broadcast_to(self._adjust_c1, len(self.COEFF)), rows
            )

        # Default adjustments
        period = self._select_rows(self.COEFF.period, rows)
        adjust_interface = np.interp(
            np.log(np.maximum(0.01, period)),
            np.log([0.3, 0.5, 1, 2, 3]),
            [0.2, 0.1, 0, -0.1, -0.2],
            left=0.2,
            right=-0.2,
        )
        is_slab = model.as_column(model.equals(self._scenario.event_type, "intraslab"))
        return np.where(is_slab, -0.3, adjust_interface)

    def _calc_ln_resp(self, pga_ref: ArrayLike, rows=None) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Parameters
        ----------
        pga_ref : float or array_like
            peak ground acceleration [g] at the reference
            condition. If :class:`np.nan`, then the response at the reference
            condition is computed.
        rows : array_like, optional
            coefficient rows (periods) to compute. If *None*, the rows for the
            requested intensity measures are used.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response

        """
        s = self._scenario
        c = self._select_rows(self.COEFF, rows)
        adjust_c1 = self._calc_adjust_c1(rows)
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        mag = model.as_column(s.mag)
        event_type = model.as_column(s.event_type)
        is_slab = model.equals(event_type, "intraslab")

        f_event = np.where(is_slab, 1, 0)
        dist = self._calc_dist(is_slab)

        path_atten = (c.t_2 + c.t_14 * f_event + c.t_3 * (mag - 7.8)) * np.log(
            dist + (c.c_4 + self._adjust_c4) * np.exp(c.t_9 * (mag - 6))
        )

        ln_resp = (
            c.t_1
            + c.t_4 * adjust_c1
            + path_atten
            + self._scale_atten * c.t_6 * dist
            + c.t_10 * f_event
            + self._calc_f_mag(c, adjust_c1, mag)
            + self._calc_f_site(c, pga_ref)
        )

        # Forearc/backarc scaling is only applied to backarc sites
        is_backarc = model.as_column(model.equals(s.tectonic_region, "backarc"))
        if np.ndim(s.tectonic_region) > 0 or is_backarc.all():
            ln_resp = ln_resp + np.where(
                is_backarc, self._calc_f_faba(c, dist, event_type), 0
            )

        # Depth scaling is only applied to intraslab events
        if s.depth_hyp is not None:
            ln_resp = ln_resp + np.where(
                is_slab, self._calc_f_depth(c, model.as_column(s.depth_hyp)), 0
            )
        elif np.any(is_slab):
            raise ValueError("depth_hyp is required for intraslab events")

        return ln_resp

    def _calc_dist(self, is_slab: np.ndarray) -> np.ndarray:
        """Distance used by the model.

        Intraslab events use the hypocentral distance and interface events use
        the rupture distance.

        Parameters
        ----------
        is_slab : array_like
            if the event is an intraslab event

        Returns
        -------
        dist : class:`np.array`:
            distance [km]

        """
        s = self._scenario
        if s.dist_hyp is None:
            if np.any(is_slab):
                raise ValueError("dist_hyp is required for intraslab events")
            return model.as_column(s.dist_rup)
        if s.dist_rup is None:
            if not np.all(is_slab):
                raise ValueError("dist_rup is required for interface events")
            return model.as_column(s.dist_hyp)
        return np.where(
            is_slab, model.as_column(s.dist_hyp), model.as_column(s.dist_rup)
        )

    def _calc_f_mag(
        self, c, adjust_c1: float | np.ndarray, mag: float | np.ndarray
    ) -> np.ndarray:
        """Calculate the magnitude scaling. Equation 2.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute
        adjust_c1 : array_like
            adjustment to the magnitude scaling break point
        mag : array_like
            moment magnitude

        Returns
        -------
        f_mag : class:`np.array`:
            magnitude scaling

        """
        c_1 = 7.8
        coeff = np.where(mag <= (c_1 + adjust_c1), c.t_4, c.t_5)
        f_mag = coeff * (mag - (c_1 + adjust_c1)) + c.t_13 * (10 - mag) ** 2
        return f_mag

    def _calc_f_depth(self, c, depth_hyp: ArrayLike) -> np.ndarray:
        """Calculate the depth scaling. Equation 3.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute
        depth_hyp : array_like
            depth to hypocenter [km]

        Returns
        -------
        f_depth : class:`np.array`:
            depth scaling

        """
        f_depth = c.t_11 * (np.minimum(depth_hyp, 120) - 60)
        return f_depth

    def _calc_f_faba(self, c, dist: ArrayLike, event_type: ArrayLike) -> np.ndarray:
        """Calculate the forearc/backarc scaling. Equation 4.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute
        dist : array_like
            distance used by the model [km]
        event_type : array_like
            type of subduction event, either: 'interface' or 'intraslab'

        Returns
        -------
        f_faba : class:`np.array`:
            forearc/backarc scaling. *NaN* for events that are not interface or
            intraslab events.

        """
        is_slab = model.equals(event_type, "intraslab")
        is_interface = model.equals(event_type, "interface")
        if np.ndim(is_slab) == 1 and not (is_slab.any() or is_interface.any()):
            # Scalar scenario with an invalid event type
            raise NotImplementedError

        f_faba = np.select(
            [is_slab, is_interface],
            [
                c.t_7 + c.t_8 * np.log(np.maximum(dist, 85) / 40),
                c.t_15 + c.t_16 * np.log(np.maximum(dist, 100) / 40),
            ],
            default=np.nan,
        )
        return f_faba

    def _calc_f_site(self, c, pga_ref: ArrayLike) -> np.ndarray:
        """Calculate the site amplification. Equation 5.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute
        pga_ref : float or array_like
            reference PGA (V_s30 = 1000 m/sec). If 'np.nan', then the ground
            motion at the reference condition is computed.

        Returns
        -------
        f_site : class:`np.array`:
            site scaling

        """
        if np.ndim(pga_ref) == 0 and np.isnan(pga_ref):
            # Reference condition
            v_s30: float | np.ndarray = 1000.0
        else:
            v_s30 = model.as_column(self._scenario.v_s30)
            pga_ref = model.as_column(pga_ref)
        vs_ratio = np.minimum(v_s30, 1000) / c.v_lin

        f_site = np.where(
            v_s30 < c.v_lin,
            c.t_12 * np.log(vs_ratio)
            - c.b * np.log(pga_ref + c.c)
            + c.b * np.log(pga_ref + c.c * vs_ratio**c.n),
            (c.t_12 + c.b * c.n) * np.log(vs_ratio),
        )
        return f_site

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        c = self._coeff_rows(self.COEFF)

        ln_std = np.sqrt(c.phi**2 + c.tau**2)
        if np.ndim(self._resp_vector) > 1:
            # The standard deviation does not depend on the scenario, so it is
            # repeated for each scenario
            ln_std = np.broadcast_to(ln_std, np.shape(self._resp_vector))
        return ln_std
