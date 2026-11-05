"""Coppersmith and Bommer (2014) ground motion model.

The model is part of the PNNL Site Wide Hazard study, :cite:t:`coppersmith14`.
"""

import numpy as np

from .. import model
from ..registry import register
from .abrahamson_gregor_addo_2016 import AbrahamsonGregorAddo2016

__author__ = "Albert Kottke"


@register(tectonic="subduction")
class CoppersmithBommer2014(AbrahamsonGregorAddo2016):
    """Coppersmith and Bommer (2014) ground motion model.

    The model was developed as part of the PNNL hazard study for the Hanford
    site, :cite:p:`coppersmith14`.

    This model was developed for subduction regions and is a modified version
    of the BCHydro model.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``v_s30``, and ``event_type``) can be a scalar or an array, and the arrays
    are broadcast against each other. For a scalar scenario, the response and
    standard deviation have one value per period, as in other models. For
    arrays of N scenarios, they have shape (N, periods), so, for example,
    ``pga`` has shape (N,).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    scale_atten : float, optional
        scale factor applied to the anelastic attenuation term. Default is 1.
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 14 periods). PGA is the
        response at 0.01 s. The model does not provide all of the
        "psa_ngawest2_21" periods. Computing only the needed intensity
        measures is much faster for large vectorized scenarios. If *None*
        (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([7.0, 8.0]), dist_rup=np.array([50.0, 100.0]),
    ...     v_s30=400.0, event_type=np.array(["intraslab", "interface"]))
    >>> pygmm.CoppersmithBommer2014(s, ims=["pga"]).pga.shape
    (2,)
    >>> pygmm.CoppersmithBommer2014(s).spec_accels.shape
    (2, 14)

    """

    NAME = "Coppersmith and Bommer (2014)"
    ABBREV = "CB14"

    # Load the coefficients for the model
    COEFF = model.load_data_file("coppersmith_bommer_2014.csv", 1)
    PERIODS = COEFF["period"]

    INDEX_PGA = 0
    INDEX_PGV = None
    INDICES_PSA = np.arange(14)

    # FIXME
    LIMITS = dict(
        mag=(3.0, 8.5),
        dist_jb=(0.0, 300.0),
        v_s30=(150.0, 1500.0),
    )

    PARAMS = [
        model.NumericParameter("mag", True, 3, 8.5),
        model.NumericParameter("dist_rup", False),
        model.NumericParameter("v_s30", True, 150.0, 1500.0),
        model.CategoricalParameter("event_type", True, ["interface", "intraslab"]),
    ]

    def __init__(self, scenario, scale_atten=1.0, ims=None):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        scale_atten : float, optional
            scale factor applied to the anelastic attenuation term.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed.
        """
        super().__init__(
            scenario, adjust_c1=None, adjust_c4=0, scale_atten=scale_atten, ims=ims
        )

    def _calc_ln_resp(self, pga_ref, rows=None):
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
        dist_rup = model.as_column(s.dist_rup)

        path_atten = (c.t_2 + c.t_3 * (mag - 7.8)) * np.log(
            dist_rup + c.c_4 * np.exp(c.t_9 * (mag - 6))
        )

        ln_resp = (
            c.t_1
            + c.t_4 * adjust_c1
            + path_atten
            + self._scale_atten * c.t_6 * dist_rup
            + self._calc_f_mag(c, adjust_c1, mag)
            + self._calc_f_site(c, pga_ref)
            + self._calc_f_faba(c, dist_rup)
        )

        return ln_resp

    def _calc_f_faba(self, c, dist, event_type=None):
        """Calculate the forearc/backarc scaling. Equation 4.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute
        dist : array_like
            rupture distance [km]
        event_type : array_like, optional
            not used by the model.

        Returns
        -------
        f_faba : class:`np.array`:
            forearc/backarc scaling

        """
        f_faba = c.t_16 * np.log(np.maximum(dist, 40) / 40)

        return f_faba
