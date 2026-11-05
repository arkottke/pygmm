"""Akkar, Sandikkaya, & Bommer (2014, :cite:p:`akkar14`) model."""

import collections

import numpy as np

from .. import model
from ..registry import register

__author__ = "Albert Kottke"


@register()
class AkkarSandikkayaBommer2014(model.GroundMotionModel):
    """Akkar, Sandikkaya, & Bommer (2014, :cite:p:`akkar14`) model.

    The model is specified for three different distance metrics. However,
    the implementation uses only one distance metric. They are used in
    the following order:

        1. `dist_jb`

        2. `dist_hyp`

        3. `dist_epi`

    This order was selected based on evaluation of the total standard
    deviation. To compute the response for differing metrics, call the
    model multiple times with different keywords.

    The model is vectorized. Each scenario value (``mag``, ``v_s30``,
    ``mechanism``, and the distance metric ``dist_jb``, ``dist_hyp``, or
    ``dist_epi``) can be a scalar or an array, and the arrays are broadcast
    against each other. The distance metric is selected by which scenario
    value is provided, so it is the same for all scenarios. For a scalar
    scenario, the response and standard deviation have one value per period,
    as in other models. For arrays of N scenarios, they have shape
    (N, periods), so, for example, ``pga`` has shape (N,).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 62 periods, from 0.01 to
        4 s, so "psa_ngawest2_21" is not available). Computing only
        the needed intensity measures is much faster for large vectorized
        scenarios. If *None* (default), all intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_jb=np.array([10.0, 50.0]),
    ...     v_s30=300.0, mechanism=np.array(["SS", "RS"]))
    >>> pygmm.AkkarSandikkayaBommer2014(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.AkkarSandikkayaBommer2014(s, ims=["pga", "psa_1p000"])
    >>> m.spec_accels.shape
    (2, 1)

    """

    NAME = "Akkar, Sandikkaya, & Bommer (2014)"
    ABBREV = "ASB14"

    # Reference velocity (m/sec)
    V_REF = 750.0

    # Load the coefficients for the model
    COEFF = collections.OrderedDict(
        (k, model.load_data_file(f"akkar-sandikkaya-bommer-2014-{k}.csv", 2))
        for k in ["dist_jb", "dist_hyp", "dist_epi"]
    )
    PERIODS = np.array(COEFF["dist_jb"].period)

    INDICES_PSA = np.arange(2, 64)
    INDEX_PGA = 0
    INDEX_PGV = 1
    PARAMS = [
        model.NumericParameter("dist_jb", False, 0, 200),
        model.NumericParameter("dist_epi", False, 0, 200),
        model.NumericParameter("dist_hyp", False, 0, 200),
        model.NumericParameter("mag", True, 4, 8),
        model.NumericParameter("v_s30", True, 150, 1200),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def __init__(self, scenario: model.Scenario, ims=None):
        """Initialize the model.

        Args:
            scenario (:class:`pygmm.model.Scenario`): earthquake scenario.
            ims (str or sequence of str, optional): intensity measures to
                compute. If *None* (default), all intensity measures are
                computed.
        """
        super().__init__(scenario, ims)

        s = self._scenario
        for k in self.COEFF:
            if s[k] is not None:
                # Scenario values get a trailing axis to broadcast against
                # the coefficients, which have one value per period
                dist = model.as_column(s[k])
                coeff = self.COEFF[k]
                break
        else:
            raise NotImplementedError("Must provide at least one distance metric.")

        # The nonlinear site term depends on PGA at the reference condition,
        # which only needs the PGA coefficients
        c_ref = coeff[[self.INDEX_PGA]]
        pga_ref = np.exp(self._calc_ln_resp_ref(c_ref, dist)[..., 0])

        c = self._coeff_rows(coeff)
        ln_resp_ref = self._calc_ln_resp_ref(c, dist)

        # Compute the nonlinear site term
        v_s30 = model.as_column(s.v_s30)
        pga_ref = model.as_column(pga_ref)
        vs_ratio = v_s30 / self.V_REF
        vs_ratio_n = vs_ratio**c.n
        # Nonlinear site term for v_s30 at or below the reference velocity
        site_nl = c.b_1 * np.log(vs_ratio) + c.b_2 * np.log(
            (pga_ref + c.c * vs_ratio_n) / ((pga_ref + c.c) * vs_ratio_n)
        )
        # Linear site term for v_s30 above the reference velocity
        site_lin = c.b_1 * np.log(np.minimum(v_s30, c.v_con) / self.V_REF)
        site = np.where(v_s30 <= self.V_REF, site_nl, site_lin)

        self._ln_resp = ln_resp_ref + site
        # The standard deviation does not depend on the scenario
        self._ln_std = np.broadcast_to(np.array(c.sd_total), self._ln_resp.shape)

    def _calc_ln_resp_ref(self, c, dist: np.ndarray) -> np.ndarray:
        """Calculate the natural logarithm of the reference response.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute.
        dist : array_like
            distance (km) with a trailing axis for the periods.

        Returns
        -------
        ln_resp_ref : class:`np.array`:
            natural log of the response at the reference condition

        """
        s = self._scenario
        mag = model.as_column(s.mag)
        mechanism = model.as_column(s.mechanism)

        ln_resp_ref = (
            c.a_1
            + c.a_3 * (8.5 - mag) ** 2
            + (c.a_4 + c.a_5 * (mag - c.c_1)) * np.log(np.sqrt(dist**2 + c.a_6**2))
        )
        ln_resp_ref = ln_resp_ref + np.where(
            mag <= c.c_1, c.a_2 * (mag - c.c_1), c.a_7 * (mag - c.c_1)
        )

        # Strike-slip has no mechanism term
        ln_resp_ref = ln_resp_ref + np.select(
            [model.equals(mechanism, "NS"), model.equals(mechanism, "RS")],
            [c.a_8, c.a_9],
            default=0.0,
        )
        return ln_resp_ref
