"""Derras, Bard and Cotton (2014, :cite:p:`derras14`) model."""

import json
import os

import numpy as np

from .. import model
from ..registry import register

__author__ = "Albert Kottke"


@register()
class DerrasBardCotton2014(model.GroundMotionModel):
    """Derras, Bard and Cotton (2014, :cite:p:`derras14`) model.

    The model is vectorized. Each scenario value (``mag``, ``dist_jb``,
    ``v_s30``, ``depth_hyp``, and ``mechanism``) can be a scalar or an array,
    and the arrays are broadcast against each other. For a scalar scenario,
    the response and standard deviation have one value per period, as in
    other models. For arrays of N scenarios, they have shape (N, periods), so,
    for example, ``pga`` has shape (N,).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 62 periods, from 0.01 to
        4 s, so "psa_ngawest2_21" is not available). The neural network
        computes all of the intensity measures at once, but only the requested
        intensity measures are converted and kept. If *None* (default), all
        intensity measures are computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([5.0, 6.0]), dist_jb=np.array([10.0, 50.0]),
    ...     v_s30=400.0, depth_hyp=10.0, mechanism=np.array(["SS", "RS"]))
    >>> pygmm.DerrasBardCotton2014(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.DerrasBardCotton2014(s, ims=["pga", "psa_1p000"])
    >>> m.spec_accels.shape
    (2, 1)

    """

    NAME = "Derras, Bard & Cotton (2014)"
    ABBREV = "DBC13"

    # Load the coefficients for the model
    COEFF = json.load(
        open(
            os.path.join(
                os.path.dirname(os.path.dirname(__file__)),
                "data",
                "derras_bard_cotton_2014.json",
            )
        )
    )
    GRAVITY = 9.80665
    PERIODS = np.array(COEFF["period"])

    INDICES_PSA = np.arange(2, 64)
    INDEX_PGA = 1
    INDEX_PGV = 0
    PARAMS = [
        model.NumericParameter("dist_jb", True, 5, 200),
        model.NumericParameter("mag", True, 4, 7),
        model.NumericParameter("v_s30", True, 200, 800),
        model.NumericParameter("depth_hyp", True, 0, 25),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def __init__(self, scenario: model.Scenario, ims=None):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed.
        """
        super().__init__(scenario, ims)
        c = self.COEFF
        # Values modified during the calculation
        s = dict(self._scenario)

        for k in ["v_s30", "dist_jb"]:
            s["log10_" + k] = np.log10(s[k])
        # Translate to mechanism integer
        mechanism = np.select(
            [
                model.equals(s["mechanism"], "NS"),
                model.equals(s["mechanism"], "RS"),
                model.equals(s["mechanism"], "SS"),
            ],
            [1, 3, 4],
            default=-1,
        )
        if np.any(mechanism < 0):
            raise ValueError("mechanism must be one of: SS, NS, RS")
        s["mechanism"] = mechanism

        # Create the normalized parameter matrix with the parameters along
        # the last axis, so N scenarios give shape (N, 5)
        keys = ["log10_dist_jb", "mag", "log10_v_s30", "depth_hyp", "mechanism"]
        values = np.stack(np.broadcast_arrays(*[s[k] for k in keys]), axis=-1).astype(
            float
        )
        limits = np.rec.array([c["min_max"][k] for k in keys], names="min,max")
        p_n = np.array(
            2 * (values - limits["min"]) / (limits["max"] - limits["min"]) - 1
        )

        # Compute the normalized response
        b_1 = np.array(c["b_1"]).T
        b_2 = np.array(c["b_2"]).T
        w_1 = np.array(c["w_1"])
        w_2 = np.array(c["w_2"])

        log10_resp_n = b_2 + self._matvec(w_2, np.tanh(b_1 + self._matvec(w_1, p_n)))
        # The network computes all periods, so select the requested ones
        if self._indices is not None:
            log10_resp_n = model.take_periods(log10_resp_n, self._indices)

        # Convert from normalized values
        log10_resp_limits = self._coeff_rows(
            np.rec.array(c["min_max"]["log10_resp"], names="min,max")
        )
        log10_resp = (
            0.5
            * (log10_resp_n + 1)
            * (log10_resp_limits["max"] - log10_resp_limits["min"])
            + log10_resp_limits["min"]
        )
        # Convert from m/sec and m/sec² into cm/sec and g
        scale = self._coeff_rows(
            np.log10(np.r_[0.01, self.GRAVITY * np.ones(self.PERIODS.size - 1)])
        )

        self._ln_resp = np.log(10 ** (log10_resp - scale))
        # The standard deviation does not depend on the scenario
        self._ln_std = np.broadcast_to(
            self._coeff_rows(np.log(10 ** np.array(c["log10_std"]["total"]))),
            self._ln_resp.shape,
        )

    @staticmethod
    def _matvec(w: np.ndarray, x: np.ndarray) -> np.ndarray:
        """Matrix-vector product of `w` with the vectors along the last axis of `x`.

        Each scenario uses the same matrix-vector product as a scalar scenario,
        so the results are identical to computing the scenarios one at a time.
        """
        if x.ndim == 1:
            return w @ x
        return (w @ x[..., np.newaxis])[..., 0]
