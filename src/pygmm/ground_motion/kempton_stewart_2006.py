"""Kempton and Stewart (2006, :cite:p:`kempton2006`) duration model."""

import numpy as np

from .. import model
from ..contracts import Duration
from ..registry import register

__author__ = ""


@register(provides=("duration",))
class KemptonStewart2006(model.Model):
    """Kempton and Stewart (2006, :cite:p:`kempton2006`) duration model.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``, and
    ``v_s30``) can be a scalar or an array, and the arrays are broadcast against
    each other. :attr:`duration` and :attr:`std_err` are record arrays with the
    fields ``D_5t75a``, ``D_5t95a``, ``D_5t75v``, and ``D_5t95v``. For a scalar
    scenario, each field is a scalar, as before. For arrays of N scenarios, each
    field has shape (N,). Multidimensional scenario arrays keep their shape.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([5.5, 6.5, 7.0]), dist_rup=np.array([5.0, 20.0, 80.0]),
    ...     v_s30=300.0)
    >>> m = pygmm.KemptonStewart2006(s)
    >>> m.duration.D_5t95a.shape
    (3,)
    >>> m.std_err.D_5t95a.shape
    (3,)

    """

    NAME = "Kempton Stewart (2006)"
    ABBREV = "KS06"

    PARAMS = [
        model.NumericParameter("mag", True, 5, 7.6),
        model.NumericParameter("dist_rup", True, 0, 200),
        model.NumericParameter("v_s30", True, 200, 1000),
    ]

    def __init__(self, scenario):
        super().__init__(scenario)

        # scenario
        s = self._scenario
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per duration measure
        mag = model.as_column(s.mag)

        # stress drop indices for duration measures: D5-75a, D5-95a, D5-75v, and
        # D5-95v. Only the D5-95 measures scale with magnitude.
        stress_drop = np.exp(
            np.array([6.02, 2.79, 5.46, 1.53])
            + np.array([0.0, 0.82, 0.0, 1.34]) * (mag - 6)
        )

        # source duration
        f_0 = self._source_dur(stress_drop)

        # path duration
        f_1 = np.array([0.07, 0.15, 0.10, 0.15]) * model.as_column(s.dist_rup)

        # site duration
        f_2 = np.array([0.82, 3.00, 1.40, 3.99]) + np.array(
            [-0.0013, -0.0041, -0.0022, -0.0062]
        ) * model.as_column(s.v_s30)

        # total druation
        self._ln_dur = np.log(f_0 + f_1 + f_2)

        # aleatory standard deviation
        self._std_err = np.array([0.57, 0.47, 0.66, 0.52])

    @property
    def duration(self):
        """Median durations [s] as a `np.recarray`, accessed by the interval."""
        return KemptonStewart2006._as_recarray(np.exp(self._ln_dur))

    @property
    def std_err(self):
        """Standard errors of the natural log durations as a `np.recarray`."""
        std_err = self._std_err
        if np.ndim(self._ln_dur) > 1:
            std_err = np.broadcast_to(std_err, np.shape(self._ln_dur))
        return KemptonStewart2006._as_recarray(std_err)

    @staticmethod
    def _as_recarray(values):
        # The duration measures are along the last axis
        return np.rec.fromarrays(
            np.moveaxis(values, -1, 0),
            names=[
                "D_5t75a",
                "D_5t95a",
                "D_5t75v",
                "D_5t95v",
            ],
        )

    def _source_dur(self, stress_drop):
        # seismic moment
        moment = model.as_column(10 ** (1.5 * self.scenario.mag + 16.05))
        return (stress_drop / moment) ** (-1 / 3) / (4.9e6 * 3.2)

    def duration_model(self) -> Duration:
        """Return the :class:`~pygmm.contracts.Duration` contract.

        The model reports acceleration- and velocity-based durations; the
        acceleration-based 5-95% value is taken as the primary one.

        The values are floats for a scalar scenario and arrays for a vectorized
        scenario.
        """
        d, s = self.duration, self.std_err
        return Duration(
            duration=model.as_scalar_or_array(d["D_5t95a"]),
            d_5_75=model.as_scalar_or_array(d["D_5t75a"]),
            d_5_95=model.as_scalar_or_array(d["D_5t95a"]),
            ln_std=model.as_scalar_or_array(s["D_5t95a"]),
        )
