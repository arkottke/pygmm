"""Abrahamson and Silva (1996, :cite:p:`abrahamson96`) duration model."""

import numpy as np

from .. import model
from ..contracts import Duration
from ..registry import register

__author__ = "Albert Kottke"


@register(provides=("duration",))
class AbrahamsonSilva1996(model.Model):
    """Abrahamson and Silva (1996, :cite:p:`abrahamson96`) duration model.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``, and
    ``site_cond``) can be a scalar or an array, and the arrays are broadcast
    against each other. For a scalar scenario, :attr:`duration` is a float and
    :attr:`std_err` is 0.55, as before. For arrays of N scenarios,
    :attr:`duration` and :attr:`std_err` have shape (N,), and :meth:`interp`
    returns an array with shape (N, len(nias)), or (N, len(stds), len(nias))
    if `stds` is provided. Multidimensional scenario arrays keep their shape.

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
    ...     site_cond=np.array(["soil", "rock", "soil"]))
    >>> m = pygmm.AbrahamsonSilva1996(s)
    >>> m.duration.shape
    (3,)
    >>> m.interp([0.5, 0.75, 0.95]).shape
    (3, 3)
    >>> m.interp([0.5, 0.75, 0.95], stds=[-1, 0, 1]).shape
    (3, 3, 3)

    """

    NAME = "Abrahamson Silva (1996)"
    ABBREV = "AS96"

    PARAMS = [
        model.NumericParameter("mag", True, 4, 7.5),
        model.NumericParameter("dist_rup", True, 0, 250),
        # FIXME add site_cond to Scenario
        model.CategoricalParameter("site_cond", True, ["soil", "rock"]),
    ]

    # Normalized Arias intensities and the corresponding standard errors used by
    # interp()
    INTERP_NIAS = [
        0.10,
        0.15,
        0.20,
        0.25,
        0.30,
        0.35,
        0.40,
        0.45,
        0.50,
        0.55,
        0.60,
        0.65,
        0.70,
        0.75,
        0.80,
        0.85,
        0.90,
        0.95,
    ]
    INTERP_STD_ERRS = [
        0.843,
        0.759,
        0.713,
        0.691,
        0.674,
        0.660,
        0.646,
        0.636,
        0.628,
        0.616,
        0.605,
        0.594,
        0.582,
        0.565,
        0.545,
        0.528,
        0.510,
        0.493,
    ]

    def __init__(self, scenario):
        super().__init__(scenario)

        s = self._scenario

        stress_drop = np.exp(5.204 + 0.851 * (s.mag - 6))
        moment = 10 ** (1.5 * s.mag + 16.05)

        is_soil = np.where(model.equals(s.site_cond, "soil"), 1, 0)

        # Natural log of the duration without the increment for the normalized
        # Arias intensity, which is added for the 5 to 75% duration and by interp()
        self._ln_dur_base = np.log(
            (stress_drop / moment) ** (-1 / 3) / (4.9e6 * 3.2)
            + 0.805 * is_soil
            + 0.063 * np.maximum(s.dist_rup - 10, 0)
        )
        self._ln_dur = self._ln_dur_base + self.calc_ln_dur_incr(0.75)
        self._std_err = 0.55

    @property
    def duration(self):
        """Duration from 5 to 75% of the normalized Arias intensity [sec]."""
        return np.exp(self._ln_dur)

    @property
    def std_err(self):
        """Logarithmic standard error of the duration.

        A float for a scalar scenario, and an array with the shape of
        :attr:`duration` for a vectorized scenario.
        """
        if np.ndim(self._ln_dur) > 0:
            return np.full(np.shape(self._ln_dur), self._std_err)
        return self._std_err

    @staticmethod
    def calc_ln_dur_incr(nias):
        """Calculate the increment in the log duration relative to the start.

        Parameters
        ----------
        nias : float or array_like
            normalized Arias intensities (between 0.10 and 0.95). Values
            outside of this range return NaN. The input is not modified.

        Returns
        -------
        ln_dur_incr : float or :class:`np.ndarray`
            increment in the natural log of the duration, with the shape of
            `nias`
        """
        # Mask out inappropriate values on a copy, so the caller's values are not
        # modified
        nias = np.array(nias, dtype=float)
        mask = (nias < 0.10) | (0.95 < nias)
        nias[mask] = np.nan

        # Compute the increment due to the difference in the normalized Arias
        # intensity
        ln_i_ratio = np.log((nias - 0.05) / (1 - nias))
        ln_dur_incr = -0.532 + 0.552 * ln_i_ratio - 0.0262 * ln_i_ratio**2

        return ln_dur_incr

    def interp(self, nias, stds=None):
        """Duration to normalized Arias intensities.

        Parameters
        ----------
        nias : float or array_like
            normalized Arias intensities (between 0.10 and 0.95). Values
            outside of this range return NaN. The input is not modified.
        stds : array_like, optional
            number of standard deviations from the median.

        Returns
        -------
        durations : :class:`np.ndarray`
            durations [sec]. For a scalar scenario, the shape is that of `nias`,
            or (len(stds), len(nias)) if `stds` is provided. For a vectorized
            scenario with N scenarios, the shape is (N, len(nias)), or
            (N, len(stds), len(nias)) if `stds` is provided.
        """
        # Mask out inappropriate values on a copy, so the caller's values are not
        # modified
        nias = np.array(nias, dtype=float)
        mask = (nias < 0.10) | (0.95 < nias)
        nias[mask] = np.nan

        ln_dur_incr = self.calc_ln_dur_incr(nias)

        if np.ndim(self._ln_dur_base) == 0:
            ln_dur = self._ln_dur_base + ln_dur_incr
            if stds is not None:
                std_errs = np.interp(nias, self.INTERP_NIAS, self.INTERP_STD_ERRS)
                ln_dur = ln_dur + np.array(stds)[:, np.newaxis] * std_errs
            return np.exp(ln_dur)

        # Vectorized scenario: the axes of the scenarios are followed by the axes
        # of the standard deviations (if provided) and the normalized Arias
        # intensities
        ln_dur = self._ln_dur_base
        if stds is None:
            ln_dur = np.reshape(ln_dur, ln_dur.shape + (1,) * nias.ndim)
            ln_dur = ln_dur + ln_dur_incr
        else:
            stds = np.asarray(stds, dtype=float)
            std_errs = np.interp(nias, self.INTERP_NIAS, self.INTERP_STD_ERRS)
            ln_dur = np.reshape(ln_dur, ln_dur.shape + (1,) * (stds.ndim + nias.ndim))
            ln_dur = ln_dur + ln_dur_incr
            ln_dur = ln_dur + np.reshape(stds, stds.shape + (1,) * nias.ndim) * std_errs
        return np.exp(ln_dur)

    def duration_model(self) -> Duration:
        """Return the :class:`~pygmm.contracts.Duration` contract.

        The model predicts 5-75% significant duration -- ``__init__`` builds
        ``_ln_dur`` with ``calc_ln_dur_incr(0.75)`` -- with a lognormal
        standard error.

        The values are floats for a scalar scenario and arrays for a vectorized
        scenario.
        """
        dur = model.as_scalar_or_array(self.duration)
        return Duration(
            duration=dur, d_5_75=dur, ln_std=model.as_scalar_or_array(self.std_err)
        )
