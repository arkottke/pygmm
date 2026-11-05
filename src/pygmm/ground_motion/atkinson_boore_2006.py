"""Atkinson and Boore (2006, :cite:p:`atkinson06`) model."""

from typing import Any

import numpy as np

from .. import model
from .._types import ArrayLike
from ..registry import register

__author__ = "Albert Kottke"


@register(tectonic="stable_continental")
class AtkinsonBoore2006(model.GroundMotionModel):
    """Atkinson and Boore (2006, :cite:p:`atkinson06`) model.

    Developed for the Eastern North America with a reference velocity of 760
    or 2000 m/s. Hard-rock sites (``v_s30`` of 2000 m/s or greater) use the
    hard-rock coefficients without site amplification. Other sites use the B/C
    coefficients (760 m/s) with the site amplification of Eq. 7 and 8.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``, and
    ``v_s30``) can be a scalar or an array, and the arrays are broadcast
    against each other. For a scalar scenario, the response and standard
    deviation have one value per period, as in other models. For arrays of N
    scenarios, the response and standard deviation have shape (N, periods),
    so, for example, ``pga`` has shape (N,) and ``spec_accels`` has shape
    (N, 24).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", "pgd", spectral periods
        such as "psa_1p000" (1.0 s), and/or "psa_all" (all 24 periods).
        Computing only the needed intensity measures is much faster for large
        vectorized scenarios. If *None* (default), all intensity measures are
        computed.

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_rup=np.array([10.0, 150.0]),
    ...     v_s30=np.array([300.0, 760.0]))
    >>> pygmm.AtkinsonBoore2006(s).pga.shape
    (2,)
    >>> pygmm.AtkinsonBoore2006(s, ims=["pga"]).pga.shape
    (2,)
    >>> m = pygmm.AtkinsonBoore2006(s, ims=["pga", "psa_0p100", "psa_1p000"])
    >>> m.psa_ims
    ['psa_0p100', 'psa_1p000']
    >>> m.spec_accels.shape
    (2, 2)

    """

    NAME = "Atkinson and Boore (2006)"
    ABBREV = "AB06"

    # Load the coefficients for the model
    COEFF = dict(
        bc=model.load_data_file("atkinson_boore_2006-bc.csv", 2),
        rock=model.load_data_file("atkinson_boore_2006-rock.csv", 2),
    )

    PERIODS = COEFF["bc"]["period"]

    COEFF_SITE = model.load_data_file("atkinson_boore_2006-site.csv", 2)
    COEFF_SF = model.load_data_file("atkinson_boore_2006-sf.csv", 2)

    INDEX_PGD = 0
    INDEX_PGV = 1
    INDEX_PGA = 2
    INDICES_PSA = np.arange(3, 27)

    PARAMS = [
        model.NumericParameter("mag", True),
        model.NumericParameter("dist_rup", True),
        model.NumericParameter("v_s30", True),
    ]

    def __init__(self, scenario: model.Scenario, ims=None):
        """Initialize the model."""
        super().__init__(scenario, ims)
        self._ln_resp = self._calc_ln_resp()
        self._ln_std = self._calc_ln_std()

    def _calc_ln_resp(self) -> np.ndarray:
        """Calculate the natural logarithm of the response.

        Returns
        -------
        ln_resp : class:`np.array`:
            natural log of the response
        """
        s = self._scenario
        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per period
        v_s30 = model.as_column(s.v_s30)
        periods = self._coeff_rows(self.PERIODS)

        # Compute the log10 response in units of cm/sec/sec at the B/C boundary,
        # including the stress drop correction
        log10_resp = self._calc_log10_ref(self._coeff_rows(self.COEFF["bc"]), periods)

        # The hard-rock coefficients without site amplification are used for
        # hard-rock sites (v_s30 >= 2000 m/s), and for v_s30 of zero as in
        # previous versions. The site amplification is computed with the
        # reference velocity there to avoid the logarithm of zero, and then
        # replaced.
        is_rock = (v_s30 == 0) | (v_s30 >= 2000)
        if np.any(is_rock):
            v_s30 = np.where(is_rock, 760.0, v_s30)

        # Compute the site amplification, which depends on PGA at the B/C
        # boundary (cm/sec/sec). PGA is taken from the response if it was
        # computed, and otherwise computed with only the PGA coefficients.
        index_pga = [self.INDEX_PGA]
        if self._indices is None:
            log10_pga = log10_resp[..., index_pga]
        elif self.INDEX_PGA in self._indices:
            log10_pga = log10_resp[..., np.searchsorted(self._indices, index_pga)]
        else:
            log10_pga = self._calc_log10_ref(
                self.COEFF["bc"][index_pga], self.PERIODS[index_pga]
            )
        pga_bc = 10**log10_pga
        log10_site = self._calc_log10_site(pga_bc, v_s30, periods)
        log10_resp = log10_resp + log10_site

        if np.any(is_rock):
            log10_rock = self._calc_log10_ref(
                self._coeff_rows(self.COEFF["rock"]), periods
            )
            log10_resp = np.where(is_rock, log10_rock, log10_resp)

        # Convert PGA and PSA from cm/sec/sec to gravity. PGV (cm/sec) and PGD
        # (cm), which have negative periods, are not converted.
        log10_resp = log10_resp - np.where(periods >= 0, np.log10(980.665), 0.0)

        ln_resp = np.log(10**log10_resp)
        return ln_resp

    def _calc_log10_ref(self, c: Any, periods: ArrayLike) -> np.ndarray:
        """Calculate the log10 response [cm/sec/sec] at the reference condition.

        Parameters
        ----------
        c : :class:`numpy.recarray`
            coefficients for the periods to compute.
        periods : array_like
            periods of the coefficients, used for the stress factor.

        Returns
        -------
        log10_resp : class:`np.array`:
            log base 10 of the response, including the stress drop correction
        """
        s = self._scenario
        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)

        # Compute the response at the reference condition
        r0 = 10.0
        r1 = 70.0
        r2 = 140.0

        f0 = np.maximum(np.log10(r0 / dist_rup), 0)
        f1 = np.minimum(np.log10(dist_rup), np.log10(r1))
        f2 = np.maximum(np.log10(dist_rup / r2), 0)

        # Compute the log10 PSA in units of cm/sec/sec
        log10_resp = (
            c.c_1
            + c.c_2 * mag
            + c.c_3 * mag**2
            + (c.c_4 + c.c_5 * mag) * f1
            + (c.c_6 + c.c_7 * mag) * f2
            + (c.c_8 + c.c_9 * mag) * f0
            + c.c_10 * dist_rup
        )

        # Apply stress drop correction
        log10_resp = log10_resp + self._calc_stress_factor(periods)

        return log10_resp

    def _calc_ln_std(self) -> np.ndarray:
        """Calculate the logarithmic standard deviation.

        Returns
        -------
        ln_std : class:`np.array`:
            natural log standard deviation

        """
        # Constant standard deviation of 0.30 in log10 units, converted to natural
        # log units, with the same shape as the response
        ln_std = np.full(np.shape(self._resp_vector), 0.30 * np.log(10))
        return ln_std

    def _calc_stress_factor(self, periods: ArrayLike) -> np.ndarray:
        """Calculate the stress correction factor.

        The factor was proposed by Atkinson and Boore (2011)
        :cite:p:`atkinson11`.

        Parameters
        ----------
        periods : array_like
            periods at which the stress factor is computed.

        Returns
        -------
        log10_stress_factor : class:`np.array`:
            log base 10 of the stress factor
        """
        s = self._scenario
        c = self.COEFF_SF
        if self._indices is not None:
            # Only the rows needed to interpolate at the periods are computed
            c = c[interp_rows(periods, c.period)]
        mag = model.as_column(s.mag)

        stress_drop = 10.0 ** (3.45 - 0.2 * mag)
        v1 = c.delta + 0.05
        v2 = 0.05 + c.delta * np.maximum(mag - c.m_1, 0) / (c.m_h - c.m_1)

        log10_stress_factor = np.minimum(2.0, stress_drop / 140.0) * np.minimum(v1, v2)

        return interp_last_axis(periods, c.period, log10_stress_factor)

    def _calc_log10_site(
        self, pga_bc: np.ndarray, v_s30: np.ndarray, periods: ArrayLike
    ) -> np.ndarray:
        """Calculate the log10 of the site amplification.

        Parameters
        ----------
        pga_bc : array_like
            peak ground acceleration (PGA, cm/sec/sec) at the B/C boundary,
            with a trailing axis of length one.
        v_s30 : array_like
            time-averaged shear-wave velocity [m/sec], with a trailing axis of
            length one.
        periods : array_like
            periods at which the site amplification is computed.

        Returns
        -------
        log_10_site : :class:`np.ndarray`
            log base 10 of the  site amplification.
        """
        c = self.COEFF_SITE
        if self._indices is not None:
            # Only the rows needed to interpolate at the periods are computed
            c = c[interp_rows(periods, c.period)]
        VS_1 = 180.0
        VS_2 = 300.0
        VS_REF = 760.0

        b_nl = np.select(
            [v_s30 <= VS_1, v_s30 <= VS_2, v_s30 <= VS_REF],
            [
                c.b_1,
                (c.b_1 - c.b_2) * np.log(v_s30 / VS_2) / np.log(VS_1 / VS_2) + c.b_2,
                c.b_2 * np.log(v_s30 / VS_REF) / np.log(VS_2 / VS_REF),
            ],
            # Vs30 > VS_REF
            default=0,
        )

        pga_bc = np.maximum(pga_bc, 60.0)

        log10_site = np.log10(
            np.exp(c.b_lin * np.log(v_s30 / VS_REF) + b_nl * np.log(pga_bc / 100.0))
        )
        return interp_last_axis(periods, c.period, log10_site)


def _interval(x: np.ndarray, xp: np.ndarray) -> np.ndarray:
    """Index j of the interval xp[j] <= x < xp[j + 1], limited to valid intervals."""
    return np.clip(np.searchsorted(xp, x, side="right") - 1, 0, len(xp) - 2)


def interp_rows(x: ArrayLike, xp: ArrayLike) -> np.ndarray:
    """Rows of `xp` needed to interpolate at `x` with :func:`interp_last_axis`.

    Interpolating with only these rows gives the same values as with all rows,
    because the interval of each point in `x` is kept.

    Parameters
    ----------
    x : array_like
        x-coordinates at which to interpolate.
    xp : array_like
        increasing x-coordinates of the data points.

    Returns
    -------
    rows : :class:`np.ndarray`
        sorted indices of the rows
    """
    j = _interval(np.atleast_1d(np.asarray(x, dtype=float)), np.asarray(xp))
    return np.unique(np.concatenate([j, j + 1]))


def interp_last_axis(x: ArrayLike, xp: ArrayLike, fp: ArrayLike) -> np.ndarray:
    """Linearly interpolate along the last axis, like :func:`numpy.interp`.

    :func:`numpy.interp` only interpolates 1-D values. This gives the same values
    for each scenario of `fp` with shape (..., len(xp)).

    Parameters
    ----------
    x : array_like
        x-coordinates at which to interpolate.
    xp : array_like
        increasing x-coordinates of the data points.
    fp : array_like
        values at `xp` along the last axis.

    Returns
    -------
    values : :class:`np.ndarray`
        interpolated values with shape (..., len(x)).
    """
    fp = np.asarray(fp)
    if fp.ndim <= 1:
        return np.interp(x, xp, fp)

    x = np.asarray(x, dtype=float)
    xp = np.asarray(xp, dtype=float)
    j = _interval(x, xp)
    # Points outside of xp take the end values, and points on xp take the values
    # without interpolation, as in numpy.interp
    exact = np.select([x < xp[0], x >= xp[-1]], [0, len(xp) - 1], default=j)
    is_exact = (x < xp[0]) | (x >= xp[-1]) | (xp[j] == x)

    slope = (fp[..., j + 1] - fp[..., j]) / (xp[j + 1] - xp[j])
    return np.where(is_exact, fp[..., exact], slope * (x - xp[j]) + fp[..., j])
