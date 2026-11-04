"""Newmark rigid sliding-block analysis of an acceleration time series.

Uses :mod:`numba` to compile the integration when it is installed
(``pip install pygmm[slope]``); otherwise a pure-Python loop is used, which
gives the same result more slowly.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
import scipy.constants
from scipy.integrate import cumulative_trapezoid

__author__ = "Albert Kottke"

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only without numba

    def njit(func):  # type: ignore[no-redef]
        """Return the function unchanged, as numba is not installed."""
        return func


@njit
def _calc_block_velocity(time_step, accels, yield_coef):
    """Compute the velocity of a sliding block.

    The calculation is adapted from Slammer's source code [1]_. However, the
    displacement is computed outside of the numerical integration.

    .. [1] https://github.com/mjibson/slammer/blob/master/programs/slammer/analysis/Decoupled.java#L135

    Parameters
    ----------
    time_step : float
        Time step [s].
    accels : array_like
        Acceleration time series [cm/s²].
    yield_coef : float
        Yield acceleration of the slope [cm/s²].

    Returns
    -------
    :class:`numpy.ndarray`
        Velocity time series of the sliding mass [cm/s].
    """
    vels = np.zeros_like(accels)
    sliding = False

    for i in range(1, len(accels)):
        accel_avg = 0.5 * (accels[i] + accels[i - 1])

        if sliding:
            vels[i] = vels[i - 1] + (accel_avg - yield_coef) * time_step
        else:
            vels[i] = 0

        # Check if sliding
        if sliding:
            if vels[i] <= 0:
                sliding = False
                vels[i] = 0
        else:
            if accels[i] > yield_coef:
                sliding = True

    return vels


def calc_rigid_disp(
    time_step: float, accels: npt.ArrayLike, yield_coef: float, invert: bool = False
) -> dict[str, npt.NDArray[np.floating]]:
    """Compute the displacement and velocity of a rigid sliding mass.

    Parameters
    ----------
    time_step : float
        Time step [s].
    accels : array_like
        Acceleration time series [g].
    yield_coef : float
        Yield coefficient of the slope [g].
    invert : bool, default False
        If the acceleration time series should be inverted.

    Returns
    -------
    dict
        ``"disps"``: displacement time series of the sliding mass [cm], and
        ``"vels"``: velocity time series of the sliding mass [cm/s].
    """
    # Convert into cm/s/s
    scale = scipy.constants.g / scipy.constants.centi

    yield_accel = yield_coef * scale
    # Only invert the acceleration time series
    if invert:
        scale *= -1
    # Don't update the array in place
    accels = np.asarray(accels, dtype=float) * scale

    vels = _calc_block_velocity(time_step, accels, yield_accel)
    disps = np.r_[0, cumulative_trapezoid(vels, dx=time_step)]

    return {"disps": disps, "vels": vels}
