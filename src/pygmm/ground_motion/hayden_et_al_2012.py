"""Hayden et al. (2012) proportion of pulse-like records."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__author__ = "Albert Kottke"


def calc_pulse_proportion_hea12(
    dist_rup: npt.ArrayLike, epsilon: npt.ArrayLike
) -> npt.NDArray[np.floating] | float:
    """Calculate the pulse proportion from Hayden et al. (2012).

    Parameters
    ----------
    dist_rup : float or array_like
        Closest distance to the fault rupture [km].
    epsilon : float or array_like
        Number of standard deviations.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Proportion of records with pulse-like behavior.
    """
    dist_rup = np.asarray(dist_rup)
    epsilon = np.asarray(epsilon)

    # FIXME: The intercepts of the numerator (0.891) and denominator (0.981)
    # differ, which is unexpected for a logistic model -- check against the
    # publication. Carried over unchanged from pyGEEMs.
    return np.exp(0.891 - 0.188 * dist_rup + 1.230 * epsilon) / (
        1 + np.exp(0.981 - 0.188 * dist_rup + 1.230 * epsilon)
    )
