"""Cho and Rathje (2021) neural-network sliding-mass displacement model."""

from __future__ import annotations

import numpy as np

from .._utils import dist_lognorm

__author__ = "Albert Kottke"

# Means from Table 2, but the mean of the ln(PGV) is replaced based on
# Equation (2)
_LN_MEANS = np.array([5.855, -3.873, -1.498, -2.262])
# Bias from Table 3
_B_OUT = -4.687
# Weights and bias of each hidden neuron (columns 0-4), and the output weight
# (column 5), from Table 3
_COEFFS = np.array(
    [
        [5.262, 0.66, -2.709, -0.354, 0.505, 3.024],
        [1.678, 1.054, 1.686, -4.184, -2.9, 3.303],
        [4.253, 0.527, -0.876, -0.265, -2.818, 3.044],
        [0.82, 7.245, 1.074, 6.094, -3.828, 2.14],
        [4.05, -0.154, -1.238, 0.907, -2.904, 3.449],
    ]
)


@dist_lognorm
def calc_disp_cr21(
    pgv: float, yield_coef: float, period_slide: float, height_ratio: float, **kwds
):
    """Artificial neural network Model 1 from Cho and Rathje (2021).

    Parameters
    ----------
    pgv : float
        Peak ground velocity [cm/s].
    yield_coef : float
        Yield coefficient of the slope [g].
    period_slide : float
        Period of the sliding mass [s].
    height_ratio : float
        Height ratio of the sliding mass.
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float
        Median displacement [cm], or ``(ln_mean, ln_std)`` if `stats` is True.
    """
    # Median prediction, Equation 2
    x = np.r_[
        np.log([pgv, yield_coef, period_slide, height_ratio]) / _LN_MEANS,
        # Dummy value for including the bias of each neuron
        1,
    ]
    ln_mean = _B_OUT + _COEFFS[:, -1] @ (1 / (1 + np.exp(-_COEFFS[:, :-1] @ x)))

    # Standard deviation model. Equations 6 & 7.
    phi = 0.27 - 0.137 * np.log(period_slide)
    tau = 0.106
    ln_std = np.sqrt(phi**2 + tau**2)

    return ln_mean, ln_std
