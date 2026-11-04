"""Xu and Rathje (2021) :math:`\\kappa_0`--:math:`V_{S30}` relationship."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .._utils import dist_lognorm

__author__ = "Albert Kottke"


@dist_lognorm
def calc_site_atten_xr21(v_s30: float | npt.ArrayLike, **kwargs):
    r"""Site attenuation from Xu and Rathje (2021).

    :math:`V_{S30}` is clipped to the range 155 to 2000 m/s.

    Parameters
    ----------
    v_s30 : float or array_like
        Time-averaged shear-wave velocity in the top 30 m [m/s].
    stats : bool, default False
        If True, return ``(ln_mean, ln_std)`` instead of the median.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Median site attenuation, :math:`\kappa_0` [s].
    """
    ln_vs30 = np.log(np.clip(v_s30, 155, 2000))
    ln_mean = -0.18 * ln_vs30**2 + 1.816 * ln_vs30 - 7.38
    ln_std = 0.3
    return ln_mean, ln_std
