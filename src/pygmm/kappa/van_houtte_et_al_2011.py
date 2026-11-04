"""Van Houtte et al. (2011) :math:`\\kappa_0`--:math:`V_{S30}` relationship."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .._utils import dist_lognorm

__author__ = "Albert Kottke"


@dist_lognorm
def calc_site_atten_vhea11(v_s30: float | npt.ArrayLike, **kwargs):
    r"""Site attenuation from Van Houtte et al. (2011).

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
    ln_mean = 3.49 - 1.062 * np.log(v_s30)
    ln_std = 0.55
    return ln_mean, ln_std
