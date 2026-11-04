"""Models of the site attenuation parameter, :math:`\\kappa_0`.

Each relationship follows the function-style convention of
:func:`pygmm._utils.dist_lognorm`: it returns the median :math:`\\kappa_0`
[s], or ``(ln_mean, ln_std)`` when called with ``stats=True``.
"""

from .van_houtte_et_al_2011 import calc_site_atten_vhea11
from .xu_rathje_2021 import calc_site_atten_xr21

__all__ = ["calc_site_atten_vhea11", "calc_site_atten_xr21"]
