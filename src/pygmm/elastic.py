"""Relationships between the isotropic elastic constants.

See https://en.wikipedia.org/wiki/Elastic_modulus. The moduli may be in any
consistent units.
"""

from __future__ import annotations

import numpy as np

__author__ = "Albert Kottke"


def calc_mod_bulk(mod_shear: float | np.ndarray, poissons_ratio: float | np.ndarray):
    """Compute the bulk modulus from the shear modulus and Poisson's ratio.

    Parameters
    ----------
    mod_shear : float or array_like
        Shear modulus.
    poissons_ratio : float or array_like
        Poisson's ratio.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Bulk modulus, in the units of `mod_shear`.
    """
    return (2 * mod_shear * (1 + poissons_ratio)) / (3 * (1 - 2 * poissons_ratio))


def calc_mod_shear(mod_bulk: float | np.ndarray, poissons_ratio: float | np.ndarray):
    """Compute the shear modulus from the bulk modulus and Poisson's ratio.

    Parameters
    ----------
    mod_bulk : float or array_like
        Bulk modulus.
    poissons_ratio : float or array_like
        Poisson's ratio.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Shear modulus, in the units of `mod_bulk`.
    """
    return (3 * mod_bulk * (1 - 2 * poissons_ratio)) / (2 * (1 + poissons_ratio))


def calc_poissons_ratio(mod_bulk: float | np.ndarray, mod_shear: float | np.ndarray):
    """Compute Poisson's ratio from the bulk and shear moduli.

    Parameters
    ----------
    mod_bulk : float or array_like
        Bulk modulus.
    mod_shear : float or array_like
        Shear modulus, in the units of `mod_bulk`.

    Returns
    -------
    float or :class:`numpy.ndarray`
        Poisson's ratio.
    """
    return (3 * mod_bulk - 2 * mod_shear) / (2 * (3 * mod_bulk + mod_shear))
