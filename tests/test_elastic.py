"""Tests for pygmm.elastic."""

import numpy as np

from pygmm.elastic import calc_mod_bulk, calc_mod_shear, calc_poissons_ratio


def test_elastic_moduli():
    vel_shear = 2000
    poissons_ratio = 0.3
    dens = 2.3
    mod_shear = dens * vel_shear**2

    mod_bulk = calc_mod_bulk(mod_shear, poissons_ratio)

    mod_shear_calc = calc_mod_shear(mod_bulk, poissons_ratio)
    np.testing.assert_allclose(mod_shear_calc, mod_shear)

    poissons_ratio_calc = calc_poissons_ratio(mod_bulk, mod_shear)
    np.testing.assert_allclose(poissons_ratio_calc, poissons_ratio)
