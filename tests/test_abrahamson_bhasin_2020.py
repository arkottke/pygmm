"""Tests for Abrahamson and Bhasin (2020) conditional PGV model."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm import Scenario
from pygmm.ground_motion.abrahamson_bhasin_2020 import AbrahamsonBhasin2020


class TestConstruction:
    def test_instantiates(self):
        """Regression: ``__init__`` used to ``return`` a tuple, raising TypeError."""
        s = Scenario(mag=7.0, dist_rup=10.0, v_s30=760.0)
        assert isinstance(AbrahamsonBhasin2020(s, pga=0.5), AbrahamsonBhasin2020)

    def test_no_conditioning_im_raises(self):
        s = Scenario(mag=7.0, dist_rup=10.0, v_s30=760.0)
        with pytest.raises(ValueError, match="pga.*psa_1s"):
            AbrahamsonBhasin2020(s)

    def test_psa_mode_not_implemented(self):
        s = Scenario(mag=7.0, dist_rup=10.0, v_s30=760.0)
        with pytest.raises(NotImplementedError):
            AbrahamsonBhasin2020(s, psa=0.5)


class TestConditionedOnPga:
    """Mw=7.0, Rrup=10 km, Vs30=760 m/s, PGA=0.5 g.

    Reference values evaluated by hand from Eq. 3.4 with the Table 3.2 "pga"
    coefficients:

        f_1     = 0.738 + (0.484 - 0.738) * (7 - 5) / 2.5     = 0.5348
        ln_pgv  = 4.77 + 0.5348*ln(0.5) + 0.275*(7-6)
                  - 0.036*(8.5-7)**2 - 0.332*ln(10 + 5*exp(0.4))
                  - 0.44*ln(760/425)                          = 3.388089712
        phi, tau at Mw=7.0 = M_2 are the upper anchors 0.42, 0.26
        ln_std  = sqrt(0.42**2 + 0.26**2)                     = 0.493963561
    """

    @pytest.fixture
    def m(self):
        s = Scenario(mag=7.0, dist_rup=10.0, v_s30=760.0)
        return AbrahamsonBhasin2020(s, pga=0.5)

    def test_ln_pgv(self, m):
        assert_allclose(m.ln_pgv, 3.388089712, rtol=1e-9)

    def test_pgv(self, m):
        assert_allclose(m.pgv, 29.609335851, rtol=1e-9)

    def test_sigma_components(self, m):
        assert m.phi == pytest.approx(0.42)
        assert m.tau == pytest.approx(0.26)
        assert_allclose(m.ln_std, np.sqrt(0.42**2 + 0.26**2), rtol=1e-10)

    def test_plus_minus_sigma(self, m):
        assert_allclose(m.pgv_plus_sigma, np.exp(m.ln_pgv + m.ln_std))
        assert_allclose(m.pgv_minus_sigma, np.exp(m.ln_pgv - m.ln_std))


class TestConditionedOnPsa1s:
    """Mw=6.0, Rrup=20 km, Vs30=500 m/s, Sa(1s)=0.15 g.

    Hand evaluation with the Table 3.2 "psa(T=1s)" coefficients:

        f_1     = 0.82 + (0.55 - 0.82) * (6 - 5) / 2.5        = 0.712
        ln_pgv  = 4.80 + 0.712*ln(0.15) + 0.27*(6-6)
                  + 0.054*(8.5-6)**2 - 0.382*ln(20 + 5)
                  - 0.21*ln(500/425)                          = 2.523011030
        phi     = 0.28 + (0.38 - 0.28) * (6 - 5) / 2          = 0.33
        tau     = 0.12 + (0.17 - 0.12) * (6 - 5) / 2          = 0.145
    """

    @pytest.fixture
    def m(self):
        s = Scenario(mag=6.0, dist_rup=20.0, v_s30=500.0)
        return AbrahamsonBhasin2020(s, psa_1s=0.15)

    def test_ln_pgv(self, m):
        assert_allclose(m.ln_pgv, 2.523011030, rtol=1e-9)

    def test_sigma_components(self, m):
        assert m.phi == pytest.approx(0.33)
        assert m.tau == pytest.approx(0.145)


class TestScaling:
    COMMON = dict(dist_rup=20.0, v_s30=760.0)

    def test_pgv_increases_with_conditioning_im(self):
        s = Scenario(mag=6.5, **self.COMMON)
        assert (
            AbrahamsonBhasin2020(s, pga=0.4).pgv > AbrahamsonBhasin2020(s, pga=0.1).pgv
        )

    def test_pgv_decreases_with_distance(self):
        near = Scenario(mag=6.5, dist_rup=5.0, v_s30=760.0)
        far = Scenario(mag=6.5, dist_rup=100.0, v_s30=760.0)
        assert (
            AbrahamsonBhasin2020(near, pga=0.2).pgv
            > AbrahamsonBhasin2020(far, pga=0.2).pgv
        )

    def test_pgv_increases_on_softer_site(self):
        """a_7 < 0, so lower Vs30 gives larger PGV at fixed conditioning IM."""
        soft = Scenario(mag=6.5, dist_rup=20.0, v_s30=270.0)
        stiff = Scenario(mag=6.5, dist_rup=20.0, v_s30=1000.0)
        assert (
            AbrahamsonBhasin2020(soft, pga=0.2).pgv
            > AbrahamsonBhasin2020(stiff, pga=0.2).pgv
        )

    def test_sigma_saturates_outside_mag_range(self):
        """phi/tau are clamped to the anchors below M_1 and above M_2."""
        low = AbrahamsonBhasin2020(Scenario(mag=4.0, **self.COMMON), pga=0.05)
        high = AbrahamsonBhasin2020(Scenario(mag=8.0, **self.COMMON), pga=0.6)
        assert (low.phi, low.tau) == pytest.approx((0.32, 0.12))
        assert (high.phi, high.tau) == pytest.approx((0.42, 0.26))


def test_ln_period_pgv():
    """Equation 3.6."""
    assert AbrahamsonBhasin2020.ln_period_pgv(7.0) == pytest.approx(-4.09 + 0.66 * 7.0)
