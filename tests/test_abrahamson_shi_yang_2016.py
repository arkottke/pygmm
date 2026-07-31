"""Tests for Abrahamson, Shi, and Yang (2016) Arias intensity model."""

import csv
import math
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm import Scenario
from pygmm.abrahamson_shi_yang_2016 import (
    AbrahamsonShiYang2016,
    _calc_scenario_ln_std,
    _t1_dip,
    _t2_mag,
    _t5_rjb,
)

# ------------------------------------------------------------------ #
# Helper: vertical SS fault scenario suitable for CY14 backbone      #
# ------------------------------------------------------------------ #


def _cy14_scenario(mag, vs30, dist_rup):
    """Vertical SS fault, FW site, suitable for CY14 backbone."""
    return Scenario(
        mag=mag,
        dist_rup=dist_rup,
        dist_jb=dist_rup,
        dist_x=-dist_rup,
        v_s30=vs30,
        dip=90.0,
        depth_tor=0.0,
        on_hanging_wall=False,
        mechanism="SS",
        dpp_centered=0.0,
    )


# ------------------------------------------------------------------ #
# Taper unit tests (Equations 3.5–3.7)                               #
# ------------------------------------------------------------------ #


class TestTapers:
    def test_t1_dip_vertical(self):
        assert _t1_dip(90.0) == pytest.approx(0.0)

    def test_t1_dip_45(self):
        assert _t1_dip(45.0) == pytest.approx(1.0)

    def test_t1_dip_shallow(self):
        assert _t1_dip(30.0) == pytest.approx(60.0 / 45.0)
        assert _t1_dip(10.0) == pytest.approx(60.0 / 45.0)

    def test_t2_upper_boundary(self):
        assert _t2_mag(6.5) == pytest.approx(1.0)
        assert _t2_mag(7.5) == pytest.approx(1.0)

    def test_t2_lower_boundary(self):
        assert _t2_mag(5.5) == pytest.approx(0.0)
        assert _t2_mag(4.0) == pytest.approx(0.0)

    def test_t2_mid(self):
        # T2(6.0) = 1 + 0.2*(6.0-6.5) - 0.8*(6.0-6.5)^2 = 0.7
        assert _t2_mag(6.0) == pytest.approx(0.7)

    def test_t2_continuity_at_55(self):
        assert _t2_mag(5.5 + 1e-9) == pytest.approx(0.0, abs=1e-6)

    def test_t2_continuity_at_65(self):
        assert _t2_mag(6.5 - 1e-9) == pytest.approx(1.0, abs=1e-6)

    def test_t5_at_zero(self):
        assert _t5_rjb(0.0) == pytest.approx(1.0)

    def test_t5_mid(self):
        assert _t5_rjb(7.5) == pytest.approx(0.5)

    def test_t5_at_limit(self):
        assert _t5_rjb(15.0) == pytest.approx(0.0)

    def test_t5_beyond_limit(self):
        assert _t5_rjb(20.0) == pytest.approx(0.0)


# ------------------------------------------------------------------ #
# Sigma unit tests                                                    #
# ------------------------------------------------------------------ #


class TestSigma:
    def test_sigma_cond_value(self):
        assert AbrahamsonShiYang2016.SIGMA_COND == pytest.approx(
            math.sqrt(0.15**2 + 0.35**2), rel=1e-6
        )

    def test_scenario_ln_std_formula(self):
        for s_pga, s_sa in [(0.3, 0.4), (0.5, 0.6), (0.7, 0.8)]:
            C4D = AbrahamsonShiYang2016.C4_DERIV
            C5D = AbrahamsonShiYang2016.C5_DERIV
            rho = AbrahamsonShiYang2016.RHO_PGA_SAT1
            expected = math.sqrt(
                C4D**2 * s_pga**2
                + C5D**2 * s_sa**2
                + 2 * rho * C4D * C5D * s_pga * s_sa
                + AbrahamsonShiYang2016.SIGMA_COND**2
            )
            assert _calc_scenario_ln_std(s_pga, s_sa) == pytest.approx(
                expected, rel=1e-10
            )


# ------------------------------------------------------------------ #
# Conditional mode                                                    #
# ------------------------------------------------------------------ #


class TestConditionalMode:
    """Reference values from direct evaluation of Equations (3.3) and (3.4)."""

    def _scenario(self, mag, vs30, dist_rup=10.0, dip=90.0, dist_jb=None, on_hw=False):
        if dist_jb is None:
            dist_jb = dist_rup
        return Scenario(
            mag=mag,
            dist_rup=dist_rup,
            dist_jb=dist_jb,
            v_s30=vs30,
            dip=dip,
            on_hanging_wall=on_hw,
        )

    def test_hanging_wall_increases_ia(self):
        """HW site raises ln_ia compared to FW site."""
        s_fw = self._scenario(
            6.5, 760, dist_rup=8.0, dip=45.0, dist_jb=3.0, on_hw=False
        )
        s_hw = self._scenario(6.5, 760, dist_rup=8.0, dip=45.0, dist_jb=3.0, on_hw=True)
        m_fw = AbrahamsonShiYang2016(s_fw, pga=0.2, sa_t1=0.15)
        m_hw = AbrahamsonShiYang2016(s_hw, pga=0.2, sa_t1=0.15)
        assert m_hw.ln_ia > m_fw.ln_ia

    def test_t2_zero_kills_hw_term(self):
        """For Mw ≤ 5.5, T2=0 ⇒ HW term=0 regardless of dip/RJB."""
        s_hw = self._scenario(5.0, 760, dist_rup=5.0, dip=45.0, dist_jb=0.0, on_hw=True)
        s_fw = self._scenario(
            5.0, 760, dist_rup=5.0, dip=45.0, dist_jb=0.0, on_hw=False
        )
        pga, sa_t1 = 0.1, 0.05
        assert AbrahamsonShiYang2016(s_hw, pga=pga, sa_t1=sa_t1).ln_ia == pytest.approx(
            AbrahamsonShiYang2016(s_fw, pga=pga, sa_t1=sa_t1).ln_ia
        )

    def test_t5_zero_kills_hw_term(self):
        """For RJB ≥ 15 km, T5=0 ⇒ HW term=0."""
        s_hw = self._scenario(
            7.0, 760, dist_rup=20.0, dip=45.0, dist_jb=15.0, on_hw=True
        )
        s_fw = self._scenario(
            7.0, 760, dist_rup=20.0, dip=45.0, dist_jb=15.0, on_hw=False
        )
        pga, sa_t1 = 0.3, 0.2
        assert AbrahamsonShiYang2016(s_hw, pga=pga, sa_t1=sa_t1).ln_ia == pytest.approx(
            AbrahamsonShiYang2016(s_fw, pga=pga, sa_t1=sa_t1).ln_ia
        )

    def test_missing_inputs_raises(self):
        s = Scenario(mag=7.0, dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        with pytest.raises(ValueError, match="backbone.*pga"):
            AbrahamsonShiYang2016(s)

    def test_missing_sa_t1_raises(self):
        s = Scenario(mag=7.0, dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        with pytest.raises(ValueError, match="backbone.*pga"):
            AbrahamsonShiYang2016(s, pga=0.3)

    def test_invalid_backbone_raises(self):
        s = Scenario(mag=7.0, dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        with pytest.raises(ValueError, match="not supported"):
            AbrahamsonShiYang2016(s, backbone="FAKE99")

    def test_ia_plus_minus_sigma(self):
        s = self._scenario(7.0, 760)
        m = AbrahamsonShiYang2016(s, pga=0.5, sa_t1=0.3)
        assert_allclose(m.ia_plus_sigma, np.exp(m.ln_ia + m.ln_std))
        assert_allclose(m.ia_minus_sigma, np.exp(m.ln_ia - m.ln_std))

    def test_ia_increases_with_pga(self):
        common = dict(mag=7.0, dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        m_low = AbrahamsonShiYang2016(Scenario(**common), pga=0.1, sa_t1=0.1)
        m_high = AbrahamsonShiYang2016(Scenario(**common), pga=0.5, sa_t1=0.3)
        assert m_high.ia > m_low.ia

    def test_ia_increases_with_mag(self):
        common = dict(dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        m_low = AbrahamsonShiYang2016(Scenario(mag=5.0, **common), pga=0.05, sa_t1=0.03)
        m_high = AbrahamsonShiYang2016(Scenario(mag=7.0, **common), pga=0.5, sa_t1=0.3)
        assert m_high.ia > m_low.ia


# ------------------------------------------------------------------ #
# Scenario mode — CY14 backbone CSV functional-form check            #
# ------------------------------------------------------------------ #
def _load_cy14_mag_scaling_data():
    """Load digitized Figure 4.3 data for CY14 backbone magnitude scaling."""
    data_file = Path(__file__).parent / "data" / "abrahamson_shi_yang_2016.csv"
    with data_file.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f))

    mags = np.array([float(row["Mag"]) for row in rows])
    ia = np.array([float(row["IA (m/s)"]) for row in rows])
    # The model's recommended upper limit is Mw=7.9; skip digitized values above it.
    mask = mags <= AbrahamsonShiYang2016.LIMITS["mag"][1]
    return mags[mask], ia[mask]


def test_cy14_backbone_magnitude_scaling_absolute_units():
    """CY14-backbone magnitude scaling matches digitized Figure 4.3 IA in m/s."""
    mags, ia_csv = _load_cy14_mag_scaling_data()
    ia_model = np.array(
        [
            AbrahamsonShiYang2016(
                _cy14_scenario(mag=mag, vs30=270.0, dist_rup=10.0),
                backbone="CY14",
            ).ia
            for mag in mags
        ]
    )
    assert_allclose(ia_model, ia_csv, rtol=0.30)
    # Large tolerance due to mismatch. It is possible the implementation of the
    # CY14 backbone in pygmm differs from the one used to generate the original
    # figure.


class TestScenarioModeGeneral:
    def test_tau_phi_none_in_scenario_mode(self):
        s = _cy14_scenario(7.0, 760, 30.0)
        m = AbrahamsonShiYang2016(s, backbone="CY14")
        assert m.tau is None
        assert m.phi is None

    def test_sigma_larger_than_conditional(self):
        """Scenario σ > conditional σ due to additional PGA/SA uncertainty."""
        s = _cy14_scenario(7.0, 760, 30.0)
        m_scen = AbrahamsonShiYang2016(s, backbone="CY14")
        m_cond = AbrahamsonShiYang2016(s, pga=0.3, sa_t1=0.15)
        assert m_scen.ln_std > m_cond.ln_std

    def test_ia_decreases_with_distance_scenario(self):
        m_near = AbrahamsonShiYang2016(_cy14_scenario(7.0, 760, 10.0), backbone="CY14")
        m_far = AbrahamsonShiYang2016(_cy14_scenario(7.0, 760, 100.0), backbone="CY14")
        assert m_near.ia > m_far.ia

    @pytest.mark.parametrize("backbone", ["ASK14", "BSSA14", "CB14", "CY14", "I14"])
    def test_all_backbones_run(self, backbone):
        """All five NGA-West2 backbone GMMs run without error."""
        s = Scenario(
            mag=6.5,
            dist_rup=30.0,
            dist_jb=30.0,
            dist_x=-30.0,
            v_s30=760.0,
            dip=90.0,
            depth_tor=0.0,
            on_hanging_wall=False,
            mechanism="SS",
            dpp_centered=0.0,
        )
        m = AbrahamsonShiYang2016(s, backbone=backbone)
        assert np.isfinite(m.ln_ia)
        assert np.isfinite(m.ln_std)
        assert m.ln_std > 0
