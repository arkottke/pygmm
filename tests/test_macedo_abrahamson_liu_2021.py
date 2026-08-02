"""Tests for Macedo, Abrahamson, and Liu (2021) CAV model."""

import csv
from pathlib import Path

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pygmm import Scenario
from pygmm.gm_intensity.macedo_abrahamson_liu_2021 import (
    MacedoAbrahamsonLiu2021,
    _calc_scenario_ln_std,
    _t1_dip,
    _t2_mag,
    _t5_rjb,
)

# ------------------------------------------------------------------ #
# Helper: vertical SS fault scenario for Figure 9 comparisons        #
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
# Taper unit tests (Eqs. 13a–c)                                      #
# ------------------------------------------------------------------ #


class TestTapers:
    def test_t1_dip_vertical(self):
        assert _t1_dip(90.0) == pytest.approx(0.0)

    def test_t1_dip_45(self):
        assert _t1_dip(45.0) == pytest.approx(1.0)

    def test_t1_dip_shallow(self):
        # dip <= 30 → T1 = 60/45
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
# Sigma functions                                                     #
# ------------------------------------------------------------------ #


class TestSigma:
    def test_sigma_cond_value(self):
        # sqrt(0.17^2 + 0.26^2)
        assert MacedoAbrahamsonLiu2021.SIGMA_COND == pytest.approx(
            np.sqrt(0.17**2 + 0.26**2), rel=1e-6
        )

    def test_scenario_ln_std_formula(self):
        # Eq. (17): σ = sqrt(0.097 + 0.4541 * σ²_PGA)
        for s in [0.3, 0.5, 0.7]:
            assert _calc_scenario_ln_std(s) == pytest.approx(
                np.sqrt(0.097 + 0.4541 * s**2), rel=1e-10
            )


# ------------------------------------------------------------------ #
# Conditional mode                                                    #
# ------------------------------------------------------------------ #


class TestConditionalMode:
    """Reference values from direct evaluation of Eq. (12)."""

    def test_basic_mw7_r10(self):
        """Mw=7.0, Rrup=10 km, Vs30=760 m/s, PGA=0.5 g."""
        s = Scenario(
            mag=7.0,
            dist_rup=10.0,
            dist_jb=10.0,
            v_s30=760.0,
            dip=90.0,
            on_hanging_wall=False,
        )
        m = MacedoAbrahamsonLiu2021(s, pga=0.5)
        assert_allclose(m.ln_cav, 2.19194500, rtol=1e-5)
        assert_allclose(m.cav, 8.95260905, rtol=1e-5)
        assert_allclose(m.ln_std, 0.31064449, rtol=1e-5)
        assert m.tau == pytest.approx(0.17)
        assert m.phi == pytest.approx(0.26)

    def test_basic_mw5_r30(self):
        """Mw=5.0, Rrup=30 km, Vs30=760 m/s, PGA=0.05 g."""
        s = Scenario(
            mag=5.0,
            dist_rup=30.0,
            dist_jb=30.0,
            v_s30=760.0,
            dip=90.0,
            on_hanging_wall=False,
        )
        m = MacedoAbrahamsonLiu2021(s, pga=0.05)
        assert_allclose(m.ln_cav, -0.49364340, rtol=1e-5)
        assert_allclose(m.cav, 0.61039841, rtol=1e-5)

    def test_hanging_wall_on(self):
        """Hanging-wall term is active for dip=45, Mw=6.0, RJB=5 km."""
        s = Scenario(
            mag=6.0,
            dist_rup=8.0,
            dist_jb=5.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=True,
        )
        m = MacedoAbrahamsonLiu2021(s, pga=0.2)
        assert_allclose(m.ln_cav, 1.08794372, rtol=1e-5)

    def test_hanging_wall_off_equals_fw(self):
        """F_HW=0 (FW site) must zero the c6 term."""
        s_fw = Scenario(
            mag=6.0,
            dist_rup=8.0,
            dist_jb=5.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=False,
        )
        s_hw = Scenario(
            mag=6.0,
            dist_rup=8.0,
            dist_jb=5.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=True,
        )
        m_fw = MacedoAbrahamsonLiu2021(s_fw, pga=0.2)
        m_hw = MacedoAbrahamsonLiu2021(s_hw, pga=0.2)
        assert m_hw.ln_cav > m_fw.ln_cav  # HW should increase CAV

    def test_t2_zero_kills_hw_term(self):
        """For Mw ≤ 5.5, T2=0 ⇒ HW term=0 regardless of dip/RJB."""
        s_hw = Scenario(
            mag=5.0,
            dist_rup=5.0,
            dist_jb=0.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=True,
        )
        s_fw = Scenario(
            mag=5.0,
            dist_rup=5.0,
            dist_jb=0.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=False,
        )
        pga = 0.1
        assert MacedoAbrahamsonLiu2021(s_hw, pga=pga).ln_cav == pytest.approx(
            MacedoAbrahamsonLiu2021(s_fw, pga=pga).ln_cav
        )

    def test_t5_zero_kills_hw_term(self):
        """For RJB ≥ 15 km, T5=0 ⇒ HW term=0."""
        s_hw = Scenario(
            mag=7.0,
            dist_rup=20.0,
            dist_jb=15.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=True,
        )
        s_fw = Scenario(
            mag=7.0,
            dist_rup=20.0,
            dist_jb=15.0,
            v_s30=760.0,
            dip=45.0,
            on_hanging_wall=False,
        )
        pga = 0.3
        assert MacedoAbrahamsonLiu2021(s_hw, pga=pga).ln_cav == pytest.approx(
            MacedoAbrahamsonLiu2021(s_fw, pga=pga).ln_cav
        )

    def test_pga_and_pga_model_both_none_raises(self):
        s = Scenario(mag=7.0, dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        with pytest.raises(ValueError, match="pga_model.*pga"):
            MacedoAbrahamsonLiu2021(s)

    def test_invalid_pga_model_raises(self):
        s = Scenario(mag=7.0, dist_rup=10.0, dist_jb=10.0, v_s30=760.0)
        with pytest.raises(ValueError, match="not supported"):
            MacedoAbrahamsonLiu2021(s, pga_model="FAKE99")

    def test_cav_plus_minus_sigma(self):
        s = Scenario(
            mag=7.0,
            dist_rup=10.0,
            dist_jb=10.0,
            v_s30=760.0,
            dip=90.0,
            on_hanging_wall=False,
        )
        m = MacedoAbrahamsonLiu2021(s, pga=0.5)
        assert_allclose(m.cav_plus_sigma, np.exp(m.ln_cav + m.ln_std))
        assert_allclose(m.cav_minus_sigma, np.exp(m.ln_cav - m.ln_std))

    def test_cav_increases_with_mag(self):
        common = dict(dist_rup=20.0, dist_jb=20.0, v_s30=760.0)
        m5 = MacedoAbrahamsonLiu2021(Scenario(mag=5.0, **common), pga=0.2)
        m7 = MacedoAbrahamsonLiu2021(Scenario(mag=7.0, **common), pga=0.2)
        assert m7.cav > m5.cav

    def test_cav_decreases_with_distance(self):
        common = dict(mag=6.5, v_s30=760.0)
        mn = MacedoAbrahamsonLiu2021(
            Scenario(dist_rup=5.0, dist_jb=5.0, **common), pga=0.4
        )
        mf = MacedoAbrahamsonLiu2021(
            Scenario(dist_rup=100.0, dist_jb=100.0, **common), pga=0.05
        )
        assert mn.cav > mf.cav


# ------------------------------------------------------------------ #
# Scenario mode — CY14 backbone, Figure 9 CSV functional-form check  #
# ------------------------------------------------------------------ #


def _load_cy14_distance_scaling_data():
    """Load digitized Figure 9 data (Mw=4.75, Vs30=760 m/s)."""
    data_file = Path(__file__).parent / "data" / "macedo_abrahamson_liu_2021.csv"
    with data_file.open(encoding="utf-8-sig", newline="") as f:
        rows = [{k.strip(): v for k, v in row.items()} for row in csv.DictReader(f)]

    rrup = np.array([float(row["Rrup"]) for row in rows])
    cav = np.array([float(row["CAV"]) for row in rows])
    return rrup, cav


def test_cy14_scenario_distance_scaling_from_csv():
    """CY14 scenario mode follows Figure 9 distance-scaling shape from CSV."""
    rrup, cav_csv = _load_cy14_distance_scaling_data()
    ln_cav_model = np.array(
        [
            MacedoAbrahamsonLiu2021(
                _cy14_scenario(mag=4.75, vs30=760.0, dist_rup=dist_rup),
                pga_model="CY14",
            ).ln_cav
            for dist_rup in rrup
        ]
    )

    # CSV comes from figure digitization; compare normalized shape in log-space.
    # This preserves the distance-scaling trend while remaining robust to axis scaling.
    csv_norm = (np.log(cav_csv) - np.log(cav_csv[0])) / (
        np.log(cav_csv[-1]) - np.log(cav_csv[0])
    )
    model_norm = (ln_cav_model - ln_cav_model[0]) / (ln_cav_model[-1] - ln_cav_model[0])
    assert_allclose(model_norm, csv_norm, atol=0.015)


class TestScenarioModeGeneral:
    def test_tau_phi_none_in_scenario_mode(self):
        """tau and phi are not separately tracked in scenario mode."""
        s = _cy14_scenario(7.0, 760, 30.0)
        m = MacedoAbrahamsonLiu2021(s, pga_model="CY14")
        assert m.tau is None
        assert m.phi is None

    def test_sigma_lower_than_conditional(self):
        """Scenario σ > conditional σ because PGA uncertainty is additive."""
        s = _cy14_scenario(7.0, 760, 30.0)
        m_scen = MacedoAbrahamsonLiu2021(s, pga_model="CY14")
        m_cond = MacedoAbrahamsonLiu2021(s, pga=0.3)
        assert m_scen.ln_std > m_cond.ln_std

    def test_cav_increases_with_mag_scenario(self):
        """Median CAV increases with magnitude in scenario mode."""
        m_low = MacedoAbrahamsonLiu2021(
            _cy14_scenario(5.75, 760, 30.0), pga_model="CY14"
        )
        m_high = MacedoAbrahamsonLiu2021(
            _cy14_scenario(7.6, 760, 30.0), pga_model="CY14"
        )
        assert m_high.cav > m_low.cav

    def test_cav_decreases_with_distance_scenario(self):
        """Median CAV decreases with rupture distance in scenario mode."""
        m_near = MacedoAbrahamsonLiu2021(
            _cy14_scenario(7.0, 760, 10.0), pga_model="CY14"
        )
        m_far = MacedoAbrahamsonLiu2021(
            _cy14_scenario(7.0, 760, 100.0), pga_model="CY14"
        )
        assert m_near.cav > m_far.cav
