"""Smoke tests for pygmm.contracts dataclasses.

These exist primarily to (1) guarantee the module imports cleanly and
(2) lock in field names so they cannot drift away from the consumer-side
duplicates in pyrvt and pystrata without breaking a test.
"""

from dataclasses import fields

import numpy as np
import pytest

from pygmm import contracts


def _names(cls):
    return {f.name for f in fields(cls)}


def test_response_spectrum_fields():
    assert _names(contracts.ResponseSpectrum) == {
        "periods",
        "spec_accels",
        "damping",
        "duration",
    }


def test_fourier_spectrum_fields():
    assert _names(contracts.FourierSpectrum) == {"freqs", "fourier_amps", "duration"}


def test_duration_fields():
    assert _names(contracts.Duration) == {"duration", "d_5_75", "d_5_95", "d_20_80"}


def test_nonlinear_soil_curves_fields():
    assert _names(contracts.NonlinearSoilCurves) == {
        "strains",
        "mod_reduc",
        "damping",
        "damping_min",
        "unit_wt",
        "name",
    }


def test_velocity_profile_fields():
    assert _names(contracts.VelocityProfile) == {
        "depth",
        "vs_median",
        "std_vs_ln",
        "region",
        "site_class",
    }


def test_fault_displacement_fields():
    assert _names(contracts.FaultDisplacement) == {
        "mag",
        "displacement_mean",
        "sigma_ln_disp",
        "dist_from_rupture",
        "position_ratio",
    }


def test_cpt_sounding_fields():
    assert _names(contracts.CptSounding) == {
        "depth",
        "q_c",
        "f_s",
        "u_2",
        "water_table_depth",
        "unit_wts",
    }


def test_soil_behavior_profile_fields():
    assert _names(contracts.SoilBehaviorProfile) == {
        "depth",
        "ic",
        "sbt_class",
        "fines_content",
    }


def test_liquefaction_triggering_fields():
    assert _names(contracts.LiquefactionTriggering) == {
        "depth",
        "csr",
        "crr",
        "factor_of_safety",
        "prob_liquefaction",
    }


def test_dataclasses_are_frozen():
    rs = contracts.ResponseSpectrum(
        periods=np.array([0.1, 0.2]),
        spec_accels=np.array([0.5, 0.7]),
        damping=0.05,
    )
    with pytest.raises(Exception):
        rs.damping = 0.10  # type: ignore[misc]


def test_response_spectrum_construct():
    rs = contracts.ResponseSpectrum(
        periods=np.array([0.01, 0.1, 1.0]),
        spec_accels=np.array([0.4, 0.9, 0.2]),
        damping=0.05,
        duration=12.5,
    )
    assert rs.duration == 12.5
    assert rs.spec_accels.shape == (3,)


def test_fourier_spectrum_construct():
    fs = contracts.FourierSpectrum(
        freqs=np.logspace(-1, 2, 50),
        fourier_amps=np.ones(50),
        duration=8.0,
    )
    assert fs.fourier_amps.shape == fs.freqs.shape
