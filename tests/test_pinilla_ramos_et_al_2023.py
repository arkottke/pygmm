#!/usr/bin/env python
"""Test cases for Pinilla-Ramos et al. (2023) duration model."""

import gzip
import json
import os

import pytest
from numpy.testing import assert_allclose

import pygmm
from pygmm.pinilla_ramos_et_al_2023 import PinillaRamosEtAl2023

# Load test data
fpath = os.path.join(
    os.path.dirname(__file__), "data", "pinilla_ramos_et_al_2023.json.gz"
)
with gzip.open(fpath, "rt", encoding="utf-8") as fp:
    data = json.load(fp)


def idfn(case):
    """Generate test case ID from description."""
    if isinstance(case, dict):
        return case.get("description", "Unknown case")
    else:
        return str(case)


@pytest.mark.parametrize(
    "case",
    [
        case
        for case in data["test_cases"]
        if not isinstance(case["inputs"]["mag"], list)
    ],
    ids=idfn,
)
def test_model_class(case):
    """Test the PinillaRamosEtAl2023 class against expected values."""
    inputs = case["inputs"]
    expected = case["expected"]

    # Create scenario and model instance
    scenario = pygmm.model.Scenario(
        mag=inputs["mag"], dist_rup=inputs["rrup"], v_s30=inputs["vs30"]
    )
    model_inst = PinillaRamosEtAl2023(scenario)

    if inputs["energy"] == 0.75:
        # Test default D5-75 properties
        duration_val = model_inst.duration
        plus_val = model_inst.duration_plus_sigma
        minus_val = model_inst.duration_minus_sigma
    else:
        # Test custom energy threshold
        duration_val, plus_val, minus_val = model_inst.duration_for_energy(
            inputs["energy"]
        )

    # Extract expected values
    expected_duration = (
        expected["duration"][0]
        if isinstance(expected["duration"], list)
        else expected["duration"]
    )
    expected_plus = (
        expected["duration_plus_sigma"][0]
        if isinstance(expected["duration_plus_sigma"], list)
        else expected["duration_plus_sigma"]
    )
    expected_minus = (
        expected["duration_minus_sigma"][0]
        if isinstance(expected["duration_minus_sigma"], list)
        else expected["duration_minus_sigma"]
    )

    assert_allclose(duration_val, expected_duration, rtol=1e-6, atol=1e-8)
    assert_allclose(plus_val, expected_plus, rtol=1e-6, atol=1e-8)
    assert_allclose(minus_val, expected_minus, rtol=1e-6, atol=1e-8)


def test_model_properties():
    """Test model properties and metadata."""
    scenario = pygmm.model.Scenario(mag=6.5, dist_rup=20.0, v_s30=300.0)
    model_inst = PinillaRamosEtAl2023(scenario)

    # Test model metadata
    assert model_inst.NAME == "Pinilla-Ramos et al. (2023)"
    assert model_inst.ABBREV == "PR23"

    # Test that properties return reasonable values
    assert isinstance(model_inst.duration, float)
    assert isinstance(model_inst.duration_plus_sigma, float)
    assert isinstance(model_inst.duration_minus_sigma, float)

    # Test that plus sigma > duration > minus sigma
    assert model_inst.duration_plus_sigma > model_inst.duration
    assert model_inst.duration > model_inst.duration_minus_sigma

    # Test that all values are positive
    assert model_inst.duration > 0
    assert model_inst.duration_plus_sigma > 0
    assert model_inst.duration_minus_sigma > 0


def test_energy_threshold_validation():
    """Test that energy threshold validation works correctly."""
    scenario = pygmm.model.Scenario(mag=6.5, dist_rup=20.0, v_s30=300.0)
    model_inst = PinillaRamosEtAl2023(scenario)

    # Test valid energy thresholds
    for energy in [0.10, 0.25, 0.45, 0.75, 0.95]:
        duration, plus_sigma, minus_sigma = model_inst.duration_for_energy(energy)
        assert isinstance(duration, float)
        assert isinstance(plus_sigma, float)
        assert isinstance(minus_sigma, float)
        assert duration > 0
        assert plus_sigma > duration
        assert duration > minus_sigma

    # Test invalid energy thresholds
    with pytest.raises(ValueError, match="Energy threshold .* is outside valid range"):
        model_inst.duration_for_energy(0.05)

    with pytest.raises(ValueError, match="Energy threshold .* is outside valid range"):
        model_inst.duration_for_energy(0.99)


def test_parameter_bounds():
    """Test model behavior at parameter boundaries."""
    # Test minimum magnitude
    scenario = pygmm.model.Scenario(mag=4.0, dist_rup=10.0, v_s30=300.0)
    model_inst = PinillaRamosEtAl2023(scenario)
    assert model_inst.duration > 0

    # Test maximum magnitude
    scenario = pygmm.model.Scenario(mag=8.0, dist_rup=10.0, v_s30=300.0)
    model_inst = PinillaRamosEtAl2023(scenario)
    assert model_inst.duration > 0

    # Test minimum vs30
    scenario = pygmm.model.Scenario(mag=6.5, dist_rup=20.0, v_s30=150.0)
    model_inst = PinillaRamosEtAl2023(scenario)
    assert model_inst.duration > 0

    # Test maximum vs30
    scenario = pygmm.model.Scenario(mag=6.5, dist_rup=20.0, v_s30=2000.0)
    model_inst = PinillaRamosEtAl2023(scenario)
    assert model_inst.duration > 0


def test_vs30_scaling():
    """Test that duration decreases with increasing Vs30 (generally)."""
    mag, rrup = 6.5, 20.0

    # Create scenarios with different Vs30 values
    vs30_values = [200, 400, 800, 1500]
    durations = []

    for vs30 in vs30_values:
        scenario = pygmm.model.Scenario(mag=mag, dist_rup=rrup, v_s30=vs30)
        model_inst = PinillaRamosEtAl2023(scenario)
        durations.append(model_inst.duration)

    # Soft soil should generally have longer durations than hard rock
    assert durations[0] > durations[-1]  # Vs30=200 > Vs30=1500


def test_magnitude_scaling():
    """Test that duration increases with magnitude."""
    rrup, vs30 = 20.0, 300.0

    # Create scenarios with different magnitudes
    mag_values = [5.0, 6.0, 7.0, 8.0]
    durations = []

    for mag in mag_values:
        scenario = pygmm.model.Scenario(mag=mag, dist_rup=rrup, v_s30=vs30)
        model_inst = PinillaRamosEtAl2023(scenario)
        durations.append(model_inst.duration)

    # Duration should generally increase with magnitude
    assert durations[-1] > durations[0]  # M8.0 > M5.0
    assert durations[2] > durations[1]  # M7.0 > M6.0


def test_distance_scaling():
    """Test duration behavior with distance."""
    mag, vs30 = 6.5, 300.0

    # Create scenarios with different distances
    rrup_values = [5.0, 20.0, 50.0, 150.0]
    durations = []

    for rrup in rrup_values:
        scenario = pygmm.model.Scenario(mag=mag, dist_rup=rrup, v_s30=vs30)
        model_inst = PinillaRamosEtAl2023(scenario)
        durations.append(model_inst.duration)

    # Duration should generally increase with distance
    assert durations[-1] > durations[0]  # 150km > 5km


def test_energy_threshold_behavior():
    """Test that different energy thresholds produce different durations."""
    scenario = pygmm.model.Scenario(mag=6.5, dist_rup=20.0, v_s30=300.0)
    model_inst = PinillaRamosEtAl2023(scenario)

    # Get durations for different energy thresholds
    durations = {}
    for energy in [0.25, 0.45, 0.75, 0.95]:
        duration, _, _ = model_inst.duration_for_energy(energy)
        durations[energy] = duration

    # Duration should increase with higher energy percentiles
    # D5-95 (5% to 95% energy) should be longer than D5-75 (5% to 75% energy), etc.
    assert durations[0.95] > durations[0.75]
    assert durations[0.75] > durations[0.45]
    assert durations[0.45] > durations[0.25]


if __name__ == "__main__":
    pytest.main([__file__])
