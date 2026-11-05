#!/usr/bin/env python3
"""Tests for Pinilla-Ramos et al. (2024) subduction duration model."""

import functools
import itertools
import os
import warnings

import numpy as np
import pandas as pd
import pytest

from pygmm import Scenario
from pygmm.ground_motion.pinilla_ramos_et_al_2024 import (
    PinillaRamosEtAl2024,
    duration_model,
)

# Vectorized and scalar results can differ in the last bit, because NumPy's
# SIMD loops may round differently than its scalar functions on some CPUs
assert_close = functools.partial(np.testing.assert_allclose, rtol=1e-12, atol=0)


class TestPinillaRamosEtAl2024:
    """Test suite for Pinilla-Ramos et al. (2024) duration model."""

    def test_model_initialization(self):
        """Test basic model initialization."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        assert model.NAME == "Pinilla-Ramos et al. (2024)"
        assert model.ABBREV == "PR24"
        assert hasattr(model, "duration")
        assert hasattr(model, "duration_plus_sigma")
        assert hasattr(model, "duration_minus_sigma")

    def test_invalid_region_event_combinations(self):
        """Test validation of region-event type combinations."""
        # Interface earthquakes not supported in Taiwan
        with pytest.raises(
            ValueError, match="No model available for interface earthquakes in Taiwan"
        ):
            scenario = Scenario(
                mag=7.0,
                dist_rup=100.0,
                v_s30=600.0,
                region="Taiwan",
                event_type="interface",
            )
            PinillaRamosEtAl2024(scenario)

    def test_invalid_event_type(self):
        """Test validation of event type."""
        with pytest.raises(
            ValueError, match="event_type must be 'interface' or 'slab'"
        ):
            scenario = Scenario(
                mag=7.0,
                dist_rup=100.0,
                v_s30=600.0,
                region="Japan",
                event_type="crustal",  # Invalid event type
            )
            PinillaRamosEtAl2024(scenario)

    def test_invalid_region(self):
        """Test validation of region."""
        with pytest.raises(ValueError, match="region must be one of"):
            scenario = Scenario(
                mag=7.0,
                dist_rup=100.0,
                v_s30=600.0,
                region="California",  # Invalid region
                event_type="interface",
            )
            PinillaRamosEtAl2024(scenario)

    def test_interface_earthquake_japan(self):
        """Test interface earthquake in Japan."""
        scenario = Scenario(
            mag=7.5, dist_rup=50.0, v_s30=800.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        # Test D5-75 duration
        duration = model.duration
        assert isinstance(duration, float)
        assert duration > 0

        # Test properties
        assert model.d575_median == duration
        assert model.d575_sigma > 0
        assert model.duration_plus_sigma > duration
        assert model.duration_minus_sigma < duration

    def test_slab_earthquake_taiwan(self):
        """Test slab earthquake in Taiwan."""
        scenario = Scenario(
            mag=6.5, dist_rup=25.0, v_s30=600.0, region="Taiwan", event_type="slab"
        )
        model = PinillaRamosEtAl2024(scenario)

        duration = model.duration
        assert isinstance(duration, float)
        assert duration > 0

    def test_energy_threshold_methods(self):
        """Test energy threshold calculation methods."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        # Test valid energy thresholds
        d5_10 = model.d5x_median("D5-10")
        d5_75 = model.d5x_median("D5-75")
        d5_95 = model.d5x_median("D5-95")

        # D5-10 should be shorter than D5-75, which should be shorter than D5-95
        assert d5_10 < d5_75 < d5_95

        # Test corresponding sigmas
        sigma_10 = model.d5x_sigma("D5-10")
        sigma_75 = model.d5x_sigma("D5-75")
        sigma_95 = model.d5x_sigma("D5-95")

        assert all(s > 0 for s in [sigma_10, sigma_75, sigma_95])

    def test_invalid_energy_threshold(self):
        """Test invalid energy threshold."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        with pytest.raises(ValueError, match="Unsupported energy threshold"):
            model.d5x_median("D5-99")

    def test_duration_for_energy_method(self):
        """Test duration_for_energy method."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        # Test valid energy values
        median_95, plus_95, minus_95 = model.duration_for_energy(0.95)
        assert median_95 > 0
        assert plus_95 > median_95
        assert minus_95 < median_95

        # Test D5-75 comparison
        median_75, plus_75, minus_75 = model.duration_for_energy(0.75)
        assert median_75 < median_95  # D5-75 should be shorter than D5-95

    def test_duration_for_energy_invalid_range(self):
        """Test duration_for_energy with invalid range."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        with pytest.raises(ValueError, match="outside valid range"):
            model.duration_for_energy(0.05)  # Too low

        with pytest.raises(ValueError, match="outside valid range"):
            model.duration_for_energy(0.98)  # Too high

    def test_magnitude_scaling(self):
        """Test magnitude scaling behavior."""
        base_scenario = Scenario(
            dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )

        durations = []
        for mag in [6.0, 7.0, 8.0]:
            scenario = base_scenario.copy_with(mag=mag)
            model = PinillaRamosEtAl2024(scenario)
            durations.append(model.duration)

        # Duration should generally increase with magnitude
        assert durations[1] > durations[0]  # M7 > M6
        assert durations[2] > durations[1]  # M8 > M7

    def test_distance_scaling(self):
        """Test distance scaling behavior."""
        base_scenario = Scenario(
            mag=7.0, v_s30=600.0, region="Japan", event_type="interface"
        )

        durations = []
        for dist in [10.0, 100.0, 200.0]:
            scenario = base_scenario.copy_with(dist_rup=dist)
            model = PinillaRamosEtAl2024(scenario)
            durations.append(model.duration)

        # Duration should increase with distance
        assert durations[1] > durations[0]  # 100km > 10km
        assert durations[2] > durations[1]  # 200km > 100km

    def test_vs30_scaling(self):
        """Test Vs30 scaling behavior."""
        base_scenario = Scenario(
            mag=7.0, dist_rup=100.0, region="Japan", event_type="interface"
        )

        durations = []
        for vs30 in [200.0, 600.0, 1200.0]:
            scenario = base_scenario.copy_with(v_s30=vs30)
            model = PinillaRamosEtAl2024(scenario)
            durations.append(model.duration)

        # Duration should decrease with increasing Vs30 (harder sites)
        assert durations[1] < durations[0]  # 600 m/s < 200 m/s
        assert durations[2] < durations[1]  # 1200 m/s < 600 m/s

    def test_regional_differences(self):
        """Test regional differences for interface earthquakes."""
        base_scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, event_type="interface"
        )

        durations = {}
        for region in ["Japan", "New Zealand", "South America"]:
            scenario = base_scenario.copy_with(region=region)
            model = PinillaRamosEtAl2024(scenario)
            durations[region] = model.duration

        # All durations should be positive and different
        assert all(d > 0 for d in durations.values())
        assert len(set(durations.values())) == 3  # All different

    def test_slab_regional_differences(self):
        """Test regional differences for slab earthquakes."""
        base_scenario = Scenario(mag=6.5, dist_rup=50.0, v_s30=600.0, event_type="slab")

        durations = {}
        for region in ["Japan", "New Zealand", "South America", "Taiwan"]:
            scenario = base_scenario.copy_with(region=region)
            model = PinillaRamosEtAl2024(scenario)
            durations[region] = model.duration

        # All durations should be positive and different
        assert all(d > 0 for d in durations.values())
        assert len(set(durations.values())) == 4  # All different

    def test_event_type_differences(self):
        """Test differences between interface and slab events."""
        base_scenario = Scenario(mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan")

        # Interface earthquake
        interface_scenario = base_scenario.copy_with(event_type="interface")
        interface_model = PinillaRamosEtAl2024(interface_scenario)

        # Slab earthquake
        slab_scenario = base_scenario.copy_with(event_type="slab")
        slab_model = PinillaRamosEtAl2024(slab_scenario)

        # Both should give positive durations
        assert interface_model.duration > 0
        assert slab_model.duration > 0

        # Durations should be different
        assert interface_model.duration != slab_model.duration

    def test_sigma_calculations(self):
        """Test standard deviation calculations."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        sigma = model.d575_sigma
        assert isinstance(sigma, float)
        assert sigma > 0
        assert sigma < 1.0  # Reasonable range for log-space sigma

    def test_reproducibility(self):
        """Test that the model gives reproducible results."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )

        model1 = PinillaRamosEtAl2024(scenario)
        model2 = PinillaRamosEtAl2024(scenario)

        assert model1.duration == model2.duration
        assert model1.d575_sigma == model2.d575_sigma

    def test_parameter_validation(self):
        """Test parameter validation through the model."""
        # Test that model handles edge cases of parameter ranges
        scenarios = [
            # Minimum values
            Scenario(
                mag=4.5,
                dist_rup=0.0,
                v_s30=150.0,
                region="Japan",
                event_type="interface",
            ),
            # Maximum values
            Scenario(
                mag=8.5,
                dist_rup=300.0,
                v_s30=2000.0,
                region="Japan",
                event_type="interface",
            ),
        ]

        for scenario in scenarios:
            model = PinillaRamosEtAl2024(scenario)
            assert model.duration > 0

    @pytest.mark.parametrize(
        "region,event_type",
        [
            ("Japan", "interface"),
            ("Japan", "slab"),
            ("New Zealand", "interface"),
            ("New Zealand", "slab"),
            ("South America", "interface"),
            ("South America", "slab"),
            ("Taiwan", "slab"),
        ],
    )
    def test_valid_region_event_combinations(self, region, event_type):
        """Test all valid region-event type combinations."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region=region, event_type=event_type
        )
        model = PinillaRamosEtAl2024(scenario)

        assert model.duration > 0
        assert model.d575_sigma > 0

    def test_against_test_data(self):
        """Test model against generated test data for consistency."""
        try:
            # Load test data if available
            fpath = os.path.join(
                os.path.dirname(__file__), "data", "pinilla_ramos_et_al_2024.csv.gz"
            )
            df = pd.read_csv(fpath)

            # Test a sample of cases
            sample_size = min(50, len(df))
            sample_df = df.sample(n=sample_size, random_state=42)

            for _, row in sample_df.iterrows():
                scenario = Scenario(
                    mag=row["mag"],
                    dist_rup=row["dist_rup"],
                    v_s30=row["v_s30"],
                    region=row["region"],
                    event_type=row["event_type"],
                )
                model = PinillaRamosEtAl2024(scenario)

                if row["energy_threshold"] == "D5-75":
                    computed_median = model.d575_median
                    computed_sigma = model.d575_sigma
                else:
                    computed_median = model.d5x_median(row["energy_threshold"])
                    computed_sigma = model.d5x_sigma(row["energy_threshold"])

                # Check that computed values match expected values within tolerance
                np.testing.assert_allclose(computed_median, row["median"], rtol=1e-10)
                np.testing.assert_allclose(computed_sigma, row["sigma"], rtol=1e-10)

        except FileNotFoundError:
            pytest.skip("Test data file not found")

    def test_energy_threshold_monotonicity(self):
        """Test that duration increases monotonically with energy threshold."""
        scenario = Scenario(
            mag=7.0, dist_rup=100.0, v_s30=600.0, region="Japan", event_type="interface"
        )
        model = PinillaRamosEtAl2024(scenario)

        thresholds = ["D5-10", "D5-30", "D5-50", "D5-75", "D5-90", "D5-95"]
        durations = [model.d5x_median(t) for t in thresholds]

        # Check monotonic increase
        for i in range(1, len(durations)):
            assert durations[i] > durations[i - 1], (
                f"Duration not increasing: {thresholds[i - 1]} to {thresholds[i]}"
            )


# Tests of the vectorized model against scalar scenarios

COMBOS = [
    ("Japan", "interface"),
    ("New Zealand", "interface"),
    ("South America", "interface"),
    ("Japan", "slab"),
    ("New Zealand", "slab"),
    ("South America", "slab"),
    ("Taiwan", "slab"),
]
# Slab magnitude scaling changes at 6.0, path scaling at 250 km
MAGS = [4.5, 5.99, 6.0, 6.01, 7.0, 8.5]
DISTS = [0.0, 100.0, 249.99, 250.0, 250.01, 300.0]
V_S30S = [150.0, 760.0, 2000.0]
# 0.333 is not a supported threshold and uses the closest (0.35)
ENERGIES = [0.1, 0.333, 0.5, 0.75, 0.95]
THRESHOLDS = ["D5-10", "D5-50", "D5-75", "D5-95"]
OUTPUTS = [
    "d575_median",
    "d575_sigma",
    "duration",
    "duration_plus_sigma",
    "duration_minus_sigma",
]


def scalar_results(region, event_type, mag, dist_rup, v_s30):
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            v_s30=v_s30,
            region=region,
            event_type=event_type,
        )
    )
    results = {key: getattr(m, key) for key in OUTPUTS}
    for threshold in THRESHOLDS:
        results[("d5x_median", threshold)] = m.d5x_median(threshold)
        results[("d5x_sigma", threshold)] = m.d5x_sigma(threshold)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for energy in ENERGIES:
            results[("duration_for_energy", energy)] = m.duration_for_energy(energy)
    return results


@pytest.fixture(scope="module")
def grid():
    rows = [
        (region, event_type, mag, dist, v_s30)
        for (region, event_type), mag, dist, v_s30 in itertools.product(
            COMBOS, MAGS, DISTS, V_S30S
        )
    ]
    region, event_type, mag, dist_rup, v_s30 = (np.array(c) for c in zip(*rows))
    scenario = Scenario(
        mag=mag, dist_rup=dist_rup, v_s30=v_s30, region=region, event_type=event_type
    )
    expected = [scalar_results(*row) for row in rows]
    return rows, scenario, expected


def vectorized_results(m):
    results = {key: getattr(m, key) for key in OUTPUTS}
    for threshold in THRESHOLDS:
        results[("d5x_median", threshold)] = m.d5x_median(threshold)
        results[("d5x_sigma", threshold)] = m.d5x_sigma(threshold)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for energy in ENERGIES:
            results[("duration_for_energy", energy)] = m.duration_for_energy(energy)
    return results


def test_vectorized_matches_scalar(grid):
    rows, scenario, expected = grid
    actual = vectorized_results(PinillaRamosEtAl2024(scenario))
    for key, value in actual.items():
        if key[0] == "duration_for_energy":
            assert isinstance(value, tuple) and len(value) == 3
            for i in range(3):
                assert_close(
                    value[i], np.array([e[key][i] for e in expected]), err_msg=str(key)
                )
        else:
            assert_close(value, np.array([e[key] for e in expected]), err_msg=str(key))


def test_vectorized_shapes(grid):
    rows, scenario, _ = grid
    m = PinillaRamosEtAl2024(scenario)
    for key, value in vectorized_results(m).items():
        if key[0] == "duration_for_energy":
            for v in value:
                assert v.shape == (len(rows),), key
        else:
            assert value.shape == (len(rows),), key


def test_scalar_types_are_unchanged():
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=7.0, dist_rup=100.0, v_s30=400.0, region="Japan", event_type="slab"
        )
    )
    for key, value in vectorized_results(m).items():
        if key[0] == "duration_for_energy":
            assert isinstance(value, tuple) and len(value) == 3
            assert all(type(v) is np.float64 for v in value), key
        else:
            assert type(value) is np.float64, key


def test_scalars_broadcast_with_arrays():
    mag = np.array([5.5, 6.5, 7.5])
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=mag,
            dist_rup=120.0,
            v_s30=500.0,
            region=np.array(["Japan", "Taiwan", "New Zealand"]),
            event_type="slab",
        )
    )
    for key, value in vectorized_results(m).items():
        values = value if key[0] == "duration_for_energy" else (value,)
        for v in values:
            assert v.shape == (3,), key
    for i, mag_i in enumerate(mag):
        expected = scalar_results(
            ["Japan", "Taiwan", "New Zealand"][i], "slab", mag_i, 120.0, 500.0
        )
        assert_close(m.duration[i], expected["duration"])
        assert_close(m.d575_sigma[i], expected["d575_sigma"])
        assert_close(
            m.duration_for_energy(0.95)[2][i],
            expected[("duration_for_energy", 0.95)][2],
        )


def test_values_depending_on_some_inputs_have_scenario_shape():
    # The standard deviation does not depend on the region, but it has the
    # broadcast shape of all of the scenario values
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=7.0,
            dist_rup=100.0,
            v_s30=400.0,
            region=np.array(["Japan", "Taiwan"]),
            event_type="slab",
        )
    )
    assert m.d575_sigma.shape == (2,)
    assert m.d5x_sigma("D5-95").shape == (2,)
    assert m.d575_sigma[0] == m.d575_sigma[1]
    # Broadcast results are writable arrays
    assert m.d575_sigma.flags.writeable


def test_multidimensional_inputs():
    mag = np.array([[5.5], [6.5], [7.5]])
    dist_rup = np.array([20.0, 150.0, 260.0, 300.0])
    event_type = np.array([["interface"], ["slab"], ["slab"]])
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=mag,
            dist_rup=dist_rup,
            v_s30=350.0,
            region="South America",
            event_type=event_type,
        )
    )
    actual = vectorized_results(m)
    for i, j in itertools.product(range(3), range(4)):
        expected = scalar_results(
            "South America", event_type[i, 0], mag[i, 0], dist_rup[j], 350.0
        )
        for key, value in actual.items():
            if key[0] == "duration_for_energy":
                for k in range(3):
                    assert value[k].shape == (3, 4)
                    assert_close(value[k][i, j], expected[key][k], err_msg=str(key))
            else:
                assert value.shape == (3, 4)
                assert_close(value[i, j], expected[key], err_msg=str(key))


def scenario_with(region, event_type):
    n = max(np.size(region), np.size(event_type))
    return Scenario(
        mag=np.full(n, 7.0),
        dist_rup=100.0,
        v_s30=400.0,
        region=region,
        event_type=event_type,
    )


def test_invalid_interface_taiwan_entries_raise():
    with pytest.raises(
        ValueError, match=r"interface earthquakes in Taiwan. \(1 of 3 scenarios\)"
    ):
        PinillaRamosEtAl2024(
            scenario_with(
                np.array(["Japan", "Taiwan", "Taiwan"]),
                np.array(["interface", "interface", "slab"]),
            )
        )


def test_invalid_event_type_entries_raise():
    # Like a scalar scenario, the invalid value is replaced with the default
    # (None) with a warning, and then an error is raised
    with pytest.warns(UserWarning, match="event_type has 1 of 3 values"):
        with pytest.raises(
            ValueError,
            match=r"event_type must be 'interface' or 'slab', got: None "
            r"\(1 of 3 scenarios\)",
        ):
            PinillaRamosEtAl2024(
                scenario_with("Japan", np.array(["slab", "crustal", "interface"]))
            )


def test_invalid_region_entries_raise():
    with pytest.warns(UserWarning, match="region has 2 of 3 values"):
        with pytest.raises(
            ValueError, match=r"region must be one of .* got: None \(2 of 3 scenarios\)"
        ):
            PinillaRamosEtAl2024(
                scenario_with(np.array(["Japan", "California", "Mars"]), "slab")
            )


def test_invalid_entries_are_reported_for_broadcast_scenarios():
    # A scalar invalid event type applies to all of the scenarios
    with pytest.warns(UserWarning):
        with pytest.raises(ValueError, match=r"got: None \(3 of 3 scenarios\)"):
            PinillaRamosEtAl2024(scenario_with(np.array(["Japan"] * 3), "crustal"))


def test_scalar_error_messages_are_unchanged():
    with pytest.warns(UserWarning):
        with pytest.raises(
            ValueError, match=r"^event_type must be 'interface' or 'slab', got: None$"
        ):
            PinillaRamosEtAl2024(scenario_with("Japan", "crustal").copy_with(mag=7.0))
    with pytest.raises(
        ValueError, match=r"^No model available for interface earthquakes in Taiwan.$"
    ):
        PinillaRamosEtAl2024(scenario_with("Taiwan", "interface").copy_with(mag=7.0))


def test_mismatched_lengths_raise():
    with pytest.raises(ValueError):
        PinillaRamosEtAl2024(
            Scenario(
                mag=np.array([6.0, 7.0]),
                dist_rup=np.array([10.0, 20.0, 30.0]),
                v_s30=400.0,
                region="Japan",
                event_type="slab",
            )
        )


@pytest.mark.parametrize("energy", [0.5, 0.75, 0.95])
def test_duration_model_arrays(grid, energy):
    rows, scenario, expected = grid
    actual = duration_model(
        scenario.mag,
        scenario.dist_rup,
        scenario.v_s30,
        scenario.region,
        scenario.event_type,
        energy,
    )
    assert len(actual) == 3
    for i in range(3):
        assert_close(
            actual[i],
            np.array([e[("duration_for_energy", energy)][i] for e in expected]),
        )


@pytest.mark.parametrize("energy", [0.5, 0.75, 0.95])
def test_duration_model_scalars(energy):
    actual = duration_model(7.0, 100.0, 400.0, "Japan", "interface", energy)
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=7.0, dist_rup=100.0, v_s30=400.0, region="Japan", event_type="interface"
        )
    )
    expected = m.duration_for_energy(energy)
    for a, e in zip(actual, expected):
        assert a.shape == (1,)
        assert_close(a[0], e)


def test_duration_model_broadcasts():
    median, plus, minus = duration_model(
        np.array([6.0, 7.0, 8.0]), 100.0, 400.0, "Japan", "slab", 0.95
    )
    assert median.shape == plus.shape == minus.shape == (3,)
    assert np.all(np.diff(median) > 0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])


@pytest.mark.parametrize(
    "event_type,region", [("interface", "Japan"), ("slab", "Taiwan")]
)
@pytest.mark.parametrize("threshold", ["D5-10", "D5-50", "D5-75", "D5-95"])
def test_d5x_matches_duration_for_energy(event_type, region, threshold):
    # d5x_median and d5x_sigma use the conditional model, like duration_for_energy
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=np.array([5.5, 7.0, 8.0]),
            dist_rup=np.array([20.0, 100.0, 300.0]),
            v_s30=np.array([200.0, 400.0, 1000.0]),
            region=region,
            event_type=event_type,
        )
    )
    energy = int(threshold[3:]) / 100
    median, plus, _ = m.duration_for_energy(energy)
    n = 0.15 if event_type == "interface" else 0.25
    assert_close(m.d5x_median(threshold), median)
    np.testing.assert_allclose(
        (m.d5x_median(threshold) ** n + m.d5x_sigma(threshold)) ** (1 / n),
        plus,
        rtol=1e-12,
    )


@pytest.mark.parametrize("threshold", ["D5-99", "D5-12", "D5-", "5-95", 95])
def test_d5x_unsupported_thresholds_raise(threshold):
    m = PinillaRamosEtAl2024(
        Scenario(
            mag=7.0, dist_rup=100.0, v_s30=400.0, region="Japan", event_type="interface"
        )
    )
    with pytest.raises(ValueError, match="Unsupported energy threshold"):
        m.d5x_median(threshold)
