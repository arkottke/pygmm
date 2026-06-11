#!/usr/bin/env python3
"""Tests for Pinilla-Ramos et al. (2024) subduction duration model."""

import os

import numpy as np
import pandas as pd
import pytest

from pygmm import Scenario
from pygmm.pinilla_ramos_et_al_2024 import PinillaRamosEtAl2024


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
