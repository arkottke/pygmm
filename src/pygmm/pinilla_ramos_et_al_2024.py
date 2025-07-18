"""Pinilla-Ramos et al. (2024, :cite:`pinilla-ramos24`) subduction duration model."""

from typing import Optional, Tuple, Union

import numpy as np

from . import model
from .types import ArrayLike

__author__ = "Albert Kottke"


class PinillaRamosEtAl2024(model.Model):
    """Pinilla-Ramos et al. (2024, :cite:`pinilla-ramos24`) subduction duration model.

    This model predicts significant duration for different energy thresholds
    (D5-X) where X represents the percentage of cumulative energy. The model
    is specifically designed for subduction zone earthquakes and distinguishes
    between interface and intraslab (slab) events across different tectonic regions.

    The model supports energy thresholds from D5-10 to D5-95 and includes
    regional adjustments for Japan, New Zealand, South America, and Taiwan
    (Taiwan only for slab events).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario. Must include event_type ('interface' or 'slab')
        and region specifications.

    Notes
    -----
    The model requires the following scenario parameters:
    - event_type: 'interface' for interface events, 'slab' for intraslab events
    - region: 'Japan', 'New Zealand', 'South America', or 'Taiwan' (slab only)
    - mag: moment magnitude (Mw)
    - dist_rup: rupture distance (km)
    - v_s30: time-averaged shear-wave velocity in top 30m (m/s)

    References
    ----------
    .. [1] Pinilla-Ramos et al. (2024). "Subduction Zone Significant Duration
           Model for Interface and Intraslab Earthquakes."
    """

    NAME = "Pinilla-Ramos et al. (2024)"
    ABBREV = "PR24"

    PARAMS = [
        model.NumericParameter("mag", True, 4.5, 8.5),
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("v_s30", True, 150, 2000),
        model.CategoricalParameter("event_type", True, ["interface", "slab"]),
        model.CategoricalParameter(
            "region", True, ["Japan", "New Zealand", "South America", "Taiwan"]
        ),
    ]

    # Energy thresholds supported by the model
    ENERGY_THRESHOLDS = np.array(
        [
            0.10,
            0.15,
            0.20,
            0.25,
            0.30,
            0.35,
            0.40,
            0.45,
            0.50,
            0.55,
            0.60,
            0.65,
            0.70,
            0.75,
            0.80,
            0.85,
            0.90,
            0.95,
        ]
    )

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)

        # Validate region-event_type combination
        self._validate_region_event_type()

        # Calculate the duration for D5-75 by default
        self._duration, self._duration_plus_sigma, self._duration_minus_sigma = self._calc_duration(
            0.75
        )

    def _validate_region_event_type(self) -> None:
        """Validate that the region-event_type combination is supported."""
        event_type = self._scenario.event_type
        region = self._scenario.region

        if event_type == "interface" and region == "Taiwan":
            raise ValueError("No model available for interface earthquakes in Taiwan.")

        if event_type not in ["interface", "slab"]:
            raise ValueError(f"event_type must be 'interface' or 'slab', got: {event_type}")

        valid_regions = ["Japan", "New Zealand", "South America", "Taiwan"]
        if region not in valid_regions:
            raise ValueError(f"region must be one of {valid_regions}, got: {region}")

    def _calc_d575_median(self) -> float:
        """Calculate the median D5-75 duration.

        Returns
        -------
        float
            Median D5-75 duration in seconds
        """
        s = self._scenario
        mag = float(s.mag)
        rrup = float(s.dist_rup)
        vs30 = float(s.v_s30)
        region = s.region
        event_type = s.event_type

        # Initialize variables
        c3_1 = c3_12 = mag_th = c1 = 0.0
        c3_base = c4_1 = 0.0
        r1 = v3 = 0.0

        if event_type == "interface":
            # Interface earthquake coefficients
            c3_1, c3_12, mag_th, c1 = 0.031, 0.029, 7.0, 5.268

            if region == "Japan":
                c3_base, c4_1 = 0.03, -2.485
            elif region == "New Zealand":
                c3_base, c4_1 = 0.024, -2.420
            elif region == "South America":
                c3_base, c4_1 = 0.024, -2.248

            r1, v3 = 250, 3100
            c2 = 0.275

            # Source term
            source_component = c1 * 10 ** ((mag - mag_th) * c2)

        elif event_type == "slab":
            # Slab earthquake coefficients
            c3_1, c3_12, mag_th, c1 = 0.013, 0.020, 6.0, 0.340

            if region == "Japan":
                c3_base, c4_1 = 0.046, -1.480
            elif region == "New Zealand":
                c3_base, c4_1 = 0.032, -1.252
            elif region == "South America":
                c3_base, c4_1 = 0.042, -1.510
            elif region == "Taiwan":
                c3_base, c4_1 = 0.038, -1.151

            r1, v3 = 250, 3100
            d1, m0, c2_base = 2.770, 4.250, 0.5

            # Source term (magnitude-dependent)
            if mag < mag_th:
                source_component = c1 * 10 ** ((mag - m0) * c2_base)
            else:
                source_component = c1 * 10 ** ((mag_th - m0) * c2_base) + d1 * (mag - mag_th) / (
                    8.5 - mag_th
                )

        # Site term
        site_component = c4_1 * np.log(vs30 / v3)

        # Path term
        if rrup < r1:
            path_component = c3_1 * rrup
        else:
            path_component = c3_12 * (rrup - r1) + r1 * c3_1

        # Main path term
        main_path = rrup * c3_base

        return source_component + site_component + path_component + main_path

    @property
    def d575_median(self) -> float:
        """D5-75 median duration in seconds."""
        return self._calc_d575_median()

    @property
    def d575_sigma(self) -> float:
        """D5-75 standard deviation of logarithmic duration."""
        return self._calc_d575_sigma()

    def d5x_median(self, energy_threshold: str) -> float:
        """Calculate median duration for specified energy threshold.

        Parameters
        ----------
        energy_threshold : str
            Energy threshold (e.g., "D5-10", "D5-95")

        Returns
        -------
        float
            Median duration in seconds
        """
        # Convert D5-75 to other thresholds using scaling relationships
        d575_med = self.d575_median

        # Scaling factors from paper (approximate)
        scaling_factors = {
            "D5-10": 0.25,
            "D5-20": 0.35,
            "D5-30": 0.50,
            "D5-40": 0.65,
            "D5-50": 0.75,
            "D5-60": 0.85,
            "D5-70": 0.95,
            "D5-75": 1.00,  # Reference
            "D5-80": 1.10,
            "D5-85": 1.25,
            "D5-90": 1.45,
            "D5-95": 1.80,
        }

        if energy_threshold not in scaling_factors:
            raise ValueError(f"Unsupported energy threshold: {energy_threshold}")

        return d575_med * scaling_factors[energy_threshold]

    def d5x_sigma(self, energy_threshold: str) -> float:
        """Calculate sigma for specified energy threshold.

        Parameters
        ----------
        energy_threshold : str
            Energy threshold (e.g., "D5-10", "D5-95")

        Returns
        -------
        float
            Standard deviation of logarithmic duration
        """
        # For simplicity, use same sigma as D5-75
        # Could implement threshold-specific sigmas if data available
        return self.d575_sigma

    def _calc_d575_sigma(self) -> float:
        """Calculate the standard deviation for D5-75 duration.

        Returns
        -------
        float
            Standard deviation of D5-75 duration
        """
        s = self._scenario
        mag = float(s.mag)
        rrup = float(s.dist_rup)
        vs30 = float(s.v_s30)
        event_type = s.event_type

        if event_type == "interface":
            a0, a1, a2, b1, b2, c1 = 0.3107, -0.0393, 0.0062, -0.0519, 0.0041, 0.0028
        elif event_type == "slab":
            a0, a1, a2, b1, b2, c1 = 0.3035, -0.0188, 0.004, -0.0288, 0.0019, 0.0038

        return (
            a0
            + a1 * rrup / 100
            + a2 * (rrup / 100) ** 2
            + b1 * mag
            + b2 * mag**2
            + np.log(vs30) * c1
        )

    def _get_energy_coefficients(self, energy: float) -> Tuple[float, ...]:
        """Get the model coefficients for a specific energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold (e.g., 0.75 for D5-75)

        Returns
        -------
        tuple
            Model coefficients for the specified energy threshold
        """
        # Find the closest supported energy threshold
        idx_energy = np.argmin(np.abs(self.ENERGY_THRESHOLDS - energy))

        event_type = self._scenario.event_type

        if event_type == "interface":
            c_median = np.array(
                [
                    0.04020966,
                    0.08388851,
                    0.13077859,
                    0.1815641,
                    0.23651304,
                    0.29491758,
                    0.35754282,
                    0.42528554,
                    0.4974965,
                    0.57657072,
                    0.66322614,
                    0.76042762,
                    0.87127838,
                    1.0,
                    1.156691,
                    1.35981923,
                    1.64806185,
                    2.14387062,
                ]
            )[idx_energy]

            a0 = np.array(
                [
                    0.00733073,
                    0.0082492,
                    0.00374054,
                    -0.00437519,
                    -0.01719855,
                    -0.03221052,
                    -0.04872026,
                    -0.06537968,
                    -0.0775099,
                    -0.08107197,
                    -0.07769041,
                    -0.06689982,
                    -0.04362409,
                    0.0,
                    0.07075031,
                    0.19400743,
                    0.40879893,
                    0.85474174,
                ]
            )[idx_energy]

            m1 = np.array(
                [
                    0.00014309,
                    0.00117448,
                    0.00259631,
                    0.00410708,
                    0.00587722,
                    0.00790171,
                    0.01010809,
                    0.01227336,
                    0.01365048,
                    0.01360863,
                    0.01233418,
                    0.00990628,
                    0.00596924,
                    0.0,
                    -0.00864832,
                    -0.02218093,
                    -0.04374555,
                    -0.07625391,
                ]
            )[idx_energy]

            r1 = np.array(
                [
                    0.00240356,
                    0.00495213,
                    0.00759418,
                    0.0103042,
                    0.01287056,
                    0.01498933,
                    0.01646171,
                    0.01716948,
                    0.01753574,
                    0.01709783,
                    0.01588001,
                    0.01298844,
                    0.00821674,
                    0.0,
                    -0.01349763,
                    -0.03821902,
                    -0.08847199,
                    -0.21293647,
                ]
            )[idx_energy]

            v1 = np.array(
                [
                    0.00547392,
                    0.01140383,
                    0.01633181,
                    0.02028365,
                    0.023049,
                    0.02498088,
                    0.02612515,
                    0.02645767,
                    0.02556579,
                    0.02366726,
                    0.02005288,
                    0.01446657,
                    0.00779351,
                    0.0,
                    -0.00992416,
                    -0.02254522,
                    -0.04507584,
                    -0.07641776,
                ]
            )[idx_energy]

            rho_c_d575 = np.array(
                [
                    -0.04461951,
                    0.06461454,
                    0.14724378,
                    0.20758114,
                    0.2514994,
                    0.28656202,
                    0.31098974,
                    0.33129013,
                    0.34623991,
                    0.3543134,
                    0.3577682,
                    0.3535797,
                    0.33183304,
                    0.0,
                    -0.38333019,
                    -0.44116537,
                    -0.4875535,
                    -0.53284828,
                ]
            )[idx_energy]

            sigma_c = np.array(
                [
                    0.01395505,
                    0.02371422,
                    0.03224683,
                    0.03982668,
                    0.04629945,
                    0.05151767,
                    0.05540747,
                    0.05772215,
                    0.05816558,
                    0.05630814,
                    0.05154701,
                    0.04259618,
                    0.02792395,
                    0.0,
                    0.04378463,
                    0.1063117,
                    0.2196963,
                    0.48117887,
                ]
            )[idx_energy]

            n2 = 0.15

        elif event_type == "slab":
            c_median = np.array(
                [
                    0.03784117,
                    0.08352235,
                    0.13083025,
                    0.18030833,
                    0.23239066,
                    0.28891588,
                    0.34921976,
                    0.41510462,
                    0.48731871,
                    0.56681198,
                    0.65450169,
                    0.75358998,
                    0.86718004,
                    1.0,
                    1.16409534,
                    1.38075492,
                    1.69269723,
                    2.22279683,
                ]
            )[idx_energy]

            a0 = np.array(
                [
                    -0.02602886,
                    -0.01820382,
                    -0.0179457,
                    -0.02417655,
                    -0.03142039,
                    -0.03858319,
                    -0.04451747,
                    -0.04873467,
                    -0.04927577,
                    -0.04910907,
                    -0.04485417,
                    -0.03689654,
                    -0.02248334,
                    0.0,
                    0.04127725,
                    0.10055556,
                    0.20353897,
                    0.38563394,
                ]
            )[idx_energy]

            m1 = np.array(
                [
                    0.0052392,
                    0.00598378,
                    0.00764049,
                    0.00998344,
                    0.01243448,
                    0.01438996,
                    0.01589705,
                    0.01668025,
                    0.01651024,
                    0.01585809,
                    0.01396992,
                    0.01083706,
                    0.00631921,
                    0.0,
                    -0.00947123,
                    -0.02266151,
                    -0.04240792,
                    -0.07401577,
                ]
            )[idx_energy]

            r1 = np.array(
                [
                    0.00226879,
                    0.00419625,
                    0.00701796,
                    0.01021163,
                    0.0129733,
                    0.01554978,
                    0.01720173,
                    0.0182664,
                    0.01846989,
                    0.01753247,
                    0.01580923,
                    0.01256225,
                    0.00751923,
                    0.0,
                    -0.01202052,
                    -0.03231773,
                    -0.07103404,
                    -0.16362926,
                ]
            )[idx_energy]

            v1 = np.array(
                [
                    0.00409751,
                    0.01193314,
                    0.01968596,
                    0.02664136,
                    0.03281428,
                    0.03781827,
                    0.04107292,
                    0.04268421,
                    0.04266614,
                    0.04062058,
                    0.03576105,
                    0.0275821,
                    0.0162236,
                    0.0,
                    -0.02139202,
                    -0.05326305,
                    -0.10313415,
                    -0.21286673,
                ]
            )[idx_energy]

            rho_c_d575 = np.array(
                [
                    0.23703883,
                    0.2062694,
                    0.26652301,
                    0.33081508,
                    0.38018048,
                    0.41419228,
                    0.43660475,
                    0.45004969,
                    0.45471735,
                    0.45897072,
                    0.45725023,
                    0.44906822,
                    0.42848698,
                    0.0,
                    -0.46673359,
                    -0.51548056,
                    -0.55293589,
                    -0.60041303,
                ]
            )[idx_energy]

            sigma_c = np.array(
                [
                    0.01952566,
                    0.02989091,
                    0.03838868,
                    0.04680806,
                    0.05476812,
                    0.06195816,
                    0.06758683,
                    0.0712496,
                    0.07242325,
                    0.07068507,
                    0.0647861,
                    0.05349118,
                    0.03450956,
                    0.0,
                    0.05308633,
                    0.1307839,
                    0.26720298,
                    0.56903831,
                ]
            )[idx_energy]

            n2 = 0.25

        return c_median, a0, m1, r1, v1, rho_c_d575, sigma_c, n2

    def _calc_duration(self, energy: float = 0.75) -> Tuple[float, float, float]:
        """Calculate duration for specified energy threshold.

        Parameters
        ----------
        energy : float, optional
            Energy threshold (0.75 for D5-75, 0.95 for D5-95, etc.)

        Returns
        -------
        tuple
            (median_duration, duration_plus_sigma, duration_minus_sigma) in seconds
        """
        s = self._scenario
        mag = float(s.mag)
        rrup = float(s.dist_rup)
        vs30 = float(s.v_s30)
        event_type = s.event_type

        if np.isclose(energy, 0.75):
            # Use base D5-75 model
            median = self._calc_d575_median()
            sigma = self._calc_d575_sigma()

            # Transform from log space
            n = 0.15 if event_type == "interface" else 0.25
            duration_plus = (median**n + sigma) ** (1 / n)
            duration_minus = (median**n - sigma) ** (1 / n)

            return median, duration_plus, duration_minus

        else:
            # Use conditional model for other energy thresholds
            d575_median = self._calc_d575_median()
            sigma_575 = self._calc_d575_sigma()

            # Get energy-specific coefficients
            c_median, a0, m1, r1, v1, rho_c_d575, sigma_c, n2 = self._get_energy_coefficients(
                energy
            )

            # Calculate conditional model components
            c_ratio = c_median + a0 + m1 * mag + r1 * rrup / 100 + v1 * np.log(vs30 / 3100)

            # Median duration for this energy threshold
            d5x_median = d575_median * c_ratio

            # Standard deviation calculation
            var_5x = (
                sigma_575**2 * c_ratio ** (2 * n2)
                + n2**2 * sigma_c**2 * d575_median ** (2 * n2) * c_ratio ** (2 * n2 - 2)
                + 2
                * n2
                * rho_c_d575
                * c_ratio ** (2 * n2 - 1)
                * d575_median**n2
                * sigma_575
                * sigma_c
            )
            d5x_sigma = np.sqrt(var_5x)

            # Transform from log space
            n = 0.15 if event_type == "interface" else 0.25
            duration_plus = (d5x_median**n + d5x_sigma) ** (1 / n)
            duration_minus = (d5x_median**n - d5x_sigma) ** (1 / n)

            return d5x_median, duration_plus, duration_minus

    @property
    def duration(self) -> float:
        """D5-75 duration in seconds."""
        return self._duration

    @property
    def duration_plus_sigma(self) -> float:
        """D5-75 duration plus one standard deviation in seconds."""
        return self._duration_plus_sigma

    @property
    def duration_minus_sigma(self) -> float:
        """D5-75 duration minus one standard deviation in seconds."""
        return self._duration_minus_sigma

    def duration_for_energy(self, energy: float) -> Tuple[float, float, float]:
        """Calculate duration for specified energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold. For D5-95 use 0.95, for D5-45 use 0.45, etc.
            Valid range is approximately 0.10 to 0.95.

        Returns
        -------
        tuple
            (median_duration, duration_plus_sigma, duration_minus_sigma) in seconds

        Raises
        ------
        ValueError
            If energy threshold is outside valid range

        Examples
        --------
        >>> scenario = Scenario(mag=7.0, dist_rup=100, v_s30=300,
        ...                     event_type='interface', region='Japan')
        >>> model = PinillaRamosEtAl2024(scenario)
        >>> d5_95_median, d5_95_plus, d5_95_minus = model.duration_for_energy(0.95)
        """
        if energy < 0.10 or energy > 0.95:
            raise ValueError(f"Energy threshold {energy} is outside valid range [0.10, 0.95]")

        if not np.any(np.isclose(self.ENERGY_THRESHOLDS, energy, atol=1e-3)):
            # Find closest supported threshold
            closest_idx = np.argmin(np.abs(self.ENERGY_THRESHOLDS - energy))
            closest_energy = self.ENERGY_THRESHOLDS[closest_idx]
            import warnings

            warnings.warn(
                f"Energy threshold {energy} not directly supported. Using closest "
                f"available threshold {closest_energy}."
            )
            energy = closest_energy

        return self._calc_duration(energy)


# Legacy function for backward compatibility
def duration_model(
    mag: ArrayLike, rrup: ArrayLike, vs30: ArrayLike, region: str, eq_type: str, energy: float
) -> Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Legacy duration model function for backward compatibility.

    Parameters
    ----------
    mag : array_like
        Moment magnitude
    rrup : array_like
        Rupture distance in km
    vs30 : array_like
        Time-averaged shear-wave velocity in the top 30 m in m/s
    region : str
        Tectonic region ('Japan', 'New Zealand', 'South America', 'Taiwan')
    eq_type : str
        Event type ('interface' or 'slab')
    energy : float
        Energy threshold (e.g., 0.75 for D5-75)

    Returns
    -------
    tuple
        (median_duration, duration_plus_sigma, duration_minus_sigma)
    """
    from .model import Scenario

    # Convert inputs to arrays
    mag_arr = np.atleast_1d(mag)
    rrup_arr = np.atleast_1d(rrup)
    vs30_arr = np.atleast_1d(vs30)

    # Create scenario and model
    scenario = Scenario(
        mag=mag_arr[0], dist_rup=rrup_arr[0], v_s30=vs30_arr[0], event_type=eq_type, region=region
    )
    model_inst = PinillaRamosEtAl2024(scenario)

    if np.isclose(energy, 0.75):
        duration, plus_sigma, minus_sigma = (
            model_inst.duration,
            model_inst.duration_plus_sigma,
            model_inst.duration_minus_sigma,
        )
    else:
        duration, plus_sigma, minus_sigma = model_inst.duration_for_energy(energy)

    # Return as arrays to match expected return type
    return np.array([duration]), np.array([plus_sigma]), np.array([minus_sigma])


if __name__ == "__main__":
    """
    How to use:
    The model predicts significant duration for subduction earthquakes.
    Supports both interface and intraslab (slab) events.
    To calculate D5-X, define energy as X as a decimal. For example, for D5-45, use 0.45.
    """
    from .model import Scenario

    # Example usage - Interface earthquake in Japan
    scenario = Scenario(
        mag=7.0, dist_rup=100.0, v_s30=300.0, event_type="interface", region="Japan"
    )
    model_inst = PinillaRamosEtAl2024(scenario)

    print(f"D5-75 duration: {model_inst.duration:.2f} seconds")
    print(f"D5-95 duration: {model_inst.duration_for_energy(0.95)[0]:.2f} seconds")

    # Example usage - Slab earthquake in South America
    scenario_slab = Scenario(
        mag=6.5, dist_rup=80.0, v_s30=400.0, event_type="slab", region="South America"
    )
    model_slab = PinillaRamosEtAl2024(scenario_slab)

    print(f"Slab D5-75 duration: {model_slab.duration:.2f} seconds")
