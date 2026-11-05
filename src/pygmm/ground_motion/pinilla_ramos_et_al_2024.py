"""Pinilla-Ramos et al. (2024, :cite:p:`pinilla-ramos24`) subduction duration model."""

import numpy as np

from .. import model
from .._types import ArrayLike
from ..contracts import Duration
from ..registry import register

__author__ = "Albert Kottke"


@register(provides=("duration",), tectonic="subduction")
class PinillaRamosEtAl2024(model.Model):
    """Pinilla-Ramos et al. (2024, :cite:p:`pinilla-ramos24`) subduction duration model.

    This model predicts significant duration for different energy thresholds
    (D5-X) where X represents the percentage of cumulative energy. The model
    is specifically designed for subduction zone earthquakes and distinguishes
    between interface and intraslab (slab) events across different tectonic regions.

    The model supports energy thresholds from D5-10 to D5-95 and includes
    regional adjustments for Japan, New Zealand, South America, and Taiwan
    (Taiwan only for slab events).

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``v_s30``, ``event_type``, and ``region``) can be a scalar or an array, and
    the arrays are broadcast against each other. For a scalar scenario, the
    durations and standard deviations are scalars, as before. For arrays of
    scenarios, ``d575_median``, ``d575_sigma``, ``d5x_median``, ``d5x_sigma``,
    ``duration``, ``duration_plus_sigma``, ``duration_minus_sigma``, and each of
    the three values returned by ``duration_for_energy`` are arrays with the
    broadcast shape of the scenario values, e.g., shape (N,) for N scenarios.
    The energy threshold applies to all scenarios.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario. Must include event_type ('interface' or 'slab')
        and region specifications.

    Examples
    --------
    >>> import numpy as np
    >>> from pygmm import Scenario
    >>> from pygmm.ground_motion.pinilla_ramos_et_al_2024 import PinillaRamosEtAl2024
    >>> s = Scenario(
    ...     mag=np.array([6.0, 7.0, 8.0]), dist_rup=np.array([50.0, 100.0, 200.0]),
    ...     v_s30=400.0, event_type=np.array(["slab", "interface", "interface"]),
    ...     region="Japan")
    >>> m = PinillaRamosEtAl2024(s)
    >>> m.duration.shape
    (3,)
    >>> median, plus, minus = m.duration_for_energy(0.95)
    >>> median.shape
    (3,)

    Notes
    -----
    The model requires the following scenario parameters:
    - event_type: 'interface' for interface events, 'slab' for intraslab events
    - region: 'Japan', 'New Zealand', 'South America', or 'Taiwan' (slab only)
    - mag: moment magnitude (Mw)
    - dist_rup: rupture distance [km]
    - v_s30: time-averaged shear-wave velocity in top 30m [m/s]

    References
    ----------
    :cite:t:`pinilla-ramos24`

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

        s = self._scenario
        # Broadcast shape of the scenario values, () for a scalar scenario
        self._shape = np.broadcast_shapes(
            *(
                np.shape(s[name])
                for name in ("mag", "dist_rup", "v_s30", "event_type", "region")
            )
        )

        # Validate region-event_type combination
        self._validate_region_event_type()

        # After the validation, each event is either interface or slab
        self._is_interface = model.equals(s.event_type, "interface")
        self._d575_median = self._calc_d575_median()
        self._d575_sigma = self._calc_d575_sigma()

        # Calculate the duration for D5-75 by default
        self._duration, self._duration_plus_sigma, self._duration_minus_sigma = (
            self._calc_duration(0.75)
        )

    def _output(self, value) -> ArrayLike:
        """Broadcast a computed value to the scenario shape.

        Scalar scenarios give a scalar (:class:`numpy.float64`), as before.
        """
        value = np.asarray(value)
        if value.shape != self._shape:
            value = np.broadcast_to(value, self._shape).copy()
        return value[()] if value.ndim == 0 else value

    def _validate_region_event_type(self) -> None:
        """Validate that the region-event_type combination is supported.

        For arrays of scenarios, an error is raised if any scenario is not
        supported.
        """
        event_type = self._scenario.event_type
        region = self._scenario.region
        shape = self._shape
        vectorized = len(shape) > 0

        def describe(values, invalid) -> str:
            """Describe the invalid values for the error message."""
            if not vectorized:
                return str(values)
            invalid = np.broadcast_to(invalid, shape)
            values = np.broadcast_to(np.asarray(values, dtype=object), shape)
            names = ", ".join(sorted({str(v) for v in values[invalid]}))
            return f"{names} ({np.count_nonzero(invalid)} of {invalid.size} scenarios)"

        is_interface = model.equals(event_type, "interface")
        is_slab = model.equals(event_type, "slab")
        is_taiwan = model.equals(region, "Taiwan")

        interface_taiwan = is_interface & is_taiwan
        if np.any(interface_taiwan):
            msg = "No model available for interface earthquakes in Taiwan."
            if vectorized:
                count = np.count_nonzero(np.broadcast_to(interface_taiwan, shape))
                msg += f" ({count} of {int(np.prod(shape))} scenarios)"
            raise ValueError(msg)

        invalid = ~(is_interface | is_slab)
        if np.any(invalid):
            raise ValueError(
                "event_type must be 'interface' or 'slab', got: "
                + describe(event_type, invalid)
            )

        valid_regions = ["Japan", "New Zealand", "South America", "Taiwan"]
        invalid = ~np.logical_or.reduce(
            [model.equals(region, r) for r in valid_regions]
        )
        if np.any(invalid):
            raise ValueError(
                f"region must be one of {valid_regions}, got: "
                + describe(region, invalid)
            )

    def _calc_d575_median(self) -> np.ndarray:
        """Calculate the median D5-75 duration.

        Returns
        -------
        array_like
            Median D5-75 duration in seconds (unbroadcast)
        """
        s = self._scenario
        mag = np.asarray(s.mag, dtype=float)
        rrup = np.asarray(s.dist_rup, dtype=float)
        vs30 = np.asarray(s.v_s30, dtype=float)
        is_interface = self._is_interface

        # Distance coefficients
        c3_1 = np.where(is_interface, 0.031, 0.013)
        c3_12 = np.where(is_interface, 0.029, 0.020)
        r1, v3 = 250, 3100

        # Regional coefficients for each (event_type, region) combination
        regional = [
            ("interface", "Japan", 0.03, -2.485),
            ("interface", "New Zealand", 0.024, -2.420),
            ("interface", "South America", 0.024, -2.248),
            ("slab", "Japan", 0.046, -1.480),
            ("slab", "New Zealand", 0.032, -1.252),
            ("slab", "South America", 0.042, -1.510),
            ("slab", "Taiwan", 0.038, -1.151),
        ]
        is_region = {
            region: model.equals(s.region, region)
            for region in {r[1] for r in regional}
        }
        conds = [
            (is_interface if event_type == "interface" else ~is_interface)
            & is_region[region]
            for event_type, region, _, _ in regional
        ]
        c3_base = np.select(conds, [r[2] for r in regional], default=0.0)
        c4_1 = np.select(conds, [r[3] for r in regional], default=0.0)

        # Source term. Interface: c1 = 5.268, mag_th = 7.0, c2 = 0.275. Slab:
        # c1 = 0.340, mag_th = 6.0, d1 = 2.770, m0 = 4.250, c2_base = 0.5, with
        # magnitude scaling that changes at mag_th.
        source_interface = 5.268 * 10 ** ((mag - 7.0) * 0.275)
        source_slab_small = 0.340 * 10 ** ((mag - 4.250) * 0.5)
        source_slab_large = 0.340 * 10 ** ((6.0 - 4.250) * 0.5) + 2.770 * (
            mag - 6.0
        ) / (8.5 - 6.0)
        source_component = np.select(
            [is_interface, mag < 6.0],
            [source_interface, source_slab_small],
            default=source_slab_large,
        )

        # Site term
        site_component = c4_1 * np.log(vs30 / v3)

        # Path term
        path_component = np.where(
            rrup < r1, c3_1 * rrup, c3_12 * (rrup - r1) + r1 * c3_1
        )

        # Main path term
        main_path = rrup * c3_base

        return source_component + site_component + path_component + main_path

    @property
    def d575_median(self) -> ArrayLike:
        """D5-75 median duration in seconds."""
        return self._output(self._d575_median)

    @property
    def d575_sigma(self) -> ArrayLike:
        """D5-75 standard deviation of logarithmic duration."""
        return self._output(self._d575_sigma)

    def d5x_median(self, energy_threshold: str) -> ArrayLike:
        """Calculate median duration for specified energy threshold.

        The median is computed with the conditional model, as for
        :meth:`duration_for_energy`.

        Parameters
        ----------
        energy_threshold : str
            Energy threshold (e.g., "D5-10", "D5-95"), from "D5-10" to "D5-95"
            in steps of 5%.

        Returns
        -------
        float or array_like
            Median duration in seconds, with the scenario shape
        """
        median, _ = self._calc_d5x(self._parse_energy_threshold(energy_threshold))
        return self._output(median)

    def d5x_sigma(self, energy_threshold: str) -> ArrayLike:
        """Calculate sigma for specified energy threshold.

        The standard deviation is computed with the conditional model, as for
        :meth:`duration_for_energy`, in the same units as :attr:`d575_sigma`.

        Parameters
        ----------
        energy_threshold : str
            Energy threshold (e.g., "D5-10", "D5-95"), from "D5-10" to "D5-95"
            in steps of 5%.

        Returns
        -------
        float or array_like
            Standard deviation of the duration, with the scenario shape
        """
        _, sigma = self._calc_d5x(self._parse_energy_threshold(energy_threshold))
        return self._output(sigma)

    def _parse_energy_threshold(self, energy_threshold: str) -> float:
        """Energy of a threshold name, e.g., 0.95 for "D5-95"."""
        try:
            if not energy_threshold.startswith("D5-"):
                raise ValueError
            energy = int(energy_threshold[len("D5-") :]) / 100
        except (AttributeError, ValueError):
            raise ValueError(
                f"Unsupported energy threshold: {energy_threshold}"
            ) from None
        if not np.any(np.isclose(self.ENERGY_THRESHOLDS, energy)):
            raise ValueError(f"Unsupported energy threshold: {energy_threshold}")
        return energy

    def _calc_d575_sigma(self) -> np.ndarray:
        """Calculate the standard deviation for D5-75 duration.

        Returns
        -------
        array_like
            Standard deviation of D5-75 duration (unbroadcast)
        """
        s = self._scenario
        mag = np.asarray(s.mag, dtype=float)
        rrup = np.asarray(s.dist_rup, dtype=float)
        vs30 = np.asarray(s.v_s30, dtype=float)

        a0, a1, a2, b1, b2, c1 = (
            np.where(self._is_interface, interface, slab)
            for interface, slab in zip(
                (0.3107, -0.0393, 0.0062, -0.0519, 0.0041, 0.0028),
                (0.3035, -0.0188, 0.004, -0.0288, 0.0019, 0.0038),
            )
        )

        return (
            a0
            + a1 * rrup / 100
            + a2 * (rrup / 100) ** 2
            + b1 * mag
            + b2 * mag**2
            + np.log(vs30) * c1
        )

    def _get_energy_coefficients(self, energy: float) -> tuple[np.ndarray, ...]:
        """Get the model coefficients for a specific energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold (e.g., 0.75 for D5-75)

        Returns
        -------
        tuple
            Model coefficients for the specified energy threshold, selected by
            the event type of each scenario
        """
        # Find the closest supported energy threshold
        idx_energy = np.argmin(np.abs(self.ENERGY_THRESHOLDS - energy))

        interface = self._event_type_coefficients("interface", idx_energy)
        slab = self._event_type_coefficients("slab", idx_energy)
        return tuple(
            np.where(self._is_interface, i, s) for i, s in zip(interface, slab)
        )

    @staticmethod
    def _event_type_coefficients(
        event_type: str, idx_energy: int | np.integer
    ) -> tuple[float, ...]:
        """Model coefficients for an event type and energy threshold index."""
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

    def _calc_d5x(self, energy: float) -> tuple[np.ndarray, np.ndarray]:
        """Calculate the median and standard deviation for an energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold (0.75 for D5-75, 0.95 for D5-95, etc.)

        Returns
        -------
        tuple
            (median_duration, sigma), unbroadcast. D5-75 uses the base model,
            and other thresholds use the conditional model.
        """
        if np.isclose(energy, 0.75):
            # Use base D5-75 model
            return self._d575_median, self._d575_sigma

        s = self._scenario
        mag = np.asarray(s.mag, dtype=float)
        rrup = np.asarray(s.dist_rup, dtype=float)
        vs30 = np.asarray(s.v_s30, dtype=float)
        is_interface = self._is_interface

        # Use conditional model for other energy thresholds
        d575_median = self._d575_median
        sigma_575 = self._d575_sigma

        # Get energy-specific coefficients
        c_median, a0, m1, r1, v1, rho_c_d575, sigma_c, n2 = (
            self._get_energy_coefficients(energy)
        )
        # Square of n2 from the scalar values (as for a scalar scenario)
        n2_sq = np.where(is_interface, 0.15**2, 0.25**2)

        # Calculate conditional model components
        c_ratio = c_median + a0 + m1 * mag + r1 * rrup / 100 + v1 * np.log(vs30 / 3100)

        # Median duration for this energy threshold
        d5x_median = d575_median * c_ratio

        # Standard deviation calculation
        var_5x = (
            sigma_575**2 * c_ratio ** (2 * n2)
            + n2_sq * sigma_c**2 * d575_median ** (2 * n2) * c_ratio ** (2 * n2 - 2)
            + 2
            * n2
            * rho_c_d575
            * c_ratio ** (2 * n2 - 1)
            * d575_median**n2
            * sigma_575
            * sigma_c
        )
        d5x_sigma = np.sqrt(var_5x)

        return d5x_median, d5x_sigma

    def _calc_duration(self, energy: float = 0.75) -> tuple[ArrayLike, ...]:
        """Calculate duration for specified energy threshold.

        Parameters
        ----------
        energy : float, optional
            Energy threshold (0.75 for D5-75, 0.95 for D5-95, etc.)

        Returns
        -------
        tuple
            (median_duration, duration_plus_sigma, duration_minus_sigma) in
            seconds. Each is a scalar for a scalar scenario, or an array with the
            scenario shape.
        """
        median, sigma = self._calc_d5x(energy)

        # Exponent of the transformation
        n = np.where(self._is_interface, 0.15, 0.25)

        # Transform from log space
        duration_plus = (median**n + sigma) ** (1 / n)
        duration_minus = (median**n - sigma) ** (1 / n)

        return (
            self._output(median),
            self._output(duration_plus),
            self._output(duration_minus),
        )

    @property
    def duration(self) -> ArrayLike:
        """D5-75 duration in seconds."""
        return self._duration

    @property
    def duration_plus_sigma(self) -> ArrayLike:
        """D5-75 duration plus one standard deviation in seconds."""
        return self._duration_plus_sigma

    @property
    def duration_minus_sigma(self) -> ArrayLike:
        """D5-75 duration minus one standard deviation in seconds."""
        return self._duration_minus_sigma

    def duration_for_energy(self, energy: float) -> tuple[ArrayLike, ...]:
        """Calculate duration for specified energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold. For D5-95 use 0.95, for D5-45 use 0.45, etc.
            Valid range is approximately 0.10 to 0.95. The same threshold is
            used for all scenarios.

        Returns
        -------
        tuple
            (median_duration, duration_plus_sigma, duration_minus_sigma) in
            seconds. Each is a scalar for a scalar scenario, or an array with the
            scenario shape.

        Raises
        ------
        ValueError
            If energy threshold is outside valid range

        Examples
        --------
        >>> from pygmm import Scenario
        >>> scenario = Scenario(mag=7.0, dist_rup=100, v_s30=300,
        ...                     event_type='interface', region='Japan')
        >>> model = PinillaRamosEtAl2024(scenario)
        >>> d5_95_median, d5_95_plus, d5_95_minus = model.duration_for_energy(0.95)
        """
        if energy < 0.10 or energy > 0.95:
            raise ValueError(
                f"Energy threshold {energy} is outside valid range [0.10, 0.95]"
            )

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

    def duration_model(self) -> Duration:
        """Return the :class:`~pygmm.contracts.Duration` contract.

        Sigma is applied in a transformed power space rather than
        lognormally, so the percentiles are reported directly and ``ln_std``
        is left unset.

        The values are floats for a scalar scenario and arrays for a vectorized
        scenario.
        """
        dur = model.as_scalar_or_array(self.duration)
        return Duration(
            duration=dur,
            d_5_75=dur,
            plus_sigma=model.as_scalar_or_array(self.duration_plus_sigma),
            minus_sigma=model.as_scalar_or_array(self.duration_minus_sigma),
        )


# Legacy function for backward compatibility
def duration_model(
    mag: ArrayLike,
    rrup: ArrayLike,
    vs30: ArrayLike,
    region: str,
    eq_type: str,
    energy: float,
) -> tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Legacy duration model function for backward compatibility.

    Parameters
    ----------
    mag : array_like
        Moment magnitude
    rrup : array_like
        Rupture distance in km
    vs30 : array_like
        Time-averaged shear-wave velocity in the top 30 m in m/s
    region : str or array_like
        Tectonic region ('Japan', 'New Zealand', 'South America', 'Taiwan')
    eq_type : str or array_like
        Event type ('interface' or 'slab')
    energy : float
        Energy threshold (e.g., 0.75 for D5-75)

    Returns
    -------
    tuple
        (median_duration, duration_plus_sigma, duration_minus_sigma). Each is an
        array with the broadcast shape of the inputs, or shape (1,) if all of the
        inputs are scalars.
    """
    from ..model import Scenario

    # Create scenario and model. Scalars are kept as scalars so that the
    # results are the same as for a scalar scenario.
    scenario = Scenario(
        mag=mag if np.ndim(mag) == 0 else np.asarray(mag),
        dist_rup=rrup if np.ndim(rrup) == 0 else np.asarray(rrup),
        v_s30=vs30 if np.ndim(vs30) == 0 else np.asarray(vs30),
        event_type=eq_type if np.ndim(eq_type) == 0 else np.asarray(eq_type),
        region=region if np.ndim(region) == 0 else np.asarray(region),
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
    return (
        np.atleast_1d(duration),
        np.atleast_1d(plus_sigma),
        np.atleast_1d(minus_sigma),
    )


if __name__ == "__main__":
    """
    How to use:
    The model predicts significant duration for subduction earthquakes.
    Supports both interface and intraslab (slab) events.
    To calculate D5-X, define energy as X as a decimal.

    For example, for D5-45, use 0.45.
    """
    from ..model import Scenario

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
