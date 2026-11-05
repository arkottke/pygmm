"""Pinilla-Ramos et al. (2023, :cite:p:`pinilla-ramos23`) duration model."""

import numpy as np

from .. import model
from ..contracts import Duration
from ..registry import register

__author__ = "Albert Kottke"


def _conditional_model() -> np.ndarray:
    """Conditional model coefficients for the energy thresholds."""
    header = [
        ("energy", "f8"),
        ("cmed", "f8"),
        ("a0", "f8"),
        ("r1x", "f8"),
        ("v1x", "f8"),
        ("std_c", "f8"),
        ("r", "f8"),
    ]

    conditional_model = np.zeros(19, dtype=header)

    conditional_model["energy"] = [0.05 + 0.05 * i for i in range(19)]
    # fmt: off
    conditional_model["cmed"] = [
        0, 0.157, 0.264, 0.342, 0.402, 0.455, 0.505, 0.553, 0.603, 0.654, 0.71,
        0.769, 0.835, 0.912, 1, 1.114, 1.273, 1.522, 2.014,
    ]
    conditional_model["a0"] = [
        0, -0.010798, -0.016831, -0.012831, 0.002943, 0.02267, 0.047579,
        0.076718, 0.107148, 0.136351, 0.115442, 0.092914, 0.067803, 0.034992, 0,
        -0.044725, -0.112447, -0.209689, -0.38092,
    ]
    conditional_model["r1x"] = [
        0, 0.0007, 0.0012, 0.0014, 0.0015, 0.0015, 0.0014, 0.0013, 0.0012,
        0.001, 0.0008, 0.0007, 0.0005, 0.0002, 0, -0.0003, -0.0006, -0.001,
        -0.0015,
    ]
    conditional_model["v1x"] = [
        0, 0.039, 0.0656, 0.0852, 0.1001, 0.1134, 0.1259, 0.1377, 0.1501,
        0.1587, 0.1365, 0.1105, 0.08, 0.0428, 0, -0.0512, -0.1197, -0.2111,
        -0.3589,
    ]
    conditional_model["std_c"] = [
        0, 0.156, 0.192, 0.205, 0.206, 0.202, 0.195, 0.187, 0.177, 0.163, 0.146,
        0.125, 0.097, 0.06, 0, 0.089, 0.21, 0.434, 0.907,
    ]
    conditional_model["r"] = [
        0, -0.083, 0.022, 0.078, 0.113, 0.137, 0.154, 0.167, 0.178, 0.188,
        0.198, 0.206, 0.209, 0.204, 0, -0.301, -0.361, -0.403, -0.452,
    ]
    # fmt: on
    return conditional_model


def _as_output(value: np.ndarray) -> float | np.ndarray:
    """Return a float for a scalar scenario and an array otherwise."""
    if np.ndim(value) == 0:
        return float(value)
    return np.asarray(value)


@register(provides=("duration",))
class PinillaRamosEtAl2023(model.Model):
    """Pinilla-Ramos et al. (2023, :cite:p:`pinilla-ramos23`) duration model.

    This model predicts significant duration for different energy thresholds
    (D5-X) where X represents the percentage of cumulative energy. The model
    can predict D5-75 (duration between 5% and 75% of cumulative energy) as
    well as other energy thresholds between 10% and 95%.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``, and
    ``v_s30``) can be a scalar or an array, and the arrays are broadcast
    against each other. For a scalar scenario, ``duration``,
    ``duration_plus_sigma``, ``duration_minus_sigma``, and the values returned
    by :meth:`duration_for_energy` are floats. For array scenarios, they are
    arrays with the broadcast shape of the scenario values, e.g., (N,) for
    arrays of N scenarios.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> from pygmm.ground_motion.pinilla_ramos_et_al_2023 import PinillaRamosEtAl2023
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0, 7.5]), dist_rup=np.array([10.0, 50.0, 150.0]),
    ...     v_s30=400.0)
    >>> m = PinillaRamosEtAl2023(s)
    >>> m.duration.shape
    (3,)
    >>> median, plus_sigma, minus_sigma = m.duration_for_energy(0.95)
    >>> median.shape
    (3,)
    >>> s = pygmm.Scenario(mag=7.0, dist_rup=50.0, v_s30=400.0)
    >>> isinstance(PinillaRamosEtAl2023(s).duration, float)
    True

    """

    NAME = "Pinilla-Ramos et al. (2023)"
    ABBREV = "PR23"

    PARAMS = [
        model.NumericParameter("mag", True, 4.0, 8.0),
        model.NumericParameter("dist_rup", True, 0, 300),
        model.NumericParameter("v_s30", True, 150, 2000),
    ]

    #: Conditional model coefficients for different energy thresholds
    CONDITIONAL_MODEL = _conditional_model()

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)

        # Calculate the duration for D5-75 by default
        self._duration, self._duration_plus_sigma, self._duration_minus_sigma = (
            self._calc_duration(0.75)
        )

    def _calc_duration(
        self, energy: float = 0.75
    ) -> tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray]:
        """Calculate duration for specified energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold (0.75 for D5-75, 0.95 for D5-95, etc.)

        Returns
        -------
        tuple
            (median_duration, duration_plus_sigma, duration_minus_sigma), as
            floats for a scalar scenario or arrays with the broadcast shape of
            the scenario values
        """
        conditional_model = self.CONDITIONAL_MODEL

        # Get scenario parameters. Scalars stay 0-d, and arrays are broadcast
        # against each other.
        s = self._scenario
        mag, rrup, vs30 = np.broadcast_arrays(
            np.asarray(s.mag, dtype=float),
            np.asarray(s.dist_rup, dtype=float),
            np.asarray(s.v_s30, dtype=float),
        )

        # Model coefficients
        n2, n3 = 0.3, 0.3

        # Magnitude scaling coefficients
        c1 = 3.655
        c21, c22, c23, c24, c2base = 0.41, 0.455, 0.54, 0.575, 0.515
        r1_c2, r2_c2, r3_c2 = 10, 42, 200

        # Distance scaling coefficients
        c31, c312, c32 = 0.063, 0.034, 0.083
        c3_base = 0.041

        # Site scaling coefficients
        c4 = -0.619
        phi, phi_max = 0.565, 1.11
        v1, v2, v3 = 200, 275, 2000
        r1, r2 = 44, 130
        s1 = 0.278
        mth = 6.75

        # Calculate distance-dependent magnitude scaling coefficient c2
        c2 = np.select(
            [rrup <= r1_c2, rrup <= r2_c2],
            [
                c21 + (c22 - c21) * rrup / r1_c2,
                c22 + (c23 - c22) * (rrup - r1_c2) / (r2_c2 - r1_c2),
            ],
            c23 + (c24 - c23) * (rrup - r2_c2) / (r3_c2 - r2_c2),
        )

        # Source term
        source_term = np.where(
            mag < mth,
            c1 * 10 ** ((mag - 6.75) * c2base),
            c1 * 10 ** ((mag - 6.75) * c2),
        )

        # Path term
        path_term = c3_base * rrup + np.select(
            [rrup < r1, rrup < r2],
            [c31 * rrup, c31 * r1 + c312 * (rrup - r1)],
            c31 * r1 + c312 * (r2 - r1) + c32 * (rrup - r2),
        )

        # Site term - VS30 scaling
        sigma_vs30 = np.select(
            [vs30 < v1, vs30 < v2, vs30 < v3 * 0.95, vs30 < v3],
            [
                phi_max,
                phi
                + (phi_max - phi)
                * (np.log(v2) - np.log(vs30))
                / (np.log(v2) - np.log(v1)),
                phi,
                phi * (np.log(v3) - np.log(vs30)) / (np.log(v3) - np.log(v3 * 0.95)),
            ],
            0.0,
        )

        site_term = c4 * np.log(vs30 / v3) * np.exp(s1 * sigma_vs30)

        # Median duration
        median = source_term + path_term + site_term

        # Standard deviation coefficients
        a0, a1, a2, b1, b2, d1, d2, d3, v4 = (
            0.537,
            -0.093,
            0.0278,
            -0.0372,
            0.00179,
            0.0206,
            2.401,
            0.0419,
            200,
        )

        sigma_vs30_std = np.minimum(d1 * (v4 / vs30) ** d2, d3)
        sigma = (
            a0
            + a1 * (rrup / 100)
            + a2 * (rrup / 100) ** 2
            + b1 * mag
            + b2 * mag**2
            + sigma_vs30_std
        )

        if np.isclose(energy, 0.75):
            # Use base model for D5-75
            model_plus_sigma = ((median) ** n2 + sigma) ** (1 / n2)
            model_minus_sigma = ((median) ** n2 - sigma) ** (1 / n2)

            return (
                _as_output(median),
                _as_output(model_plus_sigma),
                _as_output(model_minus_sigma),
            )

        else:  # energy != 0.75 - use conditional model
            # Find closest energy threshold
            mask = np.argmin(np.abs(conditional_model["energy"] - energy))
            c_median, a0x, r1x, v1x, sigma_c, rho_c = conditional_model[
                ["cmed", "a0", "r1x", "v1x", "std_c", "r"]
            ][mask]

            # Conditional model
            c_model = c_median + a0x + r1x * rrup + v1x * np.log(vs30 / 2000)

            df_dalpha = c_model**n3
            df_dc = median**n2 * n3 * c_model ** (n3 - 1)
            sigma_cond = np.sqrt(
                (df_dalpha * sigma) ** 2
                + (df_dc * sigma_c) ** 2
                + 2 * rho_c * df_dalpha * df_dc * sigma * sigma_c
            )

            model_plus_sigma = ((c_model * median) ** n3 + sigma_cond) ** (1 / n3)
            model_minus_sigma = ((c_model * median) ** n3 - sigma_cond) ** (1 / n3)

            return (
                _as_output(c_model * median),
                _as_output(model_plus_sigma),
                _as_output(model_minus_sigma),
            )

    @property
    def duration(self) -> float | np.ndarray:
        """Duration in seconds."""
        return self._duration

    @property
    def duration_plus_sigma(self) -> float | np.ndarray:
        """Duration plus one standard deviation in seconds."""
        return self._duration_plus_sigma

    @property
    def duration_minus_sigma(self) -> float | np.ndarray:
        """Duration minus one standard deviation in seconds."""
        return self._duration_minus_sigma

    def duration_for_energy(
        self, energy: float
    ) -> tuple[float | np.ndarray, float | np.ndarray, float | np.ndarray]:
        """Calculate duration for specified energy threshold.

        Parameters
        ----------
        energy : float
            Energy threshold. For D5-95 use 0.95, for D5-45 use 0.45, etc.
            Valid range is approximately 0.10 to 0.95.

        Returns
        -------
        tuple
            (median_duration, duration_plus_sigma, duration_minus_sigma) in
            seconds, as floats for a scalar scenario or arrays with the
            broadcast shape of the scenario values

        Raises
        ------
        ValueError
            If energy threshold is outside valid range
        """
        if energy < 0.10 or energy > 0.95:
            raise ValueError(
                f"Energy threshold {energy} is outside valid range [0.10, 0.95]"
            )

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
