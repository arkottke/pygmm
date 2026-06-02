"""Single-corner source-theory FAS model (Campbell 2003 defaults).

Ported from ``pyrvt.motions.SourceTheoryMotion``. This class produces only the
Fourier amplitude spectrum (g-sec) and a duration estimate; it has no RVT peak
calculator. Combine with ``pyrvt.motions.RvtMotion.from_fas`` to evaluate peak
ground motions.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt
from scipy.constants import g as gravity
from scipy.interpolate import interp1d

from ._source_helpers import (
    calc_geometric_spreading,
    calc_stress_drop,
    log_spaced_values,
)

# Region aliases mirror pyrvt.peak_calculators.get_region behaviour at the
# values used by SourceTheoryMotion.
_REGION_ALIASES = {
    "wna": "wna",
    "wus": "wna",
    "cena": "cena",
    "ceus": "cena",
}


def _normalize_region(region: str) -> str:
    try:
        return _REGION_ALIASES[region.lower()]
    except KeyError as exc:
        raise ValueError(f"Unknown region: {region!r}") from exc


class SourceTheoryModel:
    """Single-corner Brune source-spectrum model.

    Default crustal and path parameters follow Campbell (2003) for WNA and
    CENA. Produces an acceleration FAS in g-sec satisfying the
    :class:`pygmm.contracts.FourierSpectrum` shape.

    Parameters
    ----------
    magnitude : float
        Moment magnitude.
    distance : float
        Epicentral distance [km].
    region : {"wna", "wus", "cena", "ceus"}
        Region for default crustal / path parameters.
    stress_drop : float, optional
        Stress drop [bars]. Region default used if ``None``.
    depth : float, optional
        Hypocentral depth [km]. Default 8 km.
    freqs : array_like, optional
        Frequencies [Hz]. Default ``log_spaced_values(0.05, 200)``.
    disable_site_amp : bool, optional
        Disable crustal site amplification.
    """

    NAME = "Single-corner Source Theory"
    ABBREV = "ST"

    def __init__(
        self,
        magnitude: float,
        distance: float,
        region: str,
        stress_drop: float | None = None,
        depth: float = 8.0,
        freqs: npt.ArrayLike | None = None,
        disable_site_amp: bool = False,
    ):
        self.magnitude = float(magnitude)
        self.distance = float(distance)
        self.region = _normalize_region(region)
        self.depth = float(depth)
        self._disable_site_amp = bool(disable_site_amp)

        self._set_region_defaults(stress_drop)

        self.hypo_distance = float(np.hypot(self.distance, self.depth))
        self.seismic_moment = 10.0 ** (1.5 * (self.magnitude + 10.7))
        self.corner_freq = (
            4.9e6
            * self.shear_velocity
            * (self.stress_drop / self.seismic_moment) ** (1.0 / 3.0)
        )

        if freqs is None:
            self._freqs = log_spaced_values(0.05, 200.0)
        else:
            arr = np.asarray(freqs, dtype=float)
            self._freqs = arr if arr[0] <= arr[-1] else arr[::-1]

        self._fourier_amps = self._calc_fourier_amps(self._freqs)
        self._duration = self._calc_duration()

    @property
    def freqs(self) -> np.ndarray:
        """Frequencies [Hz]."""
        return self._freqs

    @property
    def fourier_amps(self) -> np.ndarray:
        """Acceleration Fourier amplitudes [g-sec]."""
        return self._fourier_amps

    @property
    def duration(self) -> float:
        """Source + path duration [sec]."""
        return self._duration

    def _set_region_defaults(self, stress_drop: float | None) -> None:
        if self.region == "wna":
            self.shear_velocity = 3.5
            self.density = 2.8
            self.path_atten_coeff = 180.0
            self.path_atten_power = 0.45
            self.site_atten = 0.04
            self.geometric_spreading = [(1, 40), (0.5, None)]
            self.stress_drop = float(stress_drop) if stress_drop else 100.0
            self.site_amp = interp1d(
                np.log(
                    [0.01, 0.09, 0.16, 0.51, 0.84, 1.25, 2.26, 3.17, 6.05,
                     16.60, 61.20, 100.00]
                ),
                [1.00, 1.10, 1.18, 1.42, 1.58, 1.74, 2.06, 2.25, 2.58, 3.13,
                 4.00, 4.40],
                bounds_error=False,
            )
        elif self.region == "cena":
            self.shear_velocity = 3.6
            self.density = 2.8
            self.path_atten_coeff = 680.0
            self.path_atten_power = 0.36
            self.site_atten = 0.006
            self.geometric_spreading = [(1, 70), (0, 130), (0.5, None)]
            self.stress_drop = (
                float(stress_drop) if stress_drop else calc_stress_drop(self.magnitude)
            )
            self.site_amp = interp1d(
                np.log(
                    [0.01, 0.10, 0.20, 0.30, 0.50, 0.90, 1.25, 1.80, 3.00,
                     5.30, 8.00, 14.00, 30.00, 60.00, 100.00]
                ),
                [1.00, 1.02, 1.03, 1.05, 1.07, 1.09, 1.11, 1.12, 1.13, 1.14,
                 1.15, 1.15, 1.15, 1.15, 1.15],
                bounds_error=False,
                fill_value=(1.0, 1.15),
            )
        else:  # pragma: no cover - guarded by _normalize_region
            raise NotImplementedError

    def _calc_duration(self) -> float:
        duration_source = 1.0 / self.corner_freq
        if self.region == "wna":
            duration_path = 0.05 * self.hypo_distance
        else:  # cena
            duration_path = 0.0
            if self.hypo_distance > 10:
                duration_path += 0.16 * (min(self.hypo_distance, 70) - 10.0)
            if self.hypo_distance > 70:
                duration_path += -0.03 * (min(self.hypo_distance, 130) - 70.0)
            if self.hypo_distance > 130:
                duration_path += 0.04 * (self.hypo_distance - 130.0)
        return duration_source + duration_path

    def _calc_fourier_amps(self, freqs: np.ndarray) -> np.ndarray:
        const = (0.55 * 2.0) / (
            np.sqrt(2.0) * 4.0 * np.pi * self.density * self.shear_velocity**3.0
        )
        source_comp = const * self.seismic_moment / (
            1.0 + (freqs / self.corner_freq) ** 2.0
        )

        path_atten = self.path_atten_coeff * freqs**self.path_atten_power
        geo_atten = calc_geometric_spreading(
            self.hypo_distance, self.geometric_spreading
        )
        path_comp = geo_atten * np.exp(
            -np.pi * freqs * self.hypo_distance / (path_atten * self.shear_velocity)
        )

        site_dim = np.exp(-np.pi * self.site_atten * freqs)
        ln_freqs = np.log(freqs)
        site_amp = self.site_amp(ln_freqs)
        if np.any(np.isnan(site_amp)):
            mask = ln_freqs < self.site_amp.x[0]
            site_amp[mask] = self.site_amp.y[0]
            mask = self.site_amp.x[-1] < ln_freqs
            site_amp[mask] = self.site_amp.y[-1]
        site_comp = 1.0 if self._disable_site_amp else (site_amp * site_dim)

        # Convert dyne-cm to g-sec.
        conv = 1.0e-20 / (100 * gravity)
        return (
            conv * (2.0 * np.pi * freqs) ** 2.0 * source_comp * path_comp * site_comp
        )
