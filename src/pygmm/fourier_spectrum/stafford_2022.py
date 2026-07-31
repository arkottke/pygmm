"""Stafford et al. (2022) Fourier amplitude spectrum model.

Ported from ``pyrvt.motions.StaffordEtAl22Motion``. Produces only the FAS and
the duration (Boore & Thompson 2014). Combine with
``pyrvt.motions.RvtMotion.from_fas`` for peak-motion evaluation.
"""

from __future__ import annotations

import gzip
import pathlib

import numpy as np
import numpy.typing as npt
from scipy.constants import g as gravity
from scipy.interpolate import interp1d

from ._source_helpers import calc_geometric_spreading

_DATA_PATH = pathlib.Path(__file__).parent.parent / "data" / "sea22-site_amp.csv.gz"


class StaffordEtAl2022:
    """Stafford et al. (2022) point-source FAS model for Vs30 = 760 m/s.

    Parameters mirror ``pyrvt.motions.StaffordEtAl22Motion``.
    """

    NAME = "Stafford et al. (2022)"
    ABBREV = "Sea22"

    _ln_site_amp_interpolator: interp1d | None = None

    def __init__(
        self,
        mag: float,
        dist_rup: float | None = None,
        dist_jb: float | None = None,
        mechanism: str = "U",
        method: str = "continuous",
        delta_ztor: float = 0.0,
        freqs: npt.ArrayLike | None = None,
        disable_site_amp: bool = False,
    ):
        if dist_rup is None and dist_jb is None:
            raise ValueError("Either dist_rup or dist_jb must be provided.")
        if dist_rup is None:
            depth_tor = self.calc_depth_tor(mag, mechanism)
            dist_rup = float(np.hypot(depth_tor, dist_jb))
        self.mag = float(mag)
        self.dist_rup = float(dist_rup)
        self.mechanism = mechanism
        self.method = method

        if freqs is None:
            self._freqs = np.geomspace(0.05, 200, 512)
        else:
            self._freqs = np.asarray(freqs, dtype=float)

        shear_vel = 3.5
        density = 2.75
        site_atten = 0.039

        if method == "continuous":
            ln_ds0 = 4.599
            dln_ds0 = 0.4624
            dz_a = 0.0453
            dz_b = 0.109
            y_1 = 1.1611
            y_f = 0.5
            r_t = 50
            h_a, h_b, h_c, h_d, h_e = -0.8712, 0.4451, 1.1513, 5.0948, 7.2725
            Q_0 = 205.4
            n_a = 0.6884
            n_b = 0.1354
            n_c = 5.1278
        elif method == "trilinear":
            ln_ds0 = 5.07
            dln_ds0 = 0.6451
            dz_a = 0.4077
            dz_b = 0.117
            y_1 = 1.1680
            y_2 = 0.9293
            y_f = 0.5
            h_a, h_b, h_c, h_d, h_e = -0.7771, 0.4768, 1.1513, 3.418, 7.088
            Q_0 = 183.7
            n = 0.7077
        else:
            raise NotImplementedError(method)

        stress_drop = np.exp(
            ln_ds0
            + dln_ds0 * np.minimum(self.mag - 5, 0)
            + delta_ztor * (dz_a + dz_b / np.cosh(2 * np.maximum(self.mag - 4.5, 0)))
        )

        const = (0.55 / np.sqrt(2) * 2) / (4 * np.pi * density * shear_vel**3) * 1e-20
        seismic_moment = 10 ** (1.5 * (self.mag + 10.7))
        corner_freq = 4.9058e6 * shear_vel * (stress_drop / seismic_moment) ** (1 / 3)
        source_comp = (const * seismic_moment) / (1 + (self._freqs / corner_freq) ** 2)

        fault_fact = np.exp(
            h_a
            + h_b * self.mag
            + ((h_b - h_c) / h_d) * np.log(1 + np.exp(-h_d * (self.mag - h_e)))
        )
        dist_ps = self.dist_rup + fault_fact
        if method == "continuous":
            geom_spread = np.exp(
                -y_1 * np.log(dist_ps)
                + (y_1 - y_f)
                / 2
                * np.log((self.dist_rup**2 + r_t**2) / (1**2 + r_t**2))
            )
            n = n_a + n_b * np.tanh(self.mag - n_c)
            dist_ae = self.dist_rup
        else:
            geom_spread = calc_geometric_spreading(
                dist_ps, [(y_1, 25), (y_2, 85), (y_f, None)]
            )
            dist_ae = dist_ps

        anelastic_atten = np.exp(
            -(np.pi * self._freqs ** (1 - n) * dist_ae) / (Q_0 * shear_vel)
        )
        path_comp = geom_spread * anelastic_atten

        conv = (2 * np.pi * self._freqs) ** 2 / (gravity * 100)
        site_tf = 1.0 if disable_site_amp else self.site_amp(self._freqs, site_atten)

        self._fourier_amps = conv * source_comp * path_comp * site_tf
        self._dist_ps = dist_ps
        self._duration = self.calc_duration(corner_freq, dist_ps)

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
        """Boore & Thompson (2014) duration [sec]."""
        return self._duration

    @property
    def dist_ps(self) -> float:
        """Equivalent point-source distance [km]."""
        return self._dist_ps

    @classmethod
    def site_amp(cls, freqs, site_atten):
        if cls._ln_site_amp_interpolator is None:
            data = np.genfromtxt(
                gzip.open(_DATA_PATH),
                delimiter=",",
                names=True,
                skip_header=1,
            ).view(np.recarray)
            ln_amp = np.log(data["site_amp"])
            cls._ln_site_amp_interpolator = interp1d(
                data["freq"],
                ln_amp,
                kind="linear",
                bounds_error=False,
                fill_value=(ln_amp[0], ln_amp[-1]),
            )
        ln_amp = cls._ln_site_amp_interpolator(freqs)
        return np.exp(ln_amp - np.pi * site_atten * freqs)

    @staticmethod
    def calc_depth_tor(mag: float, mechanism: str) -> float:
        """Chiou & Youngs (2014) depth-to-top-of-rupture model."""
        if mechanism == "RS":
            fact = 2.704 - 1.226 * max(mag - 5.849, 0)
        else:
            fact = 2.673 - 1.136 * max(mag - 4.970, 0)
        return max(fact, 0) ** 2

    @staticmethod
    def calc_duration(corner_freq: float, dist_ps: float) -> float:
        """Boore & Thompson (2014) duration model [sec]."""
        d_s = 1.0 / corner_freq
        DISTS = [0, 7, 45, 125, 175, 270]
        D_P = [0.0, 2.4, 8.4, 10.9, 17.4, 34.2]
        if dist_ps < DISTS[-1]:
            d_p = float(np.interp(dist_ps, DISTS, D_P))
        else:
            d_p = D_P[-1] + 0.156 * (dist_ps - DISTS[-1])
        return d_s + d_p
