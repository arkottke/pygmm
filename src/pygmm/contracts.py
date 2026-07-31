"""Data contracts for empirical-model outputs.

Frozen dataclasses with field names that mirror the attribute names already
used by pygmm models, so concrete models satisfy the contracts natively with
no adapter code. Consumers (pyrvt, pystrata) duck-type against these shapes;
they do not depend on this module at runtime — each maintains a private
duplicate of the dataclasses it needs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt


@dataclass(frozen=True)
class ResponseSpectrum:
    """Pseudo-spectral acceleration response spectrum.

    Attributes
    ----------
    periods : np.ndarray
        Oscillator periods [s].
    spec_accels : np.ndarray
        Pseudo-spectral accelerations [g].
    damping : float
        Oscillator damping ratio (decimal, e.g. 0.05).
    duration : float, optional
        Ground-motion duration [s] (e.g. for RVT consumers).
    """

    periods: npt.NDArray[np.floating]
    spec_accels: npt.NDArray[np.floating]
    damping: float
    duration: float | None = None


@dataclass(frozen=True)
class FourierSpectrum:
    """One-sided Fourier amplitude spectrum of acceleration.

    Attributes
    ----------
    freqs : np.ndarray
        Frequencies [Hz].
    fourier_amps : np.ndarray
        Fourier amplitudes of acceleration [g-s].
    duration : float
        Ground-motion duration [s] used by RVT.
    """

    freqs: npt.NDArray[np.floating]
    fourier_amps: npt.NDArray[np.floating]
    duration: float


@dataclass(frozen=True)
class Duration:
    """Ground-motion duration with optional named-pair variants.

    Attributes
    ----------
    duration : float
        Primary duration [s] (e.g. 5–95% significant duration).
    d_5_75 : float, optional
        5–75% significant duration [s].
    d_5_95 : float, optional
        5–95% significant duration [s].
    d_20_80 : float, optional
        20–80% significant duration [s].
    """

    duration: float
    d_5_75: float | None = None
    d_5_95: float | None = None
    d_20_80: float | None = None


@dataclass(frozen=True)
class NonlinearSoilCurves:
    """Strain-dependent modulus-reduction and damping curves.

    Attributes
    ----------
    strains : np.ndarray
        Shear strains (decimal, e.g. 1e-4).
    mod_reduc : np.ndarray
        Shear-modulus reduction ratio G/G_max.
    damping : np.ndarray
        Material damping ratio (decimal).
    damping_min : float
        Small-strain (minimum) damping (decimal).
    unit_wt : float, optional
        Soil unit weight [kN/m^3].
    name : str, optional
        Display name for the curve set.
    """

    strains: npt.NDArray[np.floating]
    mod_reduc: npt.NDArray[np.floating]
    damping: npt.NDArray[np.floating]
    damping_min: float
    unit_wt: float | None = None
    name: str | None = None


@dataclass(frozen=True)
class VelocityProfile:
    """Shear-wave velocity profile.

    Attributes
    ----------
    depth : np.ndarray
        Depths to layer tops [m].
    vs_median : np.ndarray
        Median shear-wave velocity at each depth [m/s].
    std_vs_ln : np.ndarray
        Natural-log standard deviation of Vs.
    region : str, optional
        Geographic / tectonic region descriptor.
    site_class : str, optional
        Site classification (e.g. NEHRP letter).
    """

    depth: npt.NDArray[np.floating]
    vs_median: npt.NDArray[np.floating]
    std_vs_ln: npt.NDArray[np.floating]
    region: str | None = None
    site_class: str | None = None


@dataclass(frozen=True)
class FaultDisplacement:
    """Surface fault-displacement prediction.

    Attributes
    ----------
    mag : float
        Moment magnitude.
    displacement_mean : float
        Mean displacement [m] (in natural-log space if `sigma_ln_disp` given).
    sigma_ln_disp : float
        Natural-log standard deviation of displacement.
    dist_from_rupture : float, optional
        Along-strike distance from the rupture endpoint [km].
    position_ratio : float, optional
        Normalized along-strike position (0 = endpoint, 0.5 = center).
    """

    mag: float
    displacement_mean: float
    sigma_ln_disp: float
    dist_from_rupture: float | None = None
    position_ratio: float | None = None


@dataclass(frozen=True)
class CptSounding:
    """Cone-penetration-test sounding (input contract).

    Attributes
    ----------
    depth : np.ndarray
        Measurement depths [m].
    q_c : np.ndarray
        Cone tip resistance [MPa].
    f_s : np.ndarray
        Sleeve friction [kPa].
    u_2 : np.ndarray, optional
        Pore-pressure measurement behind the tip [kPa].
    water_table_depth : float, optional
        Depth to ground-water table [m].
    unit_wts : np.ndarray, optional
        Per-depth unit weights [kN/m^3].
    """

    depth: npt.NDArray[np.floating]
    q_c: npt.NDArray[np.floating]
    f_s: npt.NDArray[np.floating]
    u_2: npt.NDArray[np.floating] | None = None
    water_table_depth: float | None = None
    unit_wts: npt.NDArray[np.floating] | None = None


@dataclass(frozen=True)
class SoilBehaviorProfile:
    """CPT-derived soil-behavior-type profile.

    Attributes
    ----------
    depth : np.ndarray
        Depths [m].
    ic : np.ndarray
        Soil-behavior-type index Ic.
    sbt_class : np.ndarray
        Soil-behavior-type classification (integer codes).
    fines_content : np.ndarray, optional
        Estimated fines content (decimal).
    """

    depth: npt.NDArray[np.floating]
    ic: npt.NDArray[np.floating]
    sbt_class: npt.NDArray[np.integer]
    fines_content: npt.NDArray[np.floating] | None = None


@dataclass(frozen=True)
class LiquefactionTriggering:
    """CPT-based liquefaction-triggering evaluation (future use).

    Attributes
    ----------
    depth : np.ndarray
        Depths [m].
    csr : np.ndarray
        Cyclic stress ratio.
    crr : np.ndarray
        Cyclic resistance ratio.
    factor_of_safety : np.ndarray
        FS_liq = crr / csr.
    prob_liquefaction : np.ndarray, optional
        Probability of liquefaction (decimal).
    """

    depth: npt.NDArray[np.floating]
    csr: npt.NDArray[np.floating]
    crr: npt.NDArray[np.floating]
    factor_of_safety: npt.NDArray[np.floating]
    prob_liquefaction: npt.NDArray[np.floating] | None = None


__all__ = [
    "ResponseSpectrum",
    "FourierSpectrum",
    "Duration",
    "NonlinearSoilCurves",
    "VelocityProfile",
    "FaultDisplacement",
    "CptSounding",
    "SoilBehaviorProfile",
    "LiquefactionTriggering",
]
