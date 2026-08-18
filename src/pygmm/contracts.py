"""Data contracts for empirical-model outputs.

Frozen dataclasses with field names that mirror the attribute names already
used by pygmm models, so concrete models satisfy the contracts natively with
no adapter code. Consumers (pyrvt, pystrata) duck-type against these shapes;
they do not depend on this module at runtime — each maintains a private
duplicate of the dataclasses it needs.

Because that mechanism is structural, the producer side is expressed with
:class:`typing.Protocol` rather than abstract base classes. A model satisfies
as many protocols as it structurally matches, with no inheritance edge and no
MRO involvement — which is what makes a single
:class:`~pygmm.model.GroundMotionModel` able to satisfy both
:class:`SupportsResponseSpectrum` and the peak-parameter accessors at once.
An ABC could enforce nothing across the package boundary anyway, which is why
six of the seven former ``_base.py`` classes never acquired a subclass.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

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

    Two mutually exclusive ways of expressing uncertainty are carried, because
    the published models genuinely disagree. Abrahamson & Silva (1996),
    Afshari & Stewart (2016) and Kempton & Stewart (2006) report a lognormal
    standard error, so ``ln_std`` applies. Pinilla-Ramos et al. (2023, 2024)
    apply sigma in a transformed power space --
    ``(median**n ± sigma)**(1/n)`` -- which is not symmetric in log space, so
    they report ``plus_sigma``/``minus_sigma`` directly and leave ``ln_std``
    unset. Collapsing these into one field would misstate one of them.

    Attributes
    ----------
    duration : float
        Primary duration [s].
    d_5_75, d_5_95, d_20_80 : float, optional
        Significant duration [s] over the named energy interval.
    ln_std : float, optional
        Natural-log standard deviation of ``duration``, for models whose
        uncertainty is lognormal.
    plus_sigma, minus_sigma : float, optional
        84th- and 16th-percentile duration [s], for models whose uncertainty
        is not lognormal.
    """

    duration: float
    d_5_75: float | None = None
    d_5_95: float | None = None
    d_20_80: float | None = None
    ln_std: float | None = None
    plus_sigma: float | None = None
    minus_sigma: float | None = None


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


# ---------------------------------------------------------------------------
# Producer protocols
#
# Structural, not nominal: a model satisfies these by having the right
# methods, without importing or inheriting anything from this module. That is
# how pystrata and pyrvt already consume pygmm, and it lets one class satisfy
# several protocols at once.
# ---------------------------------------------------------------------------


@runtime_checkable
class SupportsResponseSpectrum(Protocol):
    """A model that can emit a :class:`ResponseSpectrum`."""

    def response_spectrum(self, damping: float = 0.05) -> ResponseSpectrum: ...


@runtime_checkable
class SupportsFourierSpectrum(Protocol):
    """A model that can emit a :class:`FourierSpectrum`.

    ``duration`` is a parameter rather than a required attribute because
    ``BaylessAbrahamson2019`` predicts EAS only and has no intrinsic duration.
    Making it a parameter states that requirement in the signature instead of
    letting a ``None`` surface deep inside a consumer's RVT integration.
    """

    def fourier_spectrum(self, duration: float | None = None) -> FourierSpectrum: ...


@runtime_checkable
class SupportsDuration(Protocol):
    """A model that can emit a :class:`Duration`."""

    def duration_model(self) -> Duration: ...


@runtime_checkable
class SupportsSoilCurves(Protocol):
    """A model that can emit :class:`NonlinearSoilCurves`."""

    def curves(self) -> NonlinearSoilCurves: ...


__all__ = [
    "Duration",
    "FourierSpectrum",
    "NonlinearSoilCurves",
    "ResponseSpectrum",
    "SupportsDuration",
    "SupportsFourierSpectrum",
    "SupportsResponseSpectrum",
    "SupportsSoilCurves",
    "VelocityProfile",
]
