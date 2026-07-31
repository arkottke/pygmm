from __future__ import annotations

from abc import abstractmethod

from ..contracts import FourierSpectrum
from ..model import Model


class FourierSpectrumModel(Model):
    """Abstract base class for Fourier amplitude spectrum models."""

    @abstractmethod
    def fourier_spectrum(self) -> FourierSpectrum:
        """Return the computed FAS as a :class:`~pygmm.contracts.FourierSpectrum`."""
        ...
