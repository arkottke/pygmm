from __future__ import annotations

from ..contracts import ResponseSpectrum
from ..model import GroundMotionModel


class ResponseSpectrumModel(GroundMotionModel):
    """Abstract base class for response-spectrum ground motion models."""

    def response_spectrum(self) -> ResponseSpectrum:
        """Return the computed PSA as a :class:`~pygmm.contracts.ResponseSpectrum`."""
        from .. import contracts

        return contracts.ResponseSpectrum(
            periods=self.periods,
            spec_accels=self.spec_accels,
            damping=0.05,
        )
