from __future__ import annotations

from abc import ABC, abstractmethod

from ..contracts import NonlinearSoilCurves


class SoilCurveModel(ABC):
    """Abstract base class for strain-dependent soil modulus and damping models."""

    @abstractmethod
    def curves(self) -> NonlinearSoilCurves:
        """Compute and return the nonlinear soil curves."""
        ...
