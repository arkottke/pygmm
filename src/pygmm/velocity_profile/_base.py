from __future__ import annotations

from abc import ABC, abstractmethod

from ..contracts import VelocityProfile


class VelocityProfileModel(ABC):
    """Abstract base class for shear-wave velocity profile models."""

    @abstractmethod
    def profile(self) -> VelocityProfile:
        """Compute and return the velocity profile."""
        ...
