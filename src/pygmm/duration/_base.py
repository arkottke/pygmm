from __future__ import annotations

from abc import abstractmethod

from ..model import Model
from ..contracts import Duration


class DurationModel(Model):
    """Abstract base class for ground-motion duration models."""

    @abstractmethod
    def duration_model(self) -> Duration:
        """Return the computed duration as a :class:`~pygmm.contracts.Duration`."""
        ...
