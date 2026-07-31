from __future__ import annotations

from abc import abstractmethod

from ..contracts import Duration
from ..model import Model


class DurationModel(Model):
    """Abstract base class for ground-motion duration models."""

    @abstractmethod
    def duration_model(self) -> Duration:
        """Return the computed duration as a :class:`~pygmm.contracts.Duration`."""
        ...
