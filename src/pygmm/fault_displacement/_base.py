"""Scaffold: Fault displacement prediction models."""

from __future__ import annotations

from abc import abstractmethod

from ..model import Model
from ..contracts import FaultDisplacement


class FaultDisplacementModel(Model):
    """Abstract base class for fault-displacement prediction models.

    TODO: No models implemented yet.
    """

    @abstractmethod
    def fault_displacement(self) -> FaultDisplacement:
        """Return the predicted fault displacement."""
        ...
