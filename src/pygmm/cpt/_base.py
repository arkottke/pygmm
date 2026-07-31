"""Scaffold: CPT-based soil classification models (Ic, SBT, future liquefaction)."""

from __future__ import annotations

from abc import abstractmethod

from ..contracts import CptSounding, SoilBehaviorProfile
from ..model import Model


class CptModel(Model):
    """Abstract base class for CPT-based soil classification models.

    TODO: Implement Ic classification, then liquefaction triggering.
    """

    PARAMS = []

    @abstractmethod
    def classify(self, sounding: CptSounding) -> SoilBehaviorProfile:
        """Classify soil behavior type from a CPT sounding."""
        ...
