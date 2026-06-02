"""Shared helpers for FAS source-theory models.

Inlined from the original pyrvt implementations so that ``pygmm.fourier_spectrum``
has no dependency on pyrvt at import or runtime.
"""

from __future__ import annotations

import numpy as np


def log_spaced_values(lower: float, upper: float, per_decade: int = 512) -> np.ndarray:
    """Log-spaced values between ``lower`` and ``upper``."""
    lo = np.log10(lower)
    hi = np.log10(upper)
    count = int(np.ceil(per_decade * (hi - lo)))
    return np.logspace(lo, hi, num=count)


def calc_stress_drop(magnitude: float) -> float:
    """Atkinson & Boore (2011) stress drop [bars]."""
    return 10 ** (3.45 - 0.2 * max(magnitude, 5.0))


def calc_geometric_spreading(
    dist: float, params: list[tuple[float, float | None]]
) -> float:
    """Piece-wise linear geometric-spreading model."""
    initial = 1.0
    coeff = 1.0
    for slope, limit in params:
        d = min(dist, limit) if limit else dist
        coeff *= (initial / d) ** slope
        if d < dist:
            initial = d
        else:
            break
    return coeff
