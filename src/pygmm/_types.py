"""Common types used across the project."""

from typing import Literal

import numpy as np

ArrayLike = list[float] | np.ndarray

#: Interpolation kinds accepted by :func:`scipy.interpolate.interp1d`.
#: Spelled out rather than typed as ``str`` so that an invalid kind is caught
#: by the type checker instead of raising from inside scipy.
InterpKind = Literal[
    "linear",
    "nearest",
    "nearest-up",
    "zero",
    "slinear",
    "quadratic",
    "cubic",
    "previous",
    "next",
]
