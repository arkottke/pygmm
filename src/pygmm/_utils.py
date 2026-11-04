"""Helpers shared by the function-style (non-``Model``) relationships.

These relationships follow a single convention: the decorated function returns
``(ln_mean, ln_std)`` and :func:`dist_lognorm` turns that into the median,
unless the caller passes ``stats=True``.
"""

from __future__ import annotations

import functools
from collections.abc import Callable, Iterable
from typing import Any

import numpy as np
import scipy.constants as SC

#: Unit weight of water [kN/m³]
UNIT_WT_WATER = SC.g
#: Atmospheric pressure [kPa]
PRESS_ATM = SC.atm / SC.kilo
#: Multiply by this to convert kPa to atm
KPA_TO_ATM = SC.kilo / SC.atm


def check_bounds(
    value: Any, lower: float | None, upper: float | None, name: str, format: str = "f"
) -> None:
    """Raise :class:`ValueError` if any of `value` is outside [`lower`, `upper`]."""
    if lower is not None and np.any(np.asarray(value) < lower):
        found = np.min(value)
        raise ValueError(
            f'"{name}" value of {found:{format}} below lower limit of {lower:{format}}'
        )

    if upper is not None and np.any(np.asarray(value) > upper):
        found = np.max(value)
        raise ValueError(
            f'"{name}" value of {found:{format}} above upper limit of {upper:{format}}'
        )


def check_options(value: Any, options: Iterable[Any], name: str) -> None:
    """Raise :class:`ValueError` if `value` is not one of `options`."""
    if isinstance(options, str):
        # Each character is an option; "bc" must not match "bcd"
        options = tuple(options)
    if value not in options:
        raise ValueError(f'"{name}" not in options: {options}')


def dist_lognorm(f: Callable[..., tuple[Any, Any]]) -> Callable[..., Any]:
    """Return the median of a lognormal relationship, or its statistics.

    The wrapped function returns ``(ln_mean, ln_std)``. The wrapper returns
    ``exp(ln_mean)``, or the ``(ln_mean, ln_std)`` pair when called with the
    keyword ``stats=True``.
    """

    @functools.wraps(f)
    def wrapper(*args, stats: bool = False, **kwds):
        ln_mean, ln_std = f(*args, **kwds)

        if stats:
            return ln_mean, ln_std
        else:
            return np.exp(ln_mean)

    return wrapper
