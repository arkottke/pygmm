"""Pint unit-conversion decorators for pygmm.soil_curves.

Mirrors the interface in pystrata.units so that soil-curve models accept
:class:`pint.Quantity` arguments regardless of which package they come from.
Plain numeric values pass through unchanged.
"""

from __future__ import annotations

import inspect
from functools import wraps

import pint


def convert_units(**unit_specs: str):
    """Decorator that converts :class:`pint.Quantity` arguments to expected units.

    For each named parameter in *unit_specs*, if the value is a
    :class:`pint.Quantity` it is converted to the target unit and its magnitude
    is extracted.  Plain numeric values pass through unchanged.  ``None`` is
    always passed through (for optional parameters).  Incompatible units raise
    :class:`pint.DimensionalityError`.
    """

    def decorator(func):
        sig = inspect.signature(func)

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            for name, target_unit in unit_specs.items():
                if name not in bound.arguments:
                    continue
                value = bound.arguments[name]
                if value is None:
                    continue
                if isinstance(value, pint.Quantity):
                    bound.arguments[name] = value.to(target_unit).magnitude

            return func(*bound.args, **bound.kwargs)

        return wrapper

    return decorator


def convert_kwds_units(**unit_specs: str):
    """Decorator that converts :class:`pint.Quantity` values inside ``**kwds``.

    Intended for ``__init__`` methods that accept ``**kwds`` where some values
    may be :class:`pint.Quantity` objects (e.g. ``WangSoilType``).
    """

    def decorator(func):
        sig = inspect.signature(func)

        @wraps(func)
        def wrapper(*args, **kwargs):
            bound = sig.bind(*args, **kwargs)
            bound.apply_defaults()

            for name, target_unit in unit_specs.items():
                if name in bound.arguments:
                    value = bound.arguments[name]
                    if isinstance(value, pint.Quantity):
                        bound.arguments[name] = value.to(target_unit).magnitude
                for arg_value in list(bound.arguments.values()):
                    if isinstance(arg_value, dict) and name in arg_value:
                        v = arg_value[name]
                        if isinstance(v, pint.Quantity):
                            arg_value[name] = v.to(target_unit).magnitude

            return func(*bound.args, **bound.kwargs)

        return wrapper

    return decorator
