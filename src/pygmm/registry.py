"""Capability-based registry of the models in pyGMM.

The package deliberately does not encode a model's capabilities in its
location on disk: a single :class:`~pygmm.model.GroundMotionModel` commonly
emits PSA, PGA, PGV *and* PGD from one underlying response vector, and a
directory can express only one classification. Capabilities live here instead,
where a model can declare as many as it has.

Register a model with :func:`register` and query with :func:`find_models`::

    >>> import pygmm
    >>> [i.abbrev for i in pygmm.find_models(provides="pgv")]  # doctest: +SKIP
    ['ASK14', 'BSSA14', 'CB14', 'CY14', ...]

This replaces the ``pygmm.models`` list, which mixed response-spectrum,
duration, FAS and correlation models with no way to tell them apart.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import dataclass

__all__ = [
    "ModelInfo",
    "find_models",
    "get_model",
    "register",
]

#: Capabilities a model may declare. ``psa``/``pga``/``pgv``/``pgd`` are
#: derived automatically for :class:`~pygmm.model.GroundMotionModel`
#: subclasses; everything else must be declared.
CAPABILITIES = frozenset(
    {
        "psa",
        "pga",
        "pgv",
        "pgd",
        "duration",
        "fas",
        "arias",
        "cav",
        "soil_curves",
        "vh_ratio",
        "correlation",
    }
)

#: How a model is constructed.
#:
#: ``scenario``
#:     Takes a :class:`~pygmm.model.Scenario` and nothing else.
#: ``conditional``
#:     Takes a ``Scenario`` *and* a conditioning intensity measure (or a
#:     backbone model that supplies one) -- ``AbrahamsonBhasin2020``,
#:     ``AbrahamsonShiYang2016``, ``MacedoAbrahamsonLiu2021``.
#: ``kwargs``
#:     Takes loose engineering parameters rather than a ``Scenario``.
#: ``stateless``
#:     A namespace of classmethods; never instantiated.
INPUTS = frozenset({"scenario", "conditional", "kwargs", "stateless"})


@dataclass(frozen=True)
class ModelInfo:
    """What the registry knows about one model."""

    #: The model class itself.
    cls: type
    #: Unique key. The class name -- deliberately *not* ``ABBREV``, which is
    #: not unique (``CampbellBozorgnia2014`` and ``CoppersmithBommer2014``
    #: both use ``"CB14"``).
    key: str
    #: Long name (``NAME``).
    name: str
    #: Short name (``ABBREV``). Display only; may collide with another model.
    abbrev: str
    #: Capabilities, a subset of :data:`CAPABILITIES`.
    provides: frozenset[str]
    #: Tectonic setting, or ``None`` where the publication does not state one
    #: (or takes it as a runtime argument, as ``SourceTheoryModel`` does).
    tectonic: str | None
    #: One of :data:`INPUTS`.
    input: str
    #: True when the model loads a large data file on import.
    lazy_data: bool = False


_REGISTRY: dict[str, ModelInfo] = {}


def _derive_provides(cls: type) -> frozenset[str]:
    """Read the intensity measures a :class:`GroundMotionModel` exposes.

    Derived rather than declared: a hand-written list can drift out of sync
    with ``INDEX_PGV``, a derived one cannot.

    Restricted to genuine :class:`~pygmm.model.GroundMotionModel` subclasses.
    Other classes reuse the same ``INDEX_*`` idiom to slice vectors that are
    *not* spectral accelerations -- ``GulerceAbrahamson2011`` indexes
    vertical-to-horizontal ratios that way -- so reading these attributes off
    an arbitrary class would claim capabilities it does not have.
    """
    from .model import GroundMotionModel

    if not (isinstance(cls, type) and issubclass(cls, GroundMotionModel)):
        return frozenset()

    caps = set()
    if getattr(cls, "INDICES_PSA", None) is not None and len(cls.INDICES_PSA):
        caps.add("psa")
    for attr, cap in (("INDEX_PGA", "pga"), ("INDEX_PGV", "pgv"), ("INDEX_PGD", "pgd")):
        if getattr(cls, attr, None) is not None:
            caps.add(cap)
    return frozenset(caps)


def register(
    *,
    provides: Iterable[str] | None = None,
    tectonic: str | None = None,
    input: str = "scenario",
    lazy_data: bool = False,
) -> Callable[[type], type]:
    """Add a model class to the registry.

    Parameters
    ----------
    provides : iterable of str, optional
        Capabilities beyond those derivable from the class. For a
        :class:`~pygmm.model.GroundMotionModel` the PSA/PGA/PGV/PGD entries
        are derived and need not be listed; for anything else ``provides`` is
        required.
    tectonic : str, optional
        Tectonic setting, where the publication states one.
    input : str
        One of :data:`INPUTS`.
    lazy_data : bool
        Set for models that load a large data file, so callers can skip them.
    """
    if input not in INPUTS:
        raise ValueError(f"input must be one of {sorted(INPUTS)}, got {input!r}")

    def decorator(cls: type) -> type:
        key = cls.__name__
        if key in _REGISTRY:
            raise ValueError(f"{key} is already registered")

        caps = set(_derive_provides(cls))
        if provides is not None:
            unknown = set(provides) - CAPABILITIES
            if unknown:
                raise ValueError(
                    f"{key}: unknown capabilities {sorted(unknown)}; "
                    f"expected a subset of {sorted(CAPABILITIES)}"
                )
            caps |= set(provides)
        if not caps:
            raise ValueError(
                f"{key} declares no capabilities. Non-GroundMotionModel classes "
                "must pass provides=(...)."
            )

        _REGISTRY[key] = ModelInfo(
            cls=cls,
            key=key,
            name=getattr(cls, "NAME", "") or key,
            abbrev=getattr(cls, "ABBREV", ""),
            provides=frozenset(caps),
            tectonic=tectonic,
            input=input,
            lazy_data=lazy_data,
        )
        return cls

    return decorator


def find_models(
    *,
    provides: str | Iterable[str] | None = None,
    tectonic: str | None = None,
    input: str | None = None,
    lazy_data: bool | None = None,
) -> list[ModelInfo]:
    """Return the registered models matching every supplied filter.

    ``provides`` may be a single capability or an iterable, in which case a
    model must provide *all* of them.

    Examples
    --------
    >>> import pygmm
    >>> psa = pygmm.find_models(provides="psa", input="scenario")
    >>> all("psa" in i.provides for i in psa)
    True
    """
    if isinstance(provides, str):
        wanted = {provides}
    elif provides is None:
        wanted = set()
    else:
        wanted = set(provides)

    unknown = wanted - CAPABILITIES
    if unknown:
        raise ValueError(
            f"unknown capabilities {sorted(unknown)}; "
            f"expected a subset of {sorted(CAPABILITIES)}"
        )

    return [
        info
        for info in sorted(_REGISTRY.values(), key=lambda i: i.key)
        if wanted <= info.provides
        and (tectonic is None or info.tectonic == tectonic)
        and (input is None or info.input == input)
        and (lazy_data is None or info.lazy_data == lazy_data)
    ]


def get_model(key: str) -> type:
    """Return a registered model class by its key (the class name)."""
    try:
        return _REGISTRY[key].cls
    except KeyError:
        raise KeyError(
            f"{key!r} is not a registered model. Known models: {sorted(_REGISTRY)}"
        ) from None
