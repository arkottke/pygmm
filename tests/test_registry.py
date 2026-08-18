"""Integrity of the model registry.

These tests are the drift guard that ``pygmm.models`` never was. That list
was a hand-maintained grab-bag whose only reference lived in an uncollected
test, so it silently accumulated four inconsistencies with ``__all__`` and
listed seven models with no ``spec_accels``. The checks below would have
caught every one of those, plus ``AbrahamsonBhasin2020.__init__`` returning a
tuple instead of ``None``.
"""

import numpy as np
import pytest

import pygmm

# Imported for its registration side effect only. `ground_motion/__init__.py`
# deliberately does not import this module, so without this line the registry
# contents -- and therefore this file's parametrization -- would depend on
# whether some other test module happened to be collected first.
import pygmm.ground_motion.hermkes_kuehn_riggelsen_2014  # noqa: F401
from pygmm.model import GroundMotionModel
from pygmm.registry import CAPABILITIES, INPUTS

ALL_MODELS = pygmm.find_models()

#: A scenario broad enough to construct most registered models. Extra keys
#: are inert: each model filters to its own ``PARAMS``.
SCENARIO_KWDS = dict(
    depth_1_0=0.05,
    depth_2_5=5,
    depth_bor=15,
    depth_hyp=9,
    depth_tor=5,
    dip=90.0,
    dist_jb=30.0,
    dist_rup=30.0,
    dist_x=30.0,
    dpp_centered=0,
    event_type="interface",
    mag=6.5,
    mechanism="SS",
    on_hanging_wall=False,
    pga_ref=0.2,
    site_cond="soil",
    v_s30=500.0,
    width=10,
)

#: Scenario keys whose categorical options differ per model, so no single
#: value satisfies everything.
SCENARIO_OVERRIDES = {
    "PinillaRamosEtAl2024": dict(region="Japan"),
}

#: Conditioning arguments for ``input="conditional"`` models.
CONSTRUCTOR_KWDS = {
    "AbrahamsonBhasin2020": dict(pga=0.3),
    "AbrahamsonShiYang2016": dict(pga=0.3, sa_t1=0.2),
    "MacedoAbrahamsonLiu2021": dict(pga=0.3),
}


#: Models that cannot currently be constructed, with the reason. Recorded as
#: xfail rather than skipped so the breakage stays visible.
KNOWN_BROKEN = {
    "HermkesKuehnRiggelsen2014": (
        "data/hermkes_kuehn_riggelsen_2014.npz holds pickled objects, which "
        "numpy has refused to load without allow_pickle=True since 1.16.3. "
        "Pre-existing; the model's own test is marked slow and deselected by "
        "default, so this has been invisible."
    ),
}


def _build(info):
    """Construct a registered model with whatever it needs."""
    kwds = dict(SCENARIO_KWDS, **SCENARIO_OVERRIDES.get(info.key, {}))
    return info.cls(pygmm.Scenario(**kwds), **CONSTRUCTOR_KWDS.get(info.key, {}))


def _all_finite(value):
    """Finiteness check that also handles the structured arrays the duration
    models return (``AfshariStewart2016`` yields a record with ``D_5t75``,
    ``D_5t95``, ``D_20t80``; ``PinillaRamosEtAl2023`` yields a plain float)."""
    arr = np.asarray(value)
    if arr.dtype.names:
        return all(np.all(np.isfinite(arr[n])) for n in arr.dtype.names)
    return np.all(np.isfinite(arr.astype(float)))


#: Attribute each capability promises. ``correlation`` and ``soil_curves``
#: are exercised separately -- they are method-based, not attribute-based.
CAPABILITY_ATTRS = {
    "psa": ("periods", "spec_accels"),
    "pga": ("pga",),
    "pgv": ("pgv",),
    "pgd": ("pgd",),
    "duration": ("duration",),
    "fas": ("freqs", "fourier_amps"),
    "arias": ("ia",),
    "cav": ("cav",),
    "vh_ratio": ("ratio",),
}


def test_registry_is_not_empty():
    assert len(ALL_MODELS) > 25


def test_keys_are_unique():
    keys = [i.key for i in ALL_MODELS]
    assert len(keys) == len(set(keys))


def test_abbrev_is_not_used_as_a_key():
    """``ABBREV`` collides, which is why the registry keys on class name.

    ``CampbellBozorgnia2014`` and ``CoppersmithBommer2014`` both use "CB14".
    """
    abbrevs = [i.abbrev for i in ALL_MODELS if i.abbrev]
    assert len(abbrevs) != len(set(abbrevs)), (
        "ABBREV is now unique; the comment explaining why keys are class "
        "names should be revisited."
    )


@pytest.mark.parametrize("info", ALL_MODELS, ids=lambda i: i.key)
def test_metadata_is_well_formed(info):
    assert info.key == info.cls.__name__
    assert info.name
    assert info.input in INPUTS
    assert info.provides <= CAPABILITIES
    assert info.provides, "every model must declare at least one capability"


@pytest.mark.parametrize(
    "info",
    [i for i in ALL_MODELS if i.input in ("scenario", "conditional")],
    ids=lambda i: i.key,
)
def test_scenario_models_construct_and_honour_provides(info):
    """Build each model once and check every capability it claims.

    This is what makes a broken ``__init__`` impossible to ship unnoticed:
    ``AbrahamsonBhasin2020`` ended its ``__init__`` with
    ``return ln_mean, ln_std``, so merely constructing it raised ``TypeError``
    -- and nothing in the suite ever constructed it.
    """
    if info.key in KNOWN_BROKEN:
        pytest.xfail(KNOWN_BROKEN[info.key])

    obj = _build(info)

    for cap in sorted(info.provides):
        for attr in CAPABILITY_ATTRS.get(cap, ()):
            assert hasattr(obj, attr), f"claims {cap!r} but has no {attr!r}"
            value = getattr(obj, attr)
            assert value is not None, f"{cap!r}: {attr!r} is None"
            assert _all_finite(value), f"{cap!r}: {attr!r} is not finite"


@pytest.mark.parametrize(
    "info", [i for i in ALL_MODELS if "psa" in i.provides], ids=lambda i: i.key
)
def test_psa_models_are_ground_motion_models(info):
    """``psa`` is derived from ``INDEX_*``, which only means PSA on a GMM.

    ``GulerceAbrahamson2011`` uses the same idiom to index vertical-to-
    horizontal ratios, so deriving from the attributes alone would claim a
    capability it does not have.
    """
    assert issubclass(info.cls, GroundMotionModel)


@pytest.mark.parametrize(
    "info", [i for i in ALL_MODELS if i.input == "kwargs"], ids=lambda i: i.key
)
def test_kwargs_models_do_not_take_a_scenario(info):
    """``input="kwargs"`` models take loose engineering parameters."""
    assert not issubclass(info.cls, pygmm.model.Model)


def test_every_registered_model_is_publicly_reachable():
    """Registered models must be importable from ``pygmm`` or a subpackage.

    Catches a model that exists and registers but was never wired into any
    ``__all__`` -- the state ``HermkesKuehnRiggelsen2014`` was already in, and
    that ``AbrahamsonSilva1996`` was in with respect to ``pygmm.__all__``.
    """
    reachable = set(pygmm.__all__)
    for pkg in (pygmm.ground_motion, pygmm.fourier_spectrum, pygmm.soil_curves):
        reachable |= set(pkg.__all__)

    unreachable = {i.key for i in ALL_MODELS} - reachable
    # Exempt: registered on direct import only, since importing its module can
    # trigger a 430 kB download.
    assert unreachable <= {"HermkesKuehnRiggelsen2014"}, unreachable


def test_every_exported_model_class_is_registered():
    """The converse: nothing public escapes the registry."""
    exported = {
        n: getattr(pygmm, n)
        for n in pygmm.__all__
        if isinstance(getattr(pygmm, n, None), type)
    }
    model_like = {
        n for n, o in exported.items() if hasattr(o, "NAME") or n.endswith("SoilType")
    }
    assert model_like - {i.key for i in ALL_MODELS} == set()


def test_get_model_roundtrip():
    for info in ALL_MODELS:
        assert pygmm.get_model(info.key) is info.cls


def test_get_model_unknown_key():
    with pytest.raises(KeyError, match="not a registered model"):
        pygmm.get_model("NoSuchModel2099")


def test_find_models_rejects_unknown_capability():
    with pytest.raises(ValueError, match="unknown capabilities"):
        pygmm.find_models(provides="not_a_capability")


def test_find_models_requires_all_capabilities():
    """An iterable of capabilities is an AND, not an OR."""
    both = pygmm.find_models(provides=["psa", "pgv"])
    assert all({"psa", "pgv"} <= i.provides for i in both)
    assert len(both) < len(pygmm.find_models(provides="psa"))


def test_register_rejects_duplicate_key():
    with pytest.raises(ValueError, match="already registered"):
        pygmm.register(provides=("psa",))(pygmm.ChiouYoungs2014)


def test_register_rejects_unknown_capability():
    class Dummy:
        NAME = "Dummy"

    with pytest.raises(ValueError, match="unknown capabilities"):
        pygmm.register(provides=("not_a_capability",))(Dummy)


def test_register_requires_capabilities_for_plain_classes():
    class Dummy:
        NAME = "Dummy"

    with pytest.raises(ValueError, match="declares no capabilities"):
        pygmm.register()(Dummy)


def test_register_rejects_unknown_input():
    with pytest.raises(ValueError, match="input must be one of"):
        pygmm.register(provides=("psa",), input="telepathy")


def test_models_list_is_deprecated():
    with pytest.deprecated_call(match="find_models"):
        legacy = pygmm.models
    assert all(isinstance(m, type) for m in legacy)
    assert {m.__name__ for m in legacy} == {i.key for i in ALL_MODELS}


def test_unknown_module_attribute_still_raises():
    """The ``models`` shim must not swallow genuine typos."""
    with pytest.raises(AttributeError):
        pygmm.no_such_attribute
