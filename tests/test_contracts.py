"""Smoke tests for pygmm.contracts dataclasses.

These exist primarily to (1) guarantee the module imports cleanly and
(2) lock in field names so they cannot drift away from the consumer-side
duplicates in pyrvt and pystrata without breaking a test.
"""

from dataclasses import fields

import numpy as np
import pytest

import pygmm
from pygmm import contracts


def _names(cls):
    return {f.name for f in fields(cls)}


def test_response_spectrum_fields():
    assert _names(contracts.ResponseSpectrum) == {
        "periods",
        "spec_accels",
        "damping",
        "duration",
    }


def test_fourier_spectrum_fields():
    assert _names(contracts.FourierSpectrum) == {"freqs", "fourier_amps", "duration"}


def test_duration_fields():
    assert _names(contracts.Duration) == {
        "duration",
        "d_5_75",
        "d_5_95",
        "d_20_80",
        "ln_std",
        "plus_sigma",
        "minus_sigma",
    }


def test_nonlinear_soil_curves_fields():
    assert _names(contracts.NonlinearSoilCurves) == {
        "strains",
        "mod_reduc",
        "damping",
        "damping_min",
        "unit_wt",
        "name",
    }


def test_velocity_profile_fields():
    assert _names(contracts.VelocityProfile) == {
        "depth",
        "vs_median",
        "std_vs_ln",
        "region",
        "site_class",
    }


# ---------------------------------------------------------------------------
# Producer protocols
#
# These replace the seven former `_base.py` ABCs, six of which never acquired
# a subclass. Structural conformance is what pystrata and pyrvt actually rely
# on, so it is what is tested.
# ---------------------------------------------------------------------------

SCENARIO = pygmm.Scenario(
    mag=6.5,
    dist_rup=20.0,
    dist_jb=20.0,
    dist_x=20.0,
    v_s30=760.0,
    depth_tor=5.0,
    dip=90.0,
    mechanism="SS",
    on_hanging_wall=False,
)


def test_response_spectrum_protocol():
    m = pygmm.ChiouYoungs2014(SCENARIO)
    assert isinstance(m, contracts.SupportsResponseSpectrum)
    rs = m.response_spectrum()
    assert isinstance(rs, contracts.ResponseSpectrum)
    assert len(rs.periods) == len(rs.spec_accels)


def test_response_spectrum_vectorized():
    """A vectorized scenario gives one spectrum per scenario."""
    mags = np.array([5.5, 6.5, 7.5])
    m = pygmm.ChiouYoungs2014(pygmm.Scenario(**dict(SCENARIO, mag=mags)))
    rs = m.response_spectrum()
    assert rs.spec_accels.shape == (len(mags), len(rs.periods))
    for i, mag in enumerate(mags):
        scalar = pygmm.ChiouYoungs2014(pygmm.Scenario(**dict(SCENARIO, mag=mag)))
        np.testing.assert_allclose(
            rs.spec_accels[i], scalar.response_spectrum().spec_accels
        )


def test_response_spectrum_with_ims():
    """``ims`` limits the spectrum to the requested periods."""
    m = pygmm.ChiouYoungs2014(SCENARIO, ims=["psa_0p200", "psa_1p000"])
    rs = m.response_spectrum()
    np.testing.assert_allclose(rs.periods, [0.2, 1.0])
    full = pygmm.ChiouYoungs2014(SCENARIO)
    np.testing.assert_allclose(
        rs.spec_accels, full.interp_spec_accels([0.2, 1.0]), rtol=1e-12
    )


@pytest.mark.parametrize(
    "cls",
    [i.cls for i in pygmm.find_models(provides="duration")],
    ids=lambda c: c.__name__,
)
def test_duration_protocol(cls):
    kwds = {"region": "Japan"} if cls.__name__ == "PinillaRamosEtAl2024" else {}
    m = cls(
        pygmm.Scenario(
            **dict(SCENARIO), site_cond="soil", event_type="interface", **kwds
        )
    )
    assert isinstance(m, contracts.SupportsDuration)
    d = m.duration_model()
    assert isinstance(d, contracts.Duration)
    assert np.isfinite(d.duration) and d.duration > 0
    # Exactly one uncertainty convention is populated.
    assert (d.ln_std is None) != (d.plus_sigma is None)


@pytest.mark.parametrize(
    "cls",
    [i.cls for i in pygmm.find_models(provides="duration")],
    ids=lambda c: c.__name__,
)
def test_duration_protocol_vectorized(cls):
    """A vectorized scenario gives one duration per scenario."""
    kwds = {"region": "Japan"} if cls.__name__ == "PinillaRamosEtAl2024" else {}
    mags = np.array([5.5, 6.5, 7.5])
    base = dict(SCENARIO, site_cond="soil", event_type="interface", **kwds)
    d = cls(pygmm.Scenario(**dict(base, mag=mags))).duration_model()
    assert np.shape(d.duration) == (len(mags),)
    for i, mag in enumerate(mags):
        scalar = cls(pygmm.Scenario(**dict(base, mag=mag))).duration_model()
        assert isinstance(scalar.duration, float)
        for name in _names(contracts.Duration):
            value = getattr(d, name)
            if value is None:
                assert getattr(scalar, name) is None
            else:
                np.testing.assert_allclose(value[i], getattr(scalar, name))


@pytest.mark.parametrize(
    "cls",
    [i.cls for i in pygmm.find_models(provides="soil_curves")],
    ids=lambda c: c.__name__,
)
def test_soil_curves_protocol(cls):
    assert issubclass(cls, contracts.SupportsSoilCurves)


def test_fourier_spectrum_protocol_with_self_duration():
    m = pygmm.fourier_spectrum.SourceTheoryModel(
        magnitude=6.5, distance=20.0, region="wna"
    )
    assert isinstance(m, contracts.SupportsFourierSpectrum)
    fs = m.fourier_spectrum()
    assert isinstance(fs, contracts.FourierSpectrum)
    assert fs.duration == pytest.approx(m.duration)


def test_fourier_spectrum_protocol_requires_supplied_duration():
    """BA19 has no intrinsic duration; the signature says so."""
    m = pygmm.BaylessAbrahamson2019(SCENARIO)
    assert isinstance(m, contracts.SupportsFourierSpectrum)
    with pytest.raises(ValueError, match="predicts EAS only"):
        m.fourier_spectrum()
    assert m.fourier_spectrum(duration=8.0).duration == pytest.approx(8.0)


def test_deleted_contracts_are_gone():
    """Four dataclasses had no producer and no consumer; they were removed."""
    for name in (
        "FaultDisplacement",
        "CptSounding",
        "SoilBehaviorProfile",
        "LiquefactionTriggering",
    ):
        assert not hasattr(contracts, name)
