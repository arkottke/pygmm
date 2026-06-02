"""Cross-package interop tests for pygmm.fourier_spectrum models.

Verifies that the FAS-only ports in pygmm produce numerically identical
spectra (and durations) to the original pyrvt RVT-motion classes they were
extracted from, and that piping them through ``RvtMotion.from_fas`` yields
equivalent PGA / PGV.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

pyrvt = pytest.importorskip("pyrvt")

from pygmm.fourier_spectrum import SourceTheoryModel, StaffordEtAl2022


@pytest.mark.parametrize("region", ["wna", "cena"])
def test_source_theory_matches_pyrvt(region):
    """pygmm SourceTheoryModel reproduces pyrvt.SourceTheoryMotion FAS + duration."""
    kw = dict(magnitude=6.5, distance=20.0, region=region)
    pgm = SourceTheoryModel(**kw)
    prv = pyrvt.motions.SourceTheoryMotion(**kw)
    assert_allclose(pgm.freqs, prv.freqs)
    assert_allclose(pgm.fourier_amps, prv.fourier_amps, rtol=1e-12)
    assert_allclose(pgm.duration, prv.duration, rtol=1e-12)


def test_stafford_2022_matches_pyrvt():
    """pygmm StaffordEtAl2022 reproduces pyrvt.StaffordEtAl22Motion FAS + duration."""
    pgm = StaffordEtAl2022(mag=6.5, dist_rup=20.0)
    prv = pyrvt.motions.StaffordEtAl22Motion(mag=6.5, dist_rup=20.0)
    assert_allclose(pgm.freqs, prv.freqs)
    assert_allclose(pgm.fourier_amps, prv.fourier_amps, rtol=1e-12)
    assert_allclose(pgm.duration, prv.duration, rtol=1e-12)


def test_from_fas_with_source_theory_pgm():
    """Source-theory FAS from pygmm runs through RvtMotion.from_fas and matches pyrvt."""
    pgm = SourceTheoryModel(magnitude=6.0, distance=15.0, region="wna")
    motion = pyrvt.motions.RvtMotion.from_fas(pgm)
    ref = pyrvt.motions.SourceTheoryMotion(magnitude=6.0, distance=15.0, region="wna")
    assert_allclose(motion.calc_pga(), ref.calc_pga(), rtol=1e-10)
    assert_allclose(motion.calc_pgv(), ref.calc_pgv(), rtol=1e-10)


def test_from_fas_with_stafford_pgm():
    """Stafford FAS from pygmm runs through RvtMotion.from_fas (BT15 calculator)."""
    pgm = StaffordEtAl2022(mag=6.0, dist_rup=15.0)
    motion = pyrvt.motions.RvtMotion.from_fas(
        pgm,
        peak_calculator="BT15",
        calc_kwds={"mag": 6.0, "dist": 15.0, "region": "wus"},
    )
    pga = motion.calc_pga()
    pgv = motion.calc_pgv()
    assert np.isfinite(pga) and pga > 0
    assert np.isfinite(pgv) and pgv > 0
