"""Cross-package interop tests for pygmm.fourier_spectrum models.

Verifies that pygmm FAS models produce physically reasonable spectra and that
piping them through ``RvtMotion.from_fas`` yields finite, positive PGA/PGV.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

pyrvt = pytest.importorskip("pyrvt")

from pygmm.fourier_spectrum import SourceTheoryModel, StaffordEtAl2022


@pytest.mark.parametrize("region", ["wna", "cena"])
def test_source_theory_physical_ranges(region):
    """SourceTheoryModel produces finite, positive FAS and duration."""
    m = SourceTheoryModel(magnitude=6.5, distance=20.0, region=region)
    assert np.all(np.isfinite(m.freqs)) and np.all(m.freqs > 0)
    assert np.all(np.isfinite(m.fourier_amps)) and np.all(m.fourier_amps > 0)
    assert np.isfinite(m.duration) and m.duration > 0


def test_stafford_2022_physical_ranges():
    """StaffordEtAl2022 produces finite, positive FAS and duration."""
    m = StaffordEtAl2022(mag=6.5, dist_rup=20.0)
    assert np.all(np.isfinite(m.freqs)) and np.all(m.freqs > 0)
    assert np.all(np.isfinite(m.fourier_amps)) and np.all(m.fourier_amps > 0)
    assert np.isfinite(m.duration) and m.duration > 0


def test_from_fas_with_source_theory_pgm():
    """pygmm SourceTheoryModel → RvtMotion.from_fas gives finite, positive PGA/PGV."""
    pgm = SourceTheoryModel(magnitude=6.0, distance=15.0, region="wna")
    motion = pyrvt.motions.RvtMotion.from_fas(pgm)
    pga = motion.calc_pga()
    pgv = motion.calc_pgv()
    # Loose physical bounds for WNA M6 @ 15 km
    assert 0.001 < pga < 10.0
    assert 0.1 < pgv < 200.0


def test_from_fas_stub():
    """RvtMotion.from_fas duck-types: SimpleNamespace stub works, no isinstance."""
    from types import SimpleNamespace

    stub = SimpleNamespace(
        freqs=np.logspace(-1, 2, 50),
        fourier_amps=np.ones(50) * 0.01,
        duration=5.0,
    )
    motion = pyrvt.motions.RvtMotion.from_fas(stub)
    assert np.all(np.isfinite(motion.fourier_amps))
    assert motion.duration == pytest.approx(5.0)


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
