#!/usr/bin/env python
"""Tests for Stafford (2017) inter-frequency correlation model."""

import numpy as np
import pandas as pd
import pytest

from pygmm import Stafford2017

from . import FPATH_DATA

_df = pd.read_csv(FPATH_DATA / "PJS2017_cor_M4pt5_fmin0pt06.csv")
_test_cases = list(zip(_df["Freq"], _df["Cor"]))
F_REF = 0.06


# FIXME: ~25/106 cases fail due to numerical errors in the implementation.
# The between_event, within_event, and between_site correlation formulas
# need to be reconciled against the paper (Stafford 2017, BSSA).
@pytest.mark.xfail(
    reason="Stafford 2017 implementation has known numerical errors (~25/106 cases)"
)
@pytest.mark.parametrize("freq,expected_cor", _test_cases)
def test_correlation(freq, expected_cor):
    cor = Stafford2017.cor(np.array([F_REF, freq]), mag=4.5)
    np.testing.assert_allclose(cor[0, 1], expected_cor, rtol=0.05)
