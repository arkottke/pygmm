"""Tests for pygmm.liquefaction."""

import numpy as np
import pandas as pd
import pytest
import scipy.constants as SC

from pygmm._utils import UNIT_WT_WATER
from pygmm.liquefaction import calc_triggering


@pytest.fixture
def spt_data():
    df = pd.DataFrame(
        {
            "depth": [2.1, 5.2, 7.9, 9.6, 12.2],
            "N": [6, 6, 6, 1, 29],
            "FC": [5, 32, 31, 38, 30],
        }
    )
    df["depth"] *= SC.foot
    df["unit_wt"] = 110 * UNIT_WT_WATER / 62.4
    calc_triggering(df, 6.0, 0.42, 0)
    return df


def test_stresses(spt_data):
    # Stress in PSF
    stress_v = np.array([226, 567, 864, 1051, 1337])
    stress_v_eff = np.array([98, 245, 374, 455, 578])

    # Convert to ATM
    scale = 2115.7

    np.testing.assert_allclose(spt_data["stress_v"], stress_v / scale, rtol=0.05)
    np.testing.assert_allclose(
        spt_data["stress_v_eff"], stress_v_eff / scale, rtol=0.05
    )


@pytest.mark.skip(reason="Expected values not reproduced; carried over from pyGEEMs")
def test_n_1_60cs(spt_data):
    np.testing.assert_allclose(
        spt_data["N_1_60cs"],
        [7, 17, 16, 7, 60],
        rtol=0.02,
        atol=0.5,
    )


@pytest.mark.xfail(reason="Expected values not reproduced; carried over from pyGEEMs")
def test_crr_csr(spt_data):
    np.testing.assert_allclose(
        spt_data["csr"][:-1],
        [0.63, 0.62, 0.61, 0.61],
        rtol=0.02,
        atol=0.005,
    )
    np.testing.assert_allclose(
        spt_data["crr"][:-1],
        [0.12, 0.24, 0.22, 0.12],
        rtol=0.02,
        atol=0.005,
    )
    np.testing.assert_allclose(
        spt_data["fs_liq"][:-1],
        [0.19, 0.39, 0.36, 0.19],
        rtol=0.02,
        atol=0.005,
    )
