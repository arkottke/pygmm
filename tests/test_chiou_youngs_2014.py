#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test calculation of CY14 static methods."""
from numpy.testing import assert_almost_equal

from pygmm import ChiouYoungs2014 as CY14


def test_depth_1_0():
    # Value calculated from NGAW2 spreadsheet
    assert_almost_equal(CY14.calc_depth_1_0(600, "california"), 0.1259, 4)
    assert_almost_equal(CY14.calc_depth_1_0(600, "japan"), 0.0331, 4)


def test_depth_tor():
    # Value calculated from NGAW2 spreadsheet
    assert_almost_equal(CY14.calc_depth_tor(6, "SS"), 2.2587685)
    assert_almost_equal(CY14.calc_depth_tor(6, "RS"), 6.3447262)
