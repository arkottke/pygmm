#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test calculation of CB14 static methods."""
from numpy.testing import assert_almost_equal

from pygmm import CampbellBozorgnia2014 as CB14


def test_depth_2_5():
    # Value calculated from NGAW2 spreadsheet
    assert_almost_equal(CB14.calc_depth_2_5(600, "japan", None), 0.1844427)
    assert_almost_equal(CB14.calc_depth_2_5(600, "california", None), 0.7952589)
