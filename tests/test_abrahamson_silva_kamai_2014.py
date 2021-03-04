#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test calculation of ASK14 static methods."""
from numpy.testing import assert_almost_equal

from pygmm import AbrahamsonSilvaKamai2014 as ASK14


def test_depth_1_0():
    # Value calculated from NGAW2 spreadsheet
    assert_almost_equal(ASK14.calc_depth_1_0(600), 0.1424470)


def test_depth_tor():
    # Value calculated from NGAW2 spreadsheet
    assert_almost_equal(ASK14.calc_depth_tor(6), 4.2545455)


def test_width():
    # Value calculated from NGAW2 spreadsheet
    assert_almost_equal(ASK14.calc_width(6, 50), 8.9125094)
