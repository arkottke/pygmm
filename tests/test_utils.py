"""Tests for the pygmm._utils helpers."""

import numpy as np
import pytest

from pygmm._utils import check_bounds, check_options


def test_check_bounds_array_raises_value_error():
    with pytest.raises(ValueError, match="below lower limit"):
        check_bounds(np.array([1.0, -1.0]), 0.0, 2.0, "x")
    with pytest.raises(ValueError, match="above upper limit"):
        check_bounds(np.array([1.0, 3.0]), 0.0, 2.0, "x")
    check_bounds(np.array([0.5, 1.5]), 0.0, 2.0, "x")


@pytest.mark.parametrize("value", ["", "bc", "bcd"])
def test_check_options_string_needs_whole_option(value):
    with pytest.raises(ValueError, match="not in options"):
        check_options(value, "bcd", "site_class")


def test_check_options_accepts_member():
    check_options("c", "bcd", "site_class")
