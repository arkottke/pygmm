#!/usr/bin/env python
"""Test model interface using C03 model."""

import warnings

import numpy as np
import pytest

from pygmm.model import CategoricalParameter, NumericParameter, Parameter


def test_parameter():
    p = Parameter("test")
    p.check(None)


@pytest.fixture
def numeric_parameter():
    return NumericParameter("test", required=True, min_=0, max_=10)


@pytest.fixture
def categorical_parameter():
    return CategoricalParameter("test", required=True, options=["spam", "eggs"])


# See https://github.com/pytest-dev/pytest/issues/349
@pytest.fixture(params=["numeric_parameter", "categorical_parameter"])
def param(request):
    return request.getfixturevalue(request.param)


def test_required(param):
    with pytest.raises(ValueError):
        param.check(None)


def test_numeric_scalar_warns(numeric_parameter):
    with pytest.warns(UserWarning, match="less than the recommended limit"):
        assert numeric_parameter.check(-1) == -1
    with pytest.warns(UserWarning, match="greater than the recommended limit"):
        assert numeric_parameter.check(11) == 11


def test_numeric_array_warns_once_per_limit(numeric_parameter):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        value = numeric_parameter.check([-2, -1, 5, 11])
    assert isinstance(value, np.ndarray)
    np.testing.assert_array_equal(value, [-2, -1, 5, 11])
    messages = [str(w.message) for w in caught]
    assert len(messages) == 2
    assert "2 of 4 values, minimum of -2" in messages[0]
    assert "1 of 4 values, maximum of 11" in messages[1]


def test_numeric_array_within_limits_does_not_warn(numeric_parameter):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        numeric_parameter.check(np.array([0.0, 5.0, 10.0]))


def test_categorical_array_replaces_invalid_entries():
    param = CategoricalParameter("test", options=["spam", "eggs"], default="spam")
    with pytest.warns(UserWarning, match="1 of 3 values"):
        value = param.check(np.array(["eggs", "ham", "spam"]))
    np.testing.assert_array_equal(value, ["eggs", "spam", "spam"])


def test_categorical_array_of_options_does_not_warn(categorical_parameter):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        value = categorical_parameter.check(["spam", "eggs"])
    np.testing.assert_array_equal(value, ["spam", "eggs"])
