#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test model interface using C03 model."""
import pytest

from pygmm.model import CategoricalParameter
from pygmm.model import NumericParameter
from pygmm.model import Parameter


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
