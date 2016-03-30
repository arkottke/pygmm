#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function

import os

from invoke import run, task
from invoke.util import log

use_pty = True if os.name == 'posix' else False


@task
def clean():
    """clean - remove build artifacts."""
    run('rm -rf build/')
    run('rm -rf dist/')
    run('rm -rf pygmm.egg-info')
    run('find . -name __pycache__ -delete')
    run('find . -name *.pyc -delete')
    run('find . -name *.pyo -delete')
    run('find . -name *~ -delete')

    log.info('cleaned up')


@task
def test():
    """test - run the test runner."""
    run('py.test --flake8 --cov-report html --cov tests/', pty=use_pty)
    run('open htmlcov/index.html')


@task
def lint():
    """lint - check style with flake8."""
    run('flake8 pygmm tests')


@task(clean)
def publish():
    """publish - package and upload a release to the cheeseshop."""
    run('python setup.py sdist upload', pty=use_pty)
    run('python setup.py bdist_wheel upload', pty=use_pty)

    log.info('published new release')
