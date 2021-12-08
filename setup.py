#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import find_packages
from setuptools import setup

with open("README.rst") as fp:
    readme = fp.read()

with open("HISTORY.rst") as fp:
    history = fp.read()

setup(
    name="pygmm",
    version="0.6.3",
    description="Ground motion models implemented in Python.",
    long_description=readme + "\n\n" + history,
    author="Albert Kottke",
    author_email="albert.kottke@gmail.com",
    url="https://github.com/arkottke/pygmm",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "matplotlib",
        "numpy",
        "scipy >= 0.17.0",
        "six",
    ],
    license="MIT",
    zip_safe=False,
    keywords="pygmm",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering",
        "Intended Audience :: Science/Research",
    ],
    setup_requires=[
        "pytest-runner",
    ],
    tests_require=[
        "pytest >= 2.9.0",
        "pytest-cov",
        "pytest-flake8",
    ],
)
