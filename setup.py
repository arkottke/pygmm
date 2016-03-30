#!/usr/bin/env python
# -*- coding: utf-8 -*-


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [
    'matplotlib',
    'numpy',
    'scipy >= 0.17.0'
]

setup_requirements = [
    'pytest-runner',
]

test_requirements = [
    'pytest >= 2.9.0',
    'pytest-cov',
    'pytest-flake8'
]

setup(
    name='pygmm',
    version='0.3.2',
    description="Ground motion models implemented in Python.",
    long_description=readme + '\n\n' + history,
    author="Albert Kottke",
    author_email='albert.kottke@gmail.com',
    url='https://github.com/arkottke/pygmm',
    packages=[
        'pygmm'
    ],
    package_dir={'pygmm':
                 'pygmm'},
    package_data={'pygmm': ['data/*']},
    include_package_data=True,
    install_requires=requirements,
    license="MIT",
    zip_safe=False,
    keywords='pygmm',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Intended Audience :: Science/Research',
    ],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements
)
