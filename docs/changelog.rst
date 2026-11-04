=========
Changelog
=========

This document records all notable changes to pyGMM.

The format is based on `Keep a Changelog <https://keepachangelog.com/en/1.0.0/>`_,
and this project adheres to `Semantic Versioning <https://semver.org/spec/v2.0.0.html>`_.

Unreleased
==========

Added
-----
- Incorporated all of pyGEEMs, which is no longer maintained:

  - ``pygmm.slope_displacement``: Rathje & Saygili (2008), Rathje & Antonakos
    (2011), Watson-Lamprey & Abrahamson (2006), Bray & Travasarou (2007),
    Bray et al. (2017), Cho & Rathje (2020, 2021, 2022), Park & Kishida
    (2019) dam period, Hale & Abrahamson (2019), and rigid-block (Newmark)
    analysis. The optional ``slope`` extra installs numba and pyrotd.
  - ``pygmm.kappa``: Van Houtte et al. (2011) and Xu & Rathje (2021).
  - ``pygmm.ground_motion``: Rezaeian et al. (2014) damping scaling, Rathje
    et al. (2005) periods, and Hayden et al. (2012) pulse proportion.
  - ``pygmm.velocity_profile``: Boore (2016) density and Vp from Vs, and
    Wair, DeJong & Shantz (2012) Vs from SPT blow count.
  - ``pygmm.liquefaction`` (SPT triggering), ``pygmm.elastic`` (elastic
    constants), and ``pygmm.time_series`` (AT2 reader).

- Modern documentation structure with sphinx-design
- Comprehensive examples gallery
- Improved API documentation

Changed
-------
- Migrated build system from Hatch to uv
- Updated documentation theme to Furo
- Modernized development workflow

Version 0.7.4
==============

This section would contain the actual release notes. Since this is documentation
reorganization, the specific version history would be maintained from the existing
HISTORY.md file.

.. include:: ../HISTORY.md
   :parser: myst_parser.sphinx_
