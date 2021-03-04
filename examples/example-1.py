#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot influence of V_s30 predicted by CY14 model."""
import matplotlib.pyplot as plt

import pygmm

fig, ax = plt.subplots()

for v_s30 in [300, 600, 900]:
    s = pygmm.model.Scenario(
        mag=7, dist_jb=20, dist_x=20, dist_rup=25, dip=90, v_s30=v_s30
    )
    m = pygmm.ChiouYoungs2014(s)
    ax.plot(m.periods, m.spec_accels, label=str(v_s30))

ax.set_xlabel("Period (s)")
ax.set_xscale("log")

ax.set_ylabel("5%-Damped Spectral Accel. (g)")
ax.set_yscale("log")

ax.grid()

ax.legend(title="$V_{s30}$ (m/s)")

plt.show()
