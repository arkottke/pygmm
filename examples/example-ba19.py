import os

import matplotlib.pyplot as plt

import pygmm


scenario = pygmm.Scenario(mag=6.5, dist_rup=30, v_s30=500, mechanism="SS", depth_tor=0)

models = [
    pygmm.BaylessAbrahamson2019(scenario),
    pygmm.BaylessAbrahamson2019(
        scenario.copy_with(v_s30=pygmm.BaylessAbrahamson2019.V_REF)
    ),
]

fig, ax = plt.subplots()

for m in models:
    ax.plot(m.freqs, m.eas, label=f"{m.scenario.v_s30:.0f}")

ax.legend(title="$V_{s30}$ m/s")

ax.set(
    xlabel="Frequency (Hz)",
    xscale="log",
    ylabel="Effective Ampl. (g-s)",
    yscale="log",
)

fig.savefig(os.path.splitext(__file__)[0] + ".png", dpi=120)
