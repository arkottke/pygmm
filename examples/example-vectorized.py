"""Time a vectorized ground motion model against a loop over scenarios.

Computes the Abrahamson, Silva, and Kamai (2014) PGA for 100,000 random
scenarios, first with arrays of scenario values in a single call, and then
one scenario at a time. The loop is timed on a subset and scaled up, since
running all 100,000 scenarios one at a time takes several minutes.
"""

import logging
import time

import numpy as np

import pygmm

# Silence the warnings for scenario values outside of the model limits
logging.disable(logging.WARNING)

N = 100_000
N_LOOP = 2_000

rng = np.random.default_rng(42)
mag = rng.uniform(5.0, 8.0, N)
dist_jb = rng.uniform(0.0, 200.0, N)
depth_tor = rng.uniform(0.0, 10.0, N)
dip = rng.choice([45.0, 60.0, 90.0], N)
width = np.minimum(10 ** (mag - 5.5) * 5.0, 15.0 / np.sin(np.radians(dip)))
dist_rup = np.hypot(dist_jb, depth_tor)
dist_x = np.where(rng.random(N) < 0.5, dist_jb, -dist_jb)
scenario = dict(
    mag=mag,
    dist_rup=dist_rup,
    dist_jb=dist_jb,
    dist_x=dist_x,
    dip=dip,
    width=width,
    depth_tor=depth_tor,
    v_s30=rng.uniform(200.0, 1000.0, N),
    mechanism=rng.choice(["SS", "NS", "RS"], N),
    on_hanging_wall=dist_x >= 0,
    region="california",
)

# Vectorized: all scenarios in one call, computing only PGA
start = time.perf_counter()
model = pygmm.AbrahamsonSilvaKamai2014(pygmm.Scenario(**scenario), ims=["pga"])
pga_vec = model.pga
t_vec = time.perf_counter() - start
print(f"Vectorized: {N:,} scenarios in {t_vec:.2f} s")

# Loop: one scenario at a time
start = time.perf_counter()
pga_loop = np.array(
    [
        pygmm.AbrahamsonSilvaKamai2014(
            pygmm.Scenario(
                **{k: v if np.ndim(v) == 0 else v[i] for k, v in scenario.items()}
            )
        ).pga
        for i in range(N_LOOP)
    ]
)
t_loop = (time.perf_counter() - start) * N / N_LOOP
print(f"Loop:       {N:,} scenarios in about {t_loop:.0f} s (timed on {N_LOOP:,})")
print(f"Speedup:    about {t_loop / t_vec:,.0f} times")

# Both give the same results
np.testing.assert_allclose(pga_vec[:N_LOOP], pga_loop, rtol=1e-12)
