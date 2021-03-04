#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create plots for visual comparison."""
import os

import matplotlib

matplotlib.use("agg")  # NOQA
import matplotlib.pyplot as plt
import pytest

import pygmm

DEFAULT_PROPS = dict(
    depth_2_5=5,
    depth_bor=15,
    depth_hyp=9,
    depth_tor=5,
    dip=90.0,
    dist=20,
    dist_jb=30.0,
    dist_rup=30.0,
    dist_x=30.0,
    dpp_centered=0,
    flag_hw=0,
    flag_meas=0,
    mag=6,
    mechanism="SS",
    region="california",
    v_s30=500.0,
    width=10,
)

# Make the figure directory if needed
if not os.path.exists("figures"):
    os.makedirs("figures")


@pytest.mark.parametrize("model", pygmm.models, ids=lambda m: m.ABBREV)
@pytest.mark.parametrize(
    "key,values,label",
    [
        ("mag", [5, 6, 7], "Magnitude"),
        (["dist", "dist_rup", "dist_jb", "dist_x"], [10, 50, 100], "Distance (km)"),
        ("v_s30", [300, 650, 1000], "$V_{s30}$ (m/s)"),
    ],
    ids=lambda a: a[0],
)
def plot_model_with_param(model, key, values, label):
    props = dict(DEFAULT_PROPS)

    fig, ax = plt.subplots()
    for v in values:
        if isinstance(key, str):
            props[key] = v
        else:
            for k in key:
                props[k] = v
        m = model(**props)
        ax.plot(m.periods, m.spec_accels, label="%g" % v)

    ax.set_xlabel("Period (s)")
    try:
        ax.set_xscale("log")
    except ValueError:
        pass

    ax.set_ylabel("5% Damped, Spectral. Accel. (g)")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 1e1)

    ax.legend(loc="upper right", title=label, fontsize="x-small")

    ax.grid()

    fig.tight_layout()

    if isinstance(key, str):
        prefix = key
    else:
        prefix = key[0]

    fig.savefig(os.path.join("figures", prefix + "-" + model.ABBREV))
    plt.close(fig)
