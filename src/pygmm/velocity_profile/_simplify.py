"""Simplification of layered shear-wave velocity profiles."""

from __future__ import annotations

import dataclasses

import numpy as np

from ..contracts import VelocityProfile


def simplify_velocity_profile[T: VelocityProfile](vp: T, tol: float = 0.02) -> T:
    """Merge adjacent layers with similar shear-wave velocity.

    Consecutive layers whose ``vs_median`` is within *tol* (relative) of the
    first layer of the group are merged. The merged layer takes the
    travel-time-averaged velocity, so the total travel time through the
    merged layers is preserved. Other per-layer quantities are averaged with
    thickness weights: ``std_vs_ln`` and, when the profile has one,
    ``density``. The last point is the half-space and has no thickness.

    Parameters
    ----------
    vp : VelocityProfile
        Profile whose ``depth`` holds layer tops [m]. Subclasses are returned
        as the same type; attributes other than the ones above are kept.
    tol : float, optional
        Relative velocity tolerance for merging (default 0.02, i.e. 2 %).

    Returns
    -------
    VelocityProfile
        Simplified profile of the same type as *vp*.
    """
    depth = np.asarray(vp.depth, dtype=float)
    vs = np.asarray(vp.vs_median, dtype=float)
    std = np.asarray(vp.std_vs_ln, dtype=float)
    density = getattr(vp, "density", None)
    if density is not None:
        density = np.asarray(density, dtype=float)

    # The last layer is the half-space and gets zero thickness
    thickness = np.diff(depth, append=depth[-1])

    groups: list[list[int]] = []
    current = [0]
    for i in range(1, len(vs)):
        if abs(vs[i] - vs[current[0]]) / vs[current[0]] <= tol:
            current.append(i)
        else:
            groups.append(current)
            current = [i]
    groups.append(current)

    new_depth, new_vs, new_std, new_density = [], [], [], []
    for group in groups:
        thick = thickness[group]
        total = thick.sum()
        # Weights fall back to uniform for the zero-thickness half-space
        weights = thick if total > 0 else np.ones(len(group))

        new_depth.append(depth[group[0]])
        # Time-averaged velocity: total thickness / sum(thickness / Vs)
        new_vs.append(total / np.sum(thick / vs[group]) if total > 0 else vs[group[0]])
        new_std.append(np.average(std[group], weights=weights))
        if density is not None:
            new_density.append(np.average(density[group], weights=weights))

    changes: dict[str, np.ndarray] = {
        "depth": np.array(new_depth),
        "vs_median": np.array(new_vs),
        "std_vs_ln": np.array(new_std),
    }
    if density is not None:
        changes["density"] = np.array(new_density)
    return dataclasses.replace(vp, **changes)
