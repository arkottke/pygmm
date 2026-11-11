"""Al Atik and Abrahamson (2021) reference shear-wave velocity profiles."""

from __future__ import annotations

import gzip
import json
from functools import cache
from pathlib import Path

import numpy as np
import numpy.typing as npt

from ..contracts import VelocityProfile

_DATA_FILE = Path(__file__).parent.parent / "data" / "aaa21-profiles.json.gz"


@cache
def _load_profiles() -> dict:
    with gzip.open(_DATA_FILE, "rt") as fp:
        return json.load(fp)


def aaa21_profile(
    model: str,
    vs30: int,
    depth: npt.ArrayLike | None = None,
) -> VelocityProfile:
    """Reference shear-wave velocity profile from Al Atik and Abrahamson (2021).

    The profiles are 1D reference Vs profiles that are compatible with the
    site terms of the NGA-West2 ground-motion models, tabulated to 6 km.

    Parameters
    ----------
    model : str
        Ground-motion model the profile is compatible with. One of
        ``'ASK14'``, ``'BSSA14'``, ``'CB14'``, or ``'CY14'`` (case
        insensitive).
    vs30 : int
        Reference time-averaged velocity in the upper 30 m [m/s]. One of
        ``620``, ``760``, or ``1100``.
    depth : array_like, optional
        Depths to layer tops [m]. If omitted, the tabulated layers are
        returned, with a layer at the surface (0 m) added. Otherwise the
        velocity of the tabulated layer that contains each depth is returned.

    Returns
    -------
    VelocityProfile
        Profile with ``depth``, ``vs_median`` and a zero ``std_vs_ln``.

    Raises
    ------
    KeyError
        If *model* or *vs30* is not tabulated.

    References
    ----------
    Al Atik, L., and N. Abrahamson (2021). A methodology for the development
    of 1D reference VS profiles compatible with ground-motion prediction
    equations: Application to NGA-West2 GMPEs. *BSSA*, 111(4), 1765-1783.
    """
    data = _load_profiles()

    model_key = model.upper()
    if model_key not in data:
        raise KeyError(f"Unknown model '{model}'. Choose from: {list(data)}")

    vs30_key = f"vs30_{vs30}_mps"
    if vs30_key not in data[model_key]:
        raise KeyError(
            f"Unknown vs30={vs30}. Choose from: "
            f"{[int(k.split('_')[1]) for k in data[model_key]]}"
        )

    table = data[model_key][vs30_key]
    depth_tab = np.asarray(table["depth_km"], dtype=float) * 1000
    vs_tab = np.asarray(table["vs_km_per_sec"], dtype=float) * 1000

    # The table starts a little below the surface; extend the first layer up
    if depth_tab[0] > 0:
        depth_tab = np.insert(depth_tab, 0, 0.0)
        vs_tab = np.insert(vs_tab, 0, vs_tab[0])

    if depth is None:
        depth_out = depth_tab
        vs_out = vs_tab
    else:
        depth_out = np.asarray(depth, dtype=float)
        if np.any(depth_out < 0):
            raise ValueError("depth must be non-negative")
        index = np.searchsorted(depth_tab, depth_out, side="right") - 1
        vs_out = vs_tab[index]

    return VelocityProfile(
        depth=depth_out,
        vs_median=vs_out,
        std_vs_ln=np.zeros_like(vs_out),
        site_class=f"AA21 {model_key} {vs30}",
    )
