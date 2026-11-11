"""USGS National Crustal Model (NCM) geophysical profiles."""

from __future__ import annotations

import gzip
import json
import os
import urllib.request
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt

from ..contracts import VelocityProfile

NCM_BASE_URL = "https://earthquake.usgs.gov/ws/nshmp/ncm/geophysical"


@dataclass(frozen=True, kw_only=True)
class NcmProfile(VelocityProfile):
    """Shear-wave velocity profile from the National Crustal Model.

    Adds the quantities the NCM reports at the site to
    :class:`~pygmm.contracts.VelocityProfile`, so it can be used anywhere a
    ``VelocityProfile`` is accepted.

    Attributes
    ----------
    density : np.ndarray
        Mass density at each depth [kg/m³].
    water_table_depth : float
        Depth to the water table [m].
    """

    density: npt.NDArray[np.floating]
    water_table_depth: float


def fetch_ncm_profile(
    latitude: float,
    longitude: float,
    fpath: str | os.PathLike,
    depth_start: float = 0,
    depth_step: float = 5,
    depth_end: float = 1e4,
    gzip_output: bool = True,
    timeout: float = 60,
) -> None:
    """Fetch a geophysical profile from the USGS National Crustal Model.

    Queries the NCM geophysical endpoint for a single location and saves the
    JSON response to *fpath*, which :func:`load_ncm_profile` reads.

    Parameters
    ----------
    latitude : float
        Latitude of the site [degrees].
    longitude : float
        Longitude of the site [degrees].
    fpath : str or path-like
        Output file path.
    depth_start : float, optional
        Minimum depth [m] (default 0).
    depth_step : float, optional
        Depth increment [m] (default 5).
    depth_end : float, optional
        Maximum depth [m] (default 10 000).
    gzip_output : bool, optional
        If True (default), write a gzip-compressed JSON file.
    timeout : float, optional
        Request timeout [s] (default 60).

    Raises
    ------
    urllib.error.HTTPError
        If the request to the NCM service fails.
    RuntimeError
        If the response status is not successful.
    """
    url = (
        f"{NCM_BASE_URL}"
        f"?location={latitude},{longitude}"
        f"&depths={depth_start},{depth_step},{depth_end}"
    )
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        data = json.loads(resp.read().decode())

    # The service has returned the misspelled "sucess"; accept both
    if data.get("status") not in ("sucess", "success"):
        raise RuntimeError(f"NCM request failed with status: {data.get('status')}")

    opener = gzip.open if gzip_output else open
    with opener(os.fspath(fpath), "wt") as fp:
        json.dump(data, fp)


def load_ncm_profile(fpath: str | os.PathLike) -> NcmProfile:
    """Load a National Crustal Model profile saved by :func:`fetch_ncm_profile`.

    Parameters
    ----------
    fpath : str or path-like
        Path to the NCM JSON file, plain or gzip-compressed.

    Returns
    -------
    NcmProfile
        Profile with ``depth`` [m], ``vs_median`` [m/s], ``density`` [kg/m³]
        and ``water_table_depth`` [m]. ``std_vs_ln`` is zero because the NCM
        gives no uncertainty.

    Raises
    ------
    ValueError
        If the response contains no results.
    """
    fpath = os.fspath(fpath)

    # Try gzip first, fall back to plain JSON
    try:
        with gzip.open(fpath, "rt") as fp:
            data = json.load(fp)
    except gzip.BadGzipFile:
        with open(fpath) as fp:
            data = json.load(fp)

    results = data["response"]["results"]
    if not results:
        raise ValueError("NCM response contains no results")
    result = results[0]

    depth = np.array(data["request"]["depths"]["depth_vector"], dtype=float)
    vs = np.array(result["profile"]["vs"], dtype=float)
    density = np.array(result["profile"]["density"], dtype=float)

    return NcmProfile(
        depth=depth,
        vs_median=vs,
        std_vs_ln=np.zeros_like(vs),
        density=density,
        water_table_depth=float(result["location"]["water_table_depth"]),
    )
