"""Acceleration time series container and readers."""

from __future__ import annotations

import os
import re

import numpy as np
import numpy.typing as npt

__author__ = "Albert Kottke"


class TimeSeries:
    """Acceleration time series.

    Parameters
    ----------
    time_step : float
        Time step [s].
    accels : array_like
        Accelerations [g].
    info : str, optional
        Description of the record.
    """

    def __init__(self, time_step: float, accels: npt.ArrayLike, info: str = ""):
        self._time_step = time_step
        self._accels = np.asarray(accels, dtype=float)
        self._info = info

    @property
    def info(self) -> str:
        return self._info

    @property
    def time_step(self) -> float:
        return self._time_step

    @property
    def accels(self) -> npt.NDArray[np.floating]:
        return self._accels

    @classmethod
    def read_at2(cls, filename: str | os.PathLike) -> TimeSeries:
        """Read a PEER NGA ``.AT2`` file."""
        with open(filename) as fp:
            next(fp)
            info = next(fp).strip()
            next(fp)
            header = next(fp)
            # NGA-West2: "NPTS= 6000, DT= .0100 SEC"; NGA-West1: "6000 0.0100 NPTS, DT"
            match = re.search(r"DT\s*=\s*([0-9.eE+-]+)", header, re.IGNORECASE)
            if match:
                time_step = float(match.group(1))
            else:
                time_step = float(header.replace(",", " ").split()[1])
            accels = np.array([p for line in fp for p in line.split()]).astype(float)
        return cls(time_step, accels, info)
