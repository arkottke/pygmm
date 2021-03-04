#!/usr/bin/env python3
# encoding: utf-8
"""Hermkes, Kuehn, Riggelsen (2014, :cite:`hermkes14`) model."""
import logging
import pathlib

import numpy as np
from scipy.interpolate import NearestNDInterpolator

from . import model

__author__ = "Albert Kottke"

fname_data = pathlib.Path(__file__).parent.joinpath(
    "data", "hermkes_kuehn_riggelsen_2014.npz"
)

if not fname_data.exists():
    # Download the model data if not found.
    import urllib.request
    import shutil

    url = (
        "https://www.dropbox.com/s/1tu9ss1s3inctej/"
        "hermkes_kuehn_riggelsen_2014.npz?dl=0"
    )

    try:
        urllib.request.urlretrieve(url, str(fname_data))
    except urllib.request.URLError:
        logging.critical(
            "Hermkes, Kuehn, and Riggelsen (2013) model data required, "
            "which cannot be downloaded. Download the file from %s"
            "to this location: %s",
            url,
            fname_data,
        )

INTERPOLATOR = None


class HermkesKuehnRiggelsen2014(model.GroundMotionModel):
    """Hermkes, Kuehn, Riggelsen (2014, :cite:`hermkes14`) model.

    Only the *GPSELinCorr* model is implemented. This model must be imported
    directly by::

        from pygmm.hermkes_kuehn_riggelsen_2014 import
        HermkesKuehnRiggelsen2014

    This is to due to the large file size of the model data, which takes
    time to load.

    Note that this model was developed using a Bayesian non-parametric
    method, which means it is should only be used over the data range
    used to develop the model. See the paper for more details.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = "Hermkes, Kuehn, Riggelsen (2014)"
    ABBREV = "HKR14"

    # Reference velocity (m/sec)
    V_REF = None

    PERIODS = np.array([-1, 0.01, 0.1, 0.5, 1.0, 4.0])
    INDICES_PSA = np.arange(1, 6)
    INDEX_PGA = 1
    INDEX_PGV = 0
    PARAMS = [
        model.NumericParameter("depth_hyp", False, 0, 40, 15),
        model.NumericParameter("dist_jb", False, 0, 200),
        model.NumericParameter("mag", True, 4, 8),
        model.NumericParameter("v_s30", True, 100, 1200),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
    ]

    def __init__(self, scenario: model.Scenario):
        """Initialize the model."""
        super().__init__(scenario)

        s = self._scenario

        flag_rs = flag_ss = flag_ns = 0
        if s.mechanism == "SS":
            flag_ss = 1
        elif s.mechanism == "NS":
            flag_ns = 1
        elif s.mechanism == "RS":
            flag_rs = 1

        event = (s.mag, s.depth_hyp, flag_rs, flag_ss, flag_ns, s.dist_jb, s.v_s30)

        global INTERPOLATOR
        if INTERPOLATOR is None:
            with np.load(fname_data) as data:
                INTERPOLATOR = NearestNDInterpolator(
                    data["events"], data["predictions"]
                )
        prediction = INTERPOLATOR(event)
        self._ln_resp = prediction[0::2]
        self._ln_std = np.sqrt(prediction[1::2])
