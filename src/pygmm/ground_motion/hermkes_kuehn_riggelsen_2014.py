#!/usr/bin/env python3
"""Hermkes, Kuehn, Riggelsen (2014, :cite:p:`hermkes14`) model."""

import logging
import pathlib
import zipfile

import numpy as np
from scipy.interpolate import NearestNDInterpolator

from .. import model
from ..registry import register

__author__ = "Albert Kottke"

fname_data = pathlib.Path(__file__).parent.parent.joinpath(
    "data", "hermkes_kuehn_riggelsen_2014.npz"
)

if not zipfile.is_zipfile(fname_data):
    # Download the model data if it is not found. The model data (164 MB) is not
    # included in the package, and a clone without Git LFS has a small pointer
    # file instead of the data.
    import urllib.request

    # dl=1 downloads the file; dl=0 returns the Dropbox web page
    url = (
        "https://www.dropbox.com/s/1tu9ss1s3inctej/"
        "hermkes_kuehn_riggelsen_2014.npz?dl=1"
    )

    # Download to a temporary file, so that a failed download does not leave a
    # file that is not the model data
    fname_tmp = fname_data.with_suffix(".npz.download")
    try:
        urllib.request.urlretrieve(url, str(fname_tmp))
        if not zipfile.is_zipfile(fname_tmp):
            raise ValueError("the downloaded file is not the model data")
        fname_tmp.replace(fname_data)
    except (OSError, ValueError) as error:
        fname_tmp.unlink(missing_ok=True)
        logging.critical(
            "Hermkes, Kuehn, and Riggelsen (2013) model data required, "
            "which cannot be downloaded (%s). Download the file from %s "
            "to this location: %s",
            error,
            url,
            fname_data,
        )

INTERPOLATOR = None


@register(lazy_data=True)
class HermkesKuehnRiggelsen2014(model.GroundMotionModel):
    """Hermkes, Kuehn, Riggelsen (2014, :cite:p:`hermkes14`) model.

    Only the *GPSELinCorr* model is implemented. This model must be imported
    directly by::

        from pygmm.ground_motion.hermkes_kuehn_riggelsen_2014 import
        HermkesKuehnRiggelsen2014

    This is to due to the large file size of the model data, which takes
    time to load.

    Note that this model was developed using a Bayesian non-parametric
    method, which means it is should only be used over the data range
    used to develop the model. See the paper for more details.

    The model is vectorized. Each scenario value (``mag``, ``depth_hyp``,
    ``dist_jb``, ``v_s30``, and ``mechanism``) can be a scalar or an array,
    and the arrays are broadcast against each other. For a scalar scenario,
    the response and standard deviation have one value per period, as in
    other models. For arrays of N scenarios, they have shape (N, periods), so,
    for example, ``pga`` has shape (N,).

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario
    ims : str or sequence of str, optional
        intensity measures to compute: "pga", "pgv", spectral periods such as
        "psa_1p000" (1.0 s), and/or "psa_all" (all 5 periods). The nearest
        neighbor lookup provides all of the intensity measures at once, but
        only the requested intensity measures are kept. If *None* (default),
        all intensity measures are computed.

    Examples
    --------
    The examples are skipped in the doctests because of the large model data.

    >>> import numpy as np
    >>> import pygmm
    >>> from pygmm.ground_motion.hermkes_kuehn_riggelsen_2014 import (
    ...     HermkesKuehnRiggelsen2014,
    ... )
    >>> s = pygmm.Scenario(
    ...     mag=np.array([6.0, 7.0]), dist_jb=np.array([10.0, 50.0]),
    ...     v_s30=400.0, mechanism=np.array(["SS", "RS"]))
    >>> HermkesKuehnRiggelsen2014(s, ims=["pga"]).pga.shape  # doctest: +SKIP
    (2,)
    >>> HermkesKuehnRiggelsen2014(s).spec_accels.shape  # doctest: +SKIP
    (2, 5)

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

    def __init__(self, scenario: model.Scenario, ims=None):
        """Initialize the model.

        Parameters
        ----------
        scenario : pygmm.model.Scenario
            earthquake scenario.
        ims : str or sequence of str, optional
            intensity measures to compute. If *None* (default), all intensity measures
            are computed.
        """
        super().__init__(scenario, ims)

        s = self._scenario

        # Mechanism flags. Other values (e.g., an invalid mechanism replaced
        # by the default of *None*) have all flags equal to zero.
        flag_rs = model.equals(s.mechanism, "RS").astype(int)
        flag_ss = model.equals(s.mechanism, "SS").astype(int)
        flag_ns = model.equals(s.mechanism, "NS").astype(int)

        # The interpolator accepts a tuple of arrays, which are broadcast
        # against each other, and gives the predictions along the last axis
        event = (s.mag, s.depth_hyp, flag_rs, flag_ss, flag_ns, s.dist_jb, s.v_s30)

        global INTERPOLATOR
        if INTERPOLATOR is None:
            with np.load(fname_data) as data:
                INTERPOLATOR = NearestNDInterpolator(
                    data["events"], data["predictions"]
                )
        prediction = INTERPOLATOR(event)
        # The predictions are pairs of the mean and variance for each period
        ln_resp = prediction[..., 0::2]
        ln_std = np.sqrt(prediction[..., 1::2])
        # The lookup provides all periods, so select the requested ones
        if self._indices is not None:
            ln_resp = model.take_periods(ln_resp, self._indices)
            ln_std = model.take_periods(ln_std, self._indices)
        self._ln_resp = ln_resp
        self._ln_std = ln_std
