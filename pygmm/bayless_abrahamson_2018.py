# -*- coding: utf-8 -*-
"""Bayless and Abrahamson (2018, :cite:`bayless19`) correlation."""
from typing import Optional

import numpy as np
import numpy.typing as npt

from . import model

__author__ = "Albert Kottke"


class BaylessAbrahamson2018:
    """Bayless and Abrahamson (2018) model.

    This inter-frequency correlation model was developed for active tectonic regions.
    """

    #: Long name of the model
    NAME = "Bayless and Abrahamson (2018)"
    #: Short name of the model
    ABBREV = "BA18"

    # Load the coefficients for the model
    COEFF = model.load_data_file("bayless_abrahamson_2018.csv", 3)
    FREQS = COEFF["freq_hz"]

    @classmethod
    def corr(cls, freqs: npt.ArrayLike) -> np.ndarray:
        """
        Parameters
        ----------
        freq : array_like
            Frequencies

        Returns
        -------
        corr : np.ndarray
            Correlation coefficient matrix
        """

        # Create two matrices of frequencies
        n = freqs.shape[0]
        triu = np.triu(np.full((n, n), freqs))
        freqs_mat = triu + triu.T - np.diag(freqs)

        def do_interp(f_m, attr):
            fp = getattr(cls.COEFF, attr)
            return np.interp(f_m, cls.FREQS, fp, left=fp[0], right=fp[-1])

        # Compute the frequency parameters
        f_r = np.abs(np.log(freqs / freqs_mat))
        f_m = np.minimum(freqs, freqs_mat)
        a, b, c, d = [do_interp(f_m, param) for param in "ABCD"]
        # Lower matrix
        lower = np.tanh(a * np.exp(b * f_r) + c * np.exp(d * f_r)).reshape(n, n)

        corr = lower + lower.T
        np.fill_diagonal(corr, 1)

        return corr

    @classmethod
    def cov(
        cls,
        freqs,
        *,
        std: Optional[npt.ArrayLike] = None,
        component: Optional[str] = None
    ) -> np.ndarray:
        """
        Parameters
        ----------
        freq : array_like
            Frequencies [Hz]
        std : array_like, optional
            Standard deviation at the frequencies
        component : str, optional
            Component of the published standard deviation. Possible options: tau, phi_s2s, phi_s, sigma
        Returns
        -------
        corr : np.ndarray
            Correlation coefficient matrix
        """

        corr = cls.corr(freqs)
        if std is None and component is not None:
            fp = getattr(cls.COEFF, component)
            ln_std = np.interp(freqs, cls.FREQS, fp, left=fp[0], right=fp[-1])
        else:
            raise NotImplementedError

        n = len(freqs)
        ln_std = np.full((n, n), ln_std)
        return corr * ln_std * ln_std.T - np.diag(ln_std)
