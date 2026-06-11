"""Stafford (2017, :cite:`safford2017`) correlation."""

import json
import os

import numpy as np
import pandas as pd

__author__ = "Mahdi Bahrampouri"


class Stafford2017:
    """Stafford (2017) model.

    This inter-frequency correlation model was developed for active tectonic regions.
    """

    #: Long name of the model
    NAME = "Stafford (2017)"
    #: Short name of the model
    ABBREV = "PJS17"

    # Load the coefficients for the model
    COEFF = pd.Series(
        json.load(
            open(os.path.join(os.path.dirname(__file__), "data", "stafford_2017.json"))
        )
    )

    @staticmethod
    def _sigmoid(f, alpha0, alpha1, alpha2):
        """
        Compute sigmoid function used throughout the model.

        Parameters:
        -----------
        f : float or array-like
            Frequency value(s)
        alpha0, alpha1, alpha2 : float
            Sigmoid function parameters

        Returns:
        --------
        float or array-like
            Sigmoid function value(s)
        """
        return alpha0 / (1 + np.exp(-alpha2 * np.log(f / alpha1)))

    @classmethod
    def compute_corner_frequency(cls, mag):
        """
        Compute the corner frequency based on magnitude.

        This implementation follows the model in the paper using the Yenier and Atkinson
        (2015) stress parameter model with a single-corner source spectrum.
        """
        # Approximate implementation based on common scaling relationships
        # A more accurate implementation would use the full YA15 model parameters
        return 10 ** (cls.COEFF.CF_C1 - cls.COEFF.CF_C2 * mag)

    @classmethod
    def compute_gamma_E(cls, f_min, mag):
        """
        Compute gamma_E parameter for between-event correlations.
        Uses equation 12 from the paper.
        """
        # Corner frequency for given magnitude
        fc = cls.compute_corner_frequency(mag)

        # Normalized frequency
        f_min_norm = f_min / fc

        # Apply equation 12
        sigmoid_term = cls._sigmoid(
            f_min_norm, cls.COEFF.gamma_E1, cls.COEFF.gamma_E2, cls.COEFF.gamma_E3
        )
        log_term = cls.COEFF.gamma_E4 * np.log(f_min_norm / cls.COEFF.gamma_E2)

        return cls.COEFF.gamma_E0 + sigmoid_term + log_term

    @classmethod
    def compute_eta_A(cls, f_min):
        """
        Compute eta_A parameter for within-event correlations nugget effect.
        Uses equation 16 from the paper.
        """
        # Apply equation 16
        term1 = cls._sigmoid(
            f_min, cls.COEFF.eta_A0, cls.COEFF.eta_A1, cls.COEFF.eta_A2
        )
        term2 = cls._sigmoid(
            f_min, cls.COEFF.eta_A3, cls.COEFF.eta_A4, cls.COEFF.eta_A5
        )

        return term1 * (1 + term2)

    @classmethod
    def compute_gamma_A(cls, f_min):
        """
        Compute gamma_A parameter for within-event correlations.
        Uses equation 17 from the paper.
        """
        # Apply equation 17
        f_ref = cls.COEFF.gamma_A_FREF
        min_f = np.minimum(f_min, f_ref)
        S_term = min_f / (min_f + cls.COEFF.gamma_A2 * (1 - min_f / f_ref))
        log_term = cls.COEFF.gamma_A3 * np.log(np.maximum(f_min, f_ref) / f_ref) ** 2

        return cls.COEFF.gamma_A0 + S_term * cls.COEFF.gamma_A1 + log_term

    @classmethod
    def compute_eta_S(cls, f_min):
        """
        Compute eta_S parameter for between-site correlations nugget effect.
        Uses equation 18 from the paper.
        """
        f_ref1 = cls.COEFF.eta_S_FREF1  # 0.25
        f_ref2 = cls.COEFF.eta_S_FREF2  # 4.0

        # Apply equation 18
        term1 = cls.COEFF.eta_S0 * np.log(
            np.maximum(np.minimum(f_min, f_ref2), f_ref1) / f_ref1
        )
        term2 = cls.COEFF.eta_S1 * np.log(np.maximum(f_min, f_ref2) / f_ref2)

        return term1 + term2

    @classmethod
    def compute_gamma_S(cls, f_min):
        """
        Compute gamma_S parameter for between-site correlations.
        Uses equation 19 from the paper.
        """
        # Apply equation 19
        sigmoid1 = cls._sigmoid(
            f_min, cls.COEFF.gamma_S1, cls.COEFF.gamma_S2, cls.COEFF.gamma_S3
        )
        sigmoid2 = cls._sigmoid(
            f_min, cls.COEFF.gamma_S4, cls.COEFF.gamma_S5, cls.COEFF.gamma_S6
        )

        return cls.COEFF.gamma_S0 + sigmoid1 + sigmoid2

    @classmethod
    def between_event_correlation(cls, f_i, f_j, mag):
        """
        Compute the between-event correlation between two frequencies.
        Uses equation 11 from the paper.
        """
        f_min = np.minimum(f_i, f_j)
        f_max = np.maximum(f_i, f_j)

        # Corner frequency for given magnitude
        fc = cls.compute_corner_frequency(mag)

        # Normalized frequencies
        f_min_norm = f_min / fc
        f_max_norm = f_max / fc

        # Compute gamma_E
        gamma_E = cls.compute_gamma_E(f_min_norm, mag)

        # Apply equation 11
        correlation = np.exp(gamma_E * f_min_norm * np.log(f_max_norm / f_min_norm))

        return correlation

    @classmethod
    def within_event_correlation(cls, f_i, f_j):
        """
        Compute the within-event correlation between two frequencies.
        Uses equations 14 and 15 from the paper.
        """
        f_min = np.minimum(f_i, f_j)
        f_max = np.maximum(f_i, f_j)

        # Compute parameters
        eta_A = cls.compute_eta_A(f_min)
        gamma_A = cls.compute_gamma_A(f_min)

        # Base correlation (equation 14)
        rho_A0 = (1 - eta_A) * np.exp(gamma_A * np.log(f_max / f_min))

        # Full correlation including nugget effect (equation 15)
        if np.isclose(f_i, f_j):
            return 1.0
        else:
            rho_A = rho_A0 * (1 - np.exp(-cls.COEFF.nugget_exp * np.log(f_max / f_min)))
            return rho_A

    @classmethod
    def between_site_correlation(cls, f_i, f_j):
        """
        Compute the between-site correlation between two frequencies.
        Uses equations 14 and 15 from the paper with S parameters.
        """
        f_min = np.minimum(f_i, f_j)
        f_max = np.maximum(f_i, f_j)

        # Compute parameters
        eta_S = cls.compute_eta_S(f_min)
        gamma_S = cls.compute_gamma_S(f_min)

        # Base correlation (equation 14)
        rho_S0 = (1 - eta_S) * np.exp(gamma_S * np.log(f_max / f_min))

        # Full correlation including nugget effect (equation 15)
        if np.isclose(f_i, f_j):
            return 1.0
        else:
            rho_S = rho_S0 * (1 - np.exp(-cls.COEFF.nugget_exp * np.log(f_max / f_min)))
            return rho_S

    @classmethod
    def compute_variances(cls, freqs):
        """
        Compute standard deviations for between-event, between-site,
        and within-event terms.

        Parameters:
        -----------
        freqs : array-like
            Array of frequencies for which to compute standard deviations.

        Returns:
        --------
        tuple of arrays
            (sigma_E, sigma_S, sigma_A) arrays for each frequency
        """

        def compute_sigma_E(f):
            """Compute between-event standard deviation."""
            term1 = cls._sigmoid(
                f, cls.COEFF.sigma_E1, cls.COEFF.sigma_E2, cls.COEFF.sigma_E3
            )
            term2 = cls._sigmoid(
                f, cls.COEFF.sigma_E4, cls.COEFF.sigma_E5, cls.COEFF.sigma_E6
            )
            return cls.COEFF.sigma_E0 + term1 + term2

        def compute_sigma_S(f):
            """Compute between-site standard deviation."""
            term1 = cls._sigmoid(
                f, cls.COEFF.sigma_S1, cls.COEFF.sigma_S2, cls.COEFF.sigma_S3
            )
            term2 = cls._sigmoid(
                f, cls.COEFF.sigma_S4, cls.COEFF.sigma_S5, cls.COEFF.sigma_S6
            )
            return cls.COEFF.sigma_S0 + term1 + term2

        def compute_sigma_A(f):
            """Compute within-event standard deviation."""
            f_ref = cls.COEFF.sigma_A_FREF  # 5.0
            return (
                cls.COEFF.sigma_A0
                + cls.COEFF.sigma_A1 * np.log(np.maximum(f, f_ref) / f_ref) ** 2
            )

        # Compute standard deviations for all frequencies
        sigma_E = np.array([compute_sigma_E(f) for f in freqs])
        sigma_S = np.array([compute_sigma_S(f) for f in freqs])
        sigma_A = np.array([compute_sigma_A(f) for f in freqs])

        return sigma_E, sigma_S, sigma_A

    @classmethod
    def cov(cls, freqs, sigma_E=None, sigma_S=None, sigma_A=None, mag=6.0):
        """
        Compute the covariance matrix for Fourier spectral ordinates.

        Parameters:
        -----------
        freqs : array-like
            Array of frequencies for which to compute the covariance matrix.
        sigma_E : array-like, optional
            Between-event standard deviations for each frequency.
            If None, will be computed using model equations.
        sigma_S : array-like, optional
            Between-site standard deviations for each frequency.
            If None, will be computed using model equations.
        sigma_A : array-like, optional
            Within-event standard deviations for each frequency.
            If None, will be computed using model equations.
        mag: float, optional
            Earthquake magnitude (used for between-event correlations).
            Default is 6.0.

        Returns:
        --------
        cov_matrix : ndarray
            The covariance matrix for the Fourier spectral ordinates.
        """
        if sigma_E is None or sigma_S is None or sigma_A is None:
            # Compute standard deviations if not provided
            sigma_E, sigma_S, sigma_A = cls.compute_variances(freqs)

        # Check input dimensions
        n = len(freqs)
        if len(sigma_E) != n or len(sigma_S) != n or len(sigma_A) != n:
            raise ValueError("All input arrays must have the same length.")

        # Initialize the covariance matrix
        cov_matrix = np.zeros((n, n))

        # Compute the covariance matrix elements
        for i in range(n):
            for j in range(n):
                # Between-event contribution
                rho_E = cls.between_event_correlation(freqs[i], freqs[j], mag)
                cov_E = rho_E * sigma_E[i] * sigma_E[j]

                # Between-site contribution
                rho_S = cls.between_site_correlation(freqs[i], freqs[j])
                cov_S = rho_S * sigma_S[i] * sigma_S[j]

                # Within-event contribution
                rho_A = cls.within_event_correlation(freqs[i], freqs[j])
                cov_A = rho_A * sigma_A[i] * sigma_A[j]

                # Total covariance (equation 4)
                cov_matrix[i, j] = cov_E + cov_S + cov_A

        return cov_matrix

    @classmethod
    def cor(cls, freqs, sigma_E=None, sigma_S=None, sigma_A=None, mag=6.0):
        """
        Compute the correlation matrix for Fourier spectral ordinates.

        Parameters:
        -----------
        freqs : array-like
            Array of frequencies for which to compute the correlation matrix.
        sigma_E : array-like, optional
            Between-event standard deviations for each frequency.
            If None, will be computed using model equations.
        sigma_S : array-like, optional
            Between-site standard deviations for each frequency.
            If None, will be computed using model equations.
        sigma_A : array-like, optional
            Within-event standard deviations for each frequency.
            If None, will be computed using model equations.
        mag: float, optional
            Earthquake magnitude (used for between-event correlations).
            Default is 6.0.

        Returns:
        --------
        cor_matrix : ndarray
            The correlation matrix for the Fourier spectral ordinates.
        """
        cov = cls.cov(freqs, sigma_E, sigma_S, sigma_A, mag)
        stds = np.sqrt(np.diag(cov))
        cor = cov / np.outer(stds, stds)
        return cor
