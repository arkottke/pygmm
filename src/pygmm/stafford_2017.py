"""Stafford (2017, :cite:`safford2017`) correlation."""

import json
import os

import numpy as np

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
    COEFF = json.load(
        open(os.path.join(os.path.dirname(__file__), "data", "stafford_2017.json"))
    )

    @classmethod
    def compute_corner_frequency(cls, magnitude):
        """
        Compute the corner frequency based on magnitude.

        This implementation follows the model in the paper using the Yenier and Atkinson
        (2015) stress parameter model with a single-corner source spectrum.
        """
        # Approximate implementation based on common scaling relationships
        # A more accurate implementation would use the full YA15 model parameters
        return 10 ** (cls.COEFF.C1 - cls.COEFF.C2 * magnitude)

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

        # Sigmoid function (equation 13)
        def S(f, alpha0, alpha1, alpha2):
            return alpha0 / (1 + np.exp(-alpha2 * np.log(f / alpha1)))

        # Apply equation 12
        sigmoid_term = S(
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

        # Sigmoid function
        def S(f, alpha0, alpha1, alpha2):
            return alpha0 / (1 + np.exp(-alpha2 * np.log(f / alpha1)))

        # Apply equation 16
        term1 = S(f_min, cls.COEFF.eta_A0, cls.COEFF.eta_A1, cls.COEFF.eta_A2)
        term2 = S(f_min, cls.COEFF.eta_A3, cls.COEFF.eta_A4, cls.COEFF.eta_A5)

        return term1 * (1 + term2)

    @classmethod
    def compute_gamma_A(cls, f_min):
        """
        Compute gamma_A parameter for within-event correlations.
        Uses equation 17 from the paper.
        """
        # Apply equation 17
        min_f = np.minimum(f_min, 10)
        S_term = min_f / (min_f + cls.COEFF.gamma_A2 * (1 - min_f / 10))
        log_term = cls.COEFF.gamma_A3 * np.log(np.maximum(f_min, 10) / 10) ** 2

        return cls.COEFF.gamma_A0 + S_term * cls.COEFF.gamma_A1 + log_term

    @classmethod
    def compute_eta_S(cls, f_min):
        """
        Compute eta_S parameter for between-site correlations nugget effect.
        Uses equation 18 from the paper.
        """

        # Apply equation 18
        term1 = cls.COEFF.eta_S0 * np.log(np.maximum(np.minimum(f_min, 4), 0.25) / 0.25)
        term2 = cls.COEFF.eta_S1 * np.log(np.maximum(f_min, 4) / 4)

        return term1 + term2

    @classmethod
    def compute_gamma_S(cls, f_min):
        """
        Compute gamma_S parameter for between-site correlations.
        Uses equation 19 from the paper.
        """

        # Sigmoid function
        def S(f, alpha0, alpha1, alpha2):
            return alpha0 / (1 + np.exp(-alpha2 * np.log(f / alpha1)))

        # Apply equation 19
        sigmoid1 = S(f_min, cls.COEFF.gamma_S1, cls.COEFF.gamma_S2, cls.COEFF.gamma_S3)
        sigmoid2 = S(f_min, cls.COEFF.gamma_S4, cls.COEFF.gamma_S5, cls.COEFF.gamma_S6)

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
    def compute_variances(cls, frequencies):
        def compute_sigma_E(f):
            sigma_E0 = 0.8757
            sigma_E1 = -0.3309
            sigma_E2 = 0.5871
            sigma_E3 = 5.4264
            sigma_E4 = 0.5177
            sigma_E5 = 16.357
            sigma_E6 = 1.4689

            def S(f, alpha0, alpha1, alpha2):
                return alpha0 / (1 + np.exp(-alpha2 * np.log(f / alpha1)))

            return (
                sigma_E0
                + S(f, sigma_E1, sigma_E2, sigma_E3)
                + S(f, sigma_E4, sigma_E5, sigma_E6)
            )

        # Function to compute between-site standard deviation
        def compute_sigma_S(f):
            sigma_S0 = 0.6167
            sigma_S1 = -0.1495
            sigma_S2 = 0.7248
            sigma_S3 = 3.6985
            sigma_S4 = 0.3640
            sigma_S5 = 13.457
            sigma_S6 = 2.2497

            def S(f, alpha0, alpha1, alpha2):
                return alpha0 / (1 + np.exp(-alpha2 * np.log(f / alpha1)))

            return (
                sigma_S0
                + S(f, sigma_S1, sigma_S2, sigma_S3)
                + S(f, sigma_S4, sigma_S5, sigma_S6)
            )

        # Function to compute within-event standard deviation
        def compute_sigma_A(f):
            sigma_A0 = 0.7260
            sigma_A1 = 0.0328

            return sigma_A0 + sigma_A1 * np.log(np.maximum(f, 5) / 5) ** 2

        # Compute standard deviations for example frequencies
        sigma_E = np.array([compute_sigma_E(f) for f in frequencies])
        sigma_S = np.array([compute_sigma_S(f) for f in frequencies])
        sigma_A = np.array([compute_sigma_A(f) for f in frequencies])

        return sigma_E, sigma_S, sigma_A

    @classmethod
    def cov(cls, frequencies, sigma_E=None, sigma_S=None, sigma_A=None, magnitude=6.0):
        """
        Compute the covariance matrix for Fourier spectral ordinates.

        Parameters:
        -----------
        frequencies : array-like
            Array of frequencies for which to compute the covariance matrix.
        sigma_E : array-like
            Between-event standard deviations for each frequency.
        sigma_S : array-like
            Between-site standard deviations for each frequency.
        sigma_A : array-like
            Within-event standard deviations for each frequency.
        magnitude : float, optional
            Earthquake magnitude (used for between-event correlations).

        Returns:
        --------
        cov_matrix : ndarray
            The covariance matrix for the Fourier spectral ordinates.
        """
        if sigma_E is None or sigma_S is None or sigma_A is None:
            # Compute standard deviations if not provided
            sigma_E, sigma_S, sigma_A = cls.compute_variances(frequencies)
        # Check input dimensions
        n = len(frequencies)
        if len(sigma_E) != n or len(sigma_S) != n or len(sigma_A) != n:
            raise ValueError("All input arrays must have the same length.")

        # Initialize the covariance matrix
        cov_matrix = np.zeros((n, n))

        # Compute the covariance matrix elements
        for i in range(n):
            for j in range(n):
                # Between-event contribution
                rho_E = cls.between_event_correlation(
                    frequencies[i], frequencies[j], magnitude
                )
                cov_E = rho_E * sigma_E[i] * sigma_E[j]

                # Between-site contribution
                rho_S = cls.between_site_correlation(frequencies[i], frequencies[j])
                cov_S = rho_S * sigma_S[i] * sigma_S[j]

                # Within-event contribution
                rho_A = cls.within_event_correlation(frequencies[i], frequencies[j])
                cov_A = rho_A * sigma_A[i] * sigma_A[j]

                # Total covariance (equation 4)
                cov_matrix[i, j] = cov_E + cov_S + cov_A

        return cov_matrix

    def cor(cls, frequencies, sigma_E=None, sigma_S=None, sigma_A=None, magnitude=6.0):
        """
        Compute the correlation matrix for Fourier spectral ordinates.

        Parameters:
        -----------
        frequencies : array-like
            Array of frequencies for which to compute the correlation matrix.
        sigma_E : array-like
            Between-event standard deviations for each frequency.
        sigma_S : array-like
            Between-site standard deviations for each frequency.
        sigma_A : array-like
            Within-event standard deviations for each frequency.
        magnitude : float, optional
            Earthquake magnitude (used for between-event correlations).

        Returns:
        --------
        cor_matrix : ndarray
            The correlation matrix for the Fourier spectral ordinates.
        """
        cov = cls.cov(frequencies, sigma_E, sigma_S, sigma_A, magnitude)
        stds = np.sqrt(np.diag(cov))
        cor = cov / np.outer(stds, stds)
        return cor
