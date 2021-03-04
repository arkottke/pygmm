"""Baker and Jayaram (2008, :cite:`baker08`) correlation model."""
import numpy as np

from .types import ArrayLike

__author__ = "Albert Kottke"


def calc_correls(periods: ArrayLike, period_cond: float) -> np.ndarray:
    """Baker and Jayaram (2008, :cite:`baker08`) correlation model.

    Parameters
    ----------
    periods : array_like
        Periods at which the correlation should be computed.
    period_cond : float
        Conditioning period

    Returns
    -------
    correls : :class:`np.ndarray`
        Correlation coefficients
    """
    period_min = np.minimum(periods, period_cond)
    period_max = np.maximum(periods, period_cond)

    c_1 = 1 - np.cos(
        np.pi / 2 - 0.366 * np.log(period_max / np.maximum(period_min, 0.109))
    )

    # The minimum() is added to prevent an overflow issue
    c_2 = np.select(
        [period_max < 0.2, True],
        [
            1
            - 0.105
            * (1 - 1 / (1 + np.exp(100 * np.minimum(period_max, 0.2) - 5)))
            * (period_max - period_min)
            / (period_max - 0.0099),
            0,
        ],
    )

    c_3 = np.select([period_max < 0.109, True], [c_2, c_1])

    c_4 = c_1 + 0.5 * (np.sqrt(c_3) - c_3) * (1 + np.cos(np.pi * period_min / 0.109))

    correls = np.select(
        [period_max < 0.109, period_min > 0.109, period_max < 0.200, True],
        [c_2, c_1, np.minimum(c_2, c_4), c_4],
    )

    return correls


def calc_cond_mean_spectrum(
    periods: ArrayLike,
    ln_psas: ArrayLike,
    ln_stds: ArrayLike,
    period_cond: float,
    ln_psa_cond: float,
) -> (np.ndarray, np.ndarray):
    """Conditional mean spectrum by Baker & Jayaram (2008, :cite:`baker08`).

    Parameters
    ----------
    periods : array_like
        Response spectral periods.
    ln_psas : array_like
        Natural logarithm of the 5%-damped spectral accelerations.
    ln_stds : array_like
        Logarithmic standard deviations.
    period_cond : float
        Conditioning period. This period does not need to be included in
        `periods`.
    ln_psa_cond : float
        Natural logarithm of the response at the conditioning period.

    Returns
    -------
    ln_psas_cms : :class:`np.ndarray`
        Natural logarithm of the conditional mean spectral accelerations.
    ln_stds_cms : :class:`np.ndarray`
        Logarithmic standard deviation of the conditional mean spectral
        acceleration.
    """
    periods = np.asarray(periods)
    ln_psas = np.asarray(ln_psas)
    ln_stds = np.asarray(ln_stds)

    correls = calc_correls(periods, period_cond)
    epsilon = (ln_psa_cond - np.interp(period_cond, periods, ln_psas)) / np.interp(
        period_cond, periods, ln_stds
    )

    ln_psas_cms = ln_psas + ln_stds * correls * epsilon
    ln_stds_cms = np.sqrt(ln_stds ** 2 * (1 - correls ** 2))
    return ln_psas_cms, ln_stds_cms
