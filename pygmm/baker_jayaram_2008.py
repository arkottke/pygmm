import numpy as np

from . import ArrayLike


def calc_correl(periods: ArrayLike, period_cond: float) -> np.ndarray:
    """Calculate Baker and Jayaram (2008) correlation model.

    Parameters
    ----------
    periods: `array_like`
        Periods at which the correlation should be computed.
    period_cond: float
        Conditioning period

    Returns
    -------
    correl: :class:`np.ndarray`
        Correlation coefficient.
    """

    periods_min = np.minimum(periods, period_cond)
    periods_max = np.maximum(periods, period_cond)

    c_1 = (1 - np.cos(np.pi / 2 - 0.366 * np.log(periods_max / np.maximum(
        periods_min, 0.109))))

    c_2 = np.select([periods_max < 0.2, True], [
        1 - 0.105 * (1 - 1 / (1 + np.exp(100 * periods_max - 5))) *
        (periods_max - periods_min) / (periods_max - 0.0099), 0
    ])

    c_3 = np.select([periods_max < 0.109, True], [c_2, c_1])

    c_4 = (c_1 + 0.5 * (np.sqrt(c_3) - c_3) *
           (1 + np.cos(np.pi * periods_min / 0.109)))

    correl = np.select(
        [periods_max < 0.109, periods_min > 0.109, periods_max < 0.200, True],
        [c_2, c_1, np.minimum(c_2, c_4), c_4], )

    return correl


def calc_cond_mean_spectrum(periods: ArrayLike,
                            ln_psas: ArrayLike,
                            ln_stds: ArrayLike,
                            period_cond: float,
                            ln_psa_cond: float) -> (np.ndarray, np.ndarray):
    """Compute the conditional mean spectrum (CMS) and associated standard
    deviation.

    Parameters
    ----------
    periods: `array_like`
        Response spectral periods.
    ln_psas: `array_like`
        Natural logarithm of the 5%-damped spectral accelerations.
    ln_stds: `array_like`
        Logarithmic standard deviations.
    period_cond: float
        Conditioning period. This period does not need to be included in
        `periods`.
    ln_psa_cond: float
        Natural logarithm of the response at the conditioning period.

    Returns
    -------
    ln_psas_cms: :class:`np.ndarray`
        Natural logarithm of the conditional 5%-damped spectral accelerations.
        The spectral acceleration is equal to `ln_psa_cond` at `period_cond`.
    ln_stds_cms: :class:`np.ndarray`
        Logarithmic standard deviation of conditional spectrum.
    """
    periods = np.asarray(periods)
    ln_psas = np.asarray(ln_psas)
    ln_stds = np.asarray(ln_stds)

    correl = calc_correl(periods, period_cond)
    epsilon = ((ln_psa_cond - np.interp(period_cond, periods, ln_psas)) /
               np.interp(period_cond, periods, ln_stds))

    ln_psa_cms = ln_psas + ln_stds * correl * epsilon
    ln_stds_cms = np.sqrt(ln_stds ** 2 * (1 - correl ** 2))
    return ln_psa_cms, ln_stds_cms
