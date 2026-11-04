"""SPT-based liquefaction triggering.

Implements the SPT-based procedure of Boulanger and Idriss (2014): cyclic
stress ratio (CSR), clean-sand equivalent blow count
:math:`(N_1)_{60cs}`, magnitude scaling factor (MSF), and cyclic resistance
ratio (CRR). Stresses are in atmospheres throughout.
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from ._utils import PRESS_ATM, UNIT_WT_WATER

__author__ = "Albert Kottke"


def calc_stresses(
    depth_wt: float, depth: npt.ArrayLike, unit_wt: npt.ArrayLike
) -> tuple[npt.NDArray[np.floating], npt.NDArray[np.floating]]:
    """Calculate the stresses within the profile.

    Parameters
    ----------
    depth_wt : float
        Depth to the water table [m].
    depth : array_like
        Depth [m] at the top of each increment.
    unit_wt : array_like
        Unit weight [kN/m³] of each increment.

    Returns
    -------
    stress_v : :class:`numpy.ndarray`
        Vertical stress [atm].
    stress_v_eff : :class:`numpy.ndarray`
        Vertical effective stress [atm].
    """
    depth = np.asarray(depth, dtype=float)
    unit_wt = np.asarray(unit_wt, dtype=float)

    press_pore = np.maximum(depth - depth_wt, 0) * UNIT_WT_WATER / PRESS_ATM

    stress_v = np.zeros_like(depth)
    stress_v[0] = depth[0] * unit_wt[0]

    if depth_wt < 0:
        stress_v += -depth_wt * unit_wt

    # Compute the average unit weight assuming that each depth increment
    # includes half of the unit weight of the material in the two increment.
    avg_unit_wt = np.convolve(unit_wt, [0.5, 0.5], mode="valid")
    thickness = np.diff(depth)
    stress_v[1:] = stress_v[0] + np.cumsum(thickness * avg_unit_wt)
    stress_v /= PRESS_ATM

    stress_v_eff = stress_v - press_pore
    return stress_v, stress_v_eff


def calc_csr(
    mag: float,
    pga: float,
    depth: npt.ArrayLike,
    stress_v: npt.ArrayLike,
    stress_v_eff: npt.ArrayLike,
) -> npt.NDArray[np.floating]:
    """Compute the cyclic stress ratio.

    Parameters
    ----------
    mag : float
        Moment magnitude.
    pga : float
        Peak ground acceleration [g].
    depth : array_like
        Depth [m] at the top of each increment.
    stress_v : array_like
        Vertical stress [atm].
    stress_v_eff : array_like
        Vertical effective stress [atm].

    Returns
    -------
    :class:`numpy.ndarray`
        Cyclic stress ratio.
    """
    depth = np.asarray(depth, dtype=float)
    a = -1.012 - 1.126 * np.sin(depth / 11.73 + 5.133)
    b = 0.106 + 0.118 * np.sin(depth / 11.28 + 5.142)
    r_d = np.exp(a + b * mag)

    return 0.65 * np.asarray(stress_v) / np.asarray(stress_v_eff) * pga * r_d


def calc_n_1_60cs(
    depth: npt.ArrayLike,
    stress_v_eff: npt.ArrayLike,
    n_60: npt.ArrayLike,
    fc: npt.ArrayLike,
) -> npt.NDArray[np.floating]:
    """Compute the clean-sand equivalent, overburden-corrected blow count.

    Parameters
    ----------
    depth : array_like
        Depth [m].
    stress_v_eff : array_like
        Vertical effective stress [atm].
    n_60 : array_like
        Blow count corrected to 60% energy, :math:`N_{60}`.
    fc : array_like
        Fines content [%].

    Returns
    -------
    :class:`numpy.ndarray`
        :math:`(N_1)_{60cs}`.
    """
    n_60 = np.asarray(n_60, dtype=float)
    stress_v_eff = np.asarray(stress_v_eff, dtype=float)
    fc = np.asarray(fc, dtype=float)
    depth = np.asarray(depth, dtype=float)

    # Initial estimate
    n_1_60cs = n_60
    c_r = np.interp(
        depth,
        [3, 3.5, 5, 7.5, 10],
        [0.75, 0.80, 0.85, 0.95, 1.00],
        left=0.75,
        right=1.0,
    )
    # Fines content correction - Equation 2.23
    d_fc = np.exp(1.63 + 9.7 / (fc + 0.01) - (15.7 / (fc + 0.01)) ** 2)

    while True:
        # Overburden correction - Equation 2.15
        m = 0.784 - 0.0768 * np.sqrt(np.minimum(n_1_60cs, 46))
        c_n = np.minimum(stress_v_eff**-m, 1.7)

        prev = n_1_60cs
        n_1_60cs = (c_r * c_n * n_60) + d_fc
        if np.allclose(prev, n_1_60cs):
            break
    return n_1_60cs


def calc_msf(mag: float, n_1_60cs: npt.ArrayLike) -> npt.NDArray[np.floating]:
    """Compute the magnitude scaling factor."""
    n_1_60cs = np.asarray(n_1_60cs, dtype=float)
    msf_max = np.minimum(1.09 + (n_1_60cs / 31.5) ** 2, 2.2)
    return 1 + (msf_max - 1) * (8.64 * np.exp(-mag / 4) - 1.325)


def calc_crr(
    mag: float, stress_v_eff: npt.ArrayLike, n_1_60cs: npt.ArrayLike
) -> npt.NDArray[np.floating]:
    """Compute the cyclic resistance ratio, adjusted to `mag` and `stress_v_eff`.

    Parameters
    ----------
    mag : float
        Moment magnitude.
    stress_v_eff : array_like
        Vertical effective stress [atm].
    n_1_60cs : array_like
        Clean-sand equivalent blow count, :math:`(N_1)_{60cs}`.

    Returns
    -------
    :class:`numpy.ndarray`
        Cyclic resistance ratio.
    """
    n_1_60cs = np.asarray(n_1_60cs, dtype=float)
    msf = calc_msf(mag, n_1_60cs)

    c_sigma = np.minimum(1 / (18.9 - 2.55 * np.sqrt(n_1_60cs)), 0.3)
    k_sigma = np.minimum(1 - c_sigma * np.log(stress_v_eff), 1.1)

    crr_75 = np.exp(
        n_1_60cs / 14.1
        + (n_1_60cs / 126) ** 2
        - (n_1_60cs / 23.6) ** 3
        + (n_1_60cs / 25.4) ** 4
        - 2.8
    )

    return crr_75 * msf * k_sigma


def calc_triggering(df, mag: float, pga: float, depth_wt: float) -> None:
    """Evaluate liquefaction triggering for an SPT profile, in place.

    Parameters
    ----------
    df : :class:`pandas.DataFrame`
        Profile with columns ``depth`` [m], ``unit_wt`` [kN/m³], ``N``
        (:math:`N_{60}`), and ``FC`` (fines content, %). The columns
        ``stress_v``, ``stress_v_eff``, ``csr``, ``N_1_60cs``, ``crr``, and
        ``fs_liq`` (factor of safety) are added.
    mag : float
        Moment magnitude.
    pga : float
        Peak ground acceleration [g].
    depth_wt : float
        Depth to the water table [m].
    """
    df["stress_v"], df["stress_v_eff"] = calc_stresses(
        depth_wt, df["depth"], df["unit_wt"]
    )

    df["csr"] = calc_csr(mag, pga, df["depth"], df["stress_v"], df["stress_v_eff"])
    df["N_1_60cs"] = calc_n_1_60cs(df["depth"], df["stress_v_eff"], df["N"], df["FC"])

    df["crr"] = calc_crr(mag, df["stress_v_eff"], df["N_1_60cs"])

    df["fs_liq"] = df["crr"] / df["csr"]
