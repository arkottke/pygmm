"""Hale and Abrahamson (2019) flexible sliding-mass displacement of earth dams.

Requires :mod:`pyrotd` (``pip install pygmm[slope]``).
"""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

from .rigid_block import calc_rigid_disp

__author__ = "Albert Kottke"


class HaleAbrahamson19:
    """Hale and Abrahamson (2019) displacement of a flexible sliding mass.

    The input motion is filtered by a nonlinear transfer function of the dam
    before a rigid sliding-block analysis is performed.

    Provide exactly one of `shear_vel` and `freq_nat`; the other is computed
    from :math:`f_n = V_S / (3 H)`. This is the empirical relation of the
    model; it differs from the :math:`V_S / (2.61 H)` of a homogeneous shear
    wedge on a rigid base, see ``pystrata.dam.analytical.wedge_natural_freqs``.

    Parameters
    ----------
    height : float
        Height [ft].
    shear_vel : float, optional
        Shear-wave velocity [ft/s].
    freq_nat : float, optional
        Natural frequency [Hz].
    nl_model : str, default "darendeli"
        Nonlinear model of the material: ``"darendeli"`` or
        ``"vucetic_dobry"``.
    """

    C = np.rec.fromrecords(
        [
            (
                0.76,
                0.35,
                -0.00049,
                0.00329,
                -0.03558,
                0.35297,
                -0.89726,
                1.25375,
                0.00024,
                -1.14426,
                0.00090,
                -1.91824,
                0.00145,
                -0.00993,
                0.86560,
                0.57934,
                0.00172,
                0.40392,
                0.85,
            )
        ],
        names="a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,alpha",
    )[0]

    NL_MODELS = ("darendeli", "vucetic_dobry")

    def __init__(
        self,
        height: float,
        shear_vel: float | None = None,
        freq_nat: float | None = None,
        nl_model: str = "darendeli",
    ):
        if nl_model not in self.NL_MODELS:
            raise ValueError(f"nl_model must be one of: {self.NL_MODELS}")

        self._height = height
        self._nl_model = nl_model

        if shear_vel is None and freq_nat is not None:
            self._freq_nat = freq_nat
            self._shear_vel = 3 * height * freq_nat
        elif freq_nat is None and shear_vel is not None:
            self._shear_vel = shear_vel
            self._freq_nat = shear_vel / (3 * height)
        else:
            raise NotImplementedError("Provide one of `shear_vel` or `freq_nat`.")

    @property
    def freq_nat(self) -> float:
        return self._freq_nat

    @property
    def height(self) -> float:
        return self._height

    @property
    def shear_vel(self) -> float:
        return self._shear_vel

    def __call__(
        self, time_step: float, accels: npt.ArrayLike, yield_coef: float
    ) -> list[float]:
        """Compute the sliding displacements.

        Parameters
        ----------
        time_step : float
            Time step [s].
        accels : array_like
            Acceleration time series [g].
        yield_coef : float
            Yield coefficient [g].

        Returns
        -------
        list of float
            Final displacement [cm] for the inverted and the original
            polarity of the motion.
        """
        sliding_accels = self._calc_sliding_mass_accels(time_step, accels)

        return [
            calc_rigid_disp(time_step, sliding_accels, yield_coef, invert=invert)[
                "disps"
            ][-1]
            for invert in [True, False]
        ]

    def _calc_sliding_mass_accels(self, time_step: float, accels: npt.ArrayLike):
        try:
            import pyrotd
        except ImportError as err:
            raise ImportError(
                "HaleAbrahamson19 requires pyrotd: pip install pygmm[slope]"
            ) from err

        # Compute the Fourier amplitude spectra
        accels = np.asarray(accels)
        fourier_amps = np.fft.rfft(accels)
        freqs = np.fft.rfftfreq(accels.size, d=time_step)
        ln_spec_accel = np.log(
            pyrotd.calc_oscillator_resp(
                freqs,
                fourier_amps,
                0.05,
                1 / self.period_nl,
                peak_resp_only=True,
                osc_type="psa",
            )
        )

        trans_func = self._calc_trans_func(freqs, ln_spec_accel)
        fourier_amps *= trans_func

        return np.fft.irfft(fourier_amps)

    @property
    def period_nl(self) -> float:
        r"""Nonlinear period [s] of the dam (Eq. 4.3 of Hale, 2019).

        This is the period at which the 5% damped pseudo-spectral acceleration
        is computed, :math:`\ln T_{NL} = a_1 \ln T_o + a_2`, where
        :math:`T_o = 1 / f_n` is the small-strain period.
        """
        return float(np.exp(self.C.a1 * np.log(1 / self.freq_nat) + self.C.a2))

    def _calc_params(self, ln_spec_accel):
        """Compute the damping ratio and the first-mode period ratio T1 / To.

        Parameters
        ----------
        ln_spec_accel : float
            Natural log of the 5% damped pseudo-spectral acceleration [g] at
            the nonlinear period.

        Returns
        -------
        damping : float
            Damping ratio (Eq. 4.9).
        period_ratio : float
            Ratio of the first-mode period to the small-strain period
            (Eq. 4.4).
        """
        ln_height = np.log(self.height)

        if self._nl_model == "darendeli":
            c3 = 1.0
            b5 = 0.0
            b6 = 0.0
        else:
            c3 = 0.32
            b5 = -2.0
            b6 = -1.0

        C = self.C
        b1 = C.a9 * self.shear_vel + C.a10
        b2 = C.a11 * self.shear_vel + C.a12
        b3 = (C.a13 * ln_height + C.a14) * self.shear_vel + (C.a15 * ln_height + C.a16)
        b4 = C.a17 * self.shear_vel + C.a18
        damping = np.exp(b1 + b2 / (1 + np.exp(b3 + b5 + (b4 + b6) * ln_spec_accel)))

        c1 = (C.a3 * ln_height + C.a4) * self.shear_vel + (C.a5 * ln_height + C.a6)
        c2 = C.a7 * c1 + C.a8
        period_ratio = 1 + c3 * np.exp(c1 * (c2 + ln_spec_accel))

        return damping, period_ratio

    def _calc_trans_func(self, freqs: npt.ArrayLike, ln_spec_accel):
        """Transfer function to the average acceleration of the sliding mass.

        Equation 4.1 of Hale (2019), with the natural frequency from
        Equations 4.2 and 4.4.
        """
        freqs = np.asarray(freqs)
        damping, period_ratio = self._calc_params(ln_spec_accel)

        # First-mode (damped) frequency of the dam, shortened by nonlinearity
        freq_first = self.freq_nat / period_ratio
        # Natural frequency of the undamped system (Eq. 4.2)
        nat_freq_eff = freq_first / np.sqrt(1 - damping**2)

        trans_func = self.C.alpha * (
            1
            + (1j * freqs) ** 2
            / (freqs**2 - nat_freq_eff**2 - 2j * damping * nat_freq_eff * freqs)
        )
        # Amplitudes are at least one below the first-mode frequency
        mask = (freqs < freq_first) & (np.abs(trans_func) < 1.0)
        trans_func[mask] = 1.0

        return trans_func
