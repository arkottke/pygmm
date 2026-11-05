"""Afshari and Stewart (2016, :cite:p:`afshari2016`) duration model."""

import numpy as np

from .. import model
from ..contracts import Duration
from ..registry import register


@register(provides=("duration",))
class AfshariStewart2016(model.Model):
    """Afshari and Stewart (2016, :cite:p:`afshari2016`) duration model.

    The basin term uses the difference between `depth_1_0` and the average
    depth for the site's `v_s30`. The average depth uses the Japan relation for
    ``region="japan"``, and the California relation otherwise. If `depth_1_0`
    is *None* or NaN, no basin term is applied.

    The model is vectorized. Each scenario value (``mag``, ``dist_rup``,
    ``v_s30``, ``mechanism``, ``depth_1_0``, and ``region``) can be a scalar or
    an array, and the arrays are broadcast against each other. If
    ``depth_1_0`` is *None*, no basin term is applied to any of the scenarios.
    :attr:`duration` and :attr:`std_err` are record arrays with the fields
    ``D_5t75``, ``D_5t95``, and ``D_20t80``. For a scalar scenario, each field
    is a scalar, as before. For arrays of N scenarios, each field has shape
    (N,). Multidimensional scenario arrays keep their shape.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    Examples
    --------
    >>> import numpy as np
    >>> import pygmm
    >>> s = pygmm.Scenario(
    ...     mag=np.array([5.5, 6.5, 7.5]), dist_rup=np.array([5.0, 20.0, 80.0]),
    ...     v_s30=400.0, mechanism=np.array(["SS", "NS", "RS"]),
    ...     depth_1_0=np.array([0.1, 0.3, 0.6]))
    >>> m = pygmm.AfshariStewart2016(s)
    >>> m.duration.D_5t95.shape
    (3,)
    >>> m.std_err.D_5t95.shape
    (3,)

    """

    NAME = "Afshari and Stewart (2016)"
    ABBREV = "AS16"

    PARAMS = [
        model.NumericParameter("mag", True, 3, 7.9),
        model.NumericParameter("dist_rup", True, 0, 200),
        model.NumericParameter("v_s30", True, 200, 1000),
        model.CategoricalParameter("mechanism", True, ["SS", "NS", "RS"]),
        model.NumericParameter("depth_1_0", False),
        model.CategoricalParameter(
            "region", False, ["california", "global", "japan"], "california"
        ),
    ]

    def __init__(self, scenario):
        super().__init__(scenario)

        # scenario
        s = self._scenario

        # Scenario values get a trailing axis to broadcast against the
        # coefficients, which have one value per duration measure:
        # D_{5-75}, D_{5-95}, D_{20-80}
        mag = model.as_column(s.mag)
        dist_rup = model.as_column(s.dist_rup)
        v_s30 = model.as_column(s.v_s30)

        # source coefficients
        mag_1 = np.array([5.35, 5.2, 5.2])
        mag_2 = np.array([7.15, 7.4, 7.4])
        mag_star = 6
        # Coefficients by mechanism. Mechanisms that are not one of the options
        # (None after checking the scenario) use the unspecified coefficients.
        mechanism = model.as_column(s.mechanism)
        is_mech = [model.equals(mechanism, m) for m in ["NS", "RS", "SS"]]
        b_0 = np.select(
            is_mech,
            [
                np.array([1.555, 2.541, 1.409]),
                np.array([0.7806, 1.612, 0.7729]),
                np.array([1.279, 2.302, 0.8804]),
            ],
            np.array([1.280, 2.182, 0.8822]),
        )
        b_1 = np.select(
            is_mech,
            [
                np.array([4.992, 3.170, 4.778]),
                np.array([7.061, 4.536, 6.579]),
                np.array([5.578, 3.467, 6.188]),
            ],
            np.array([5.576, 3.628, 6.182]),
        )
        b_2 = np.array([0.9011, 0.9443, 0.7414])
        b_3 = np.array([-1.684, -3.911, -3.164])
        # stress drop indices for duration measures

        stress_drop = np.exp(
            b_1
            + b_2 * (np.minimum(mag, mag_2) - mag_star)
            + b_3 * np.maximum((mag - mag_2), 0)
        )
        # corner frequency
        f_0 = self._corner_freq(stress_drop)

        # source duration
        F_E = (mag <= mag_1) * b_0 + (mag > mag_1) * (1 / f_0)

        # path coeficients
        c_1 = np.array([0.1159, 0.3165, 0.0646])
        c_2 = np.array([0.1065, 0.2539, 0.0865])
        c_3 = np.array([0.0682, 0.0932, 0.0373])
        c_4 = np.array([-0.2246, -0.3183, -0.4237])
        c_5 = np.array([0.0006, 0.0006, 0.0005])
        R_1 = 10
        R_2 = 50
        # path duration
        F_P = (
            c_1 * np.minimum(dist_rup, R_1)
            + c_2 * np.maximum((np.minimum(dist_rup, R_2) - R_1), 0)
            + c_3 * np.maximum((dist_rup - R_2), 0)
        )

        # site coefficients
        V_1 = 600
        V_ref = np.array([368.2, 369.9, 369.6])
        dz_1ref = 200
        # depth to bedrock duration
        if s.depth_1_0 is not None:
            # Basin depth differential in m (Equation 10 and Figure 12 of the
            # paper). depth_1_0 and calc_depth_1_0 are in km. Table 2 lists
            # c_5 in 1/km, but Figure 12 shows that it applies to dz_1 in m,
            # consistent with dz_1ref = 200 m.
            dz_1 = model.as_column(
                1000 * (s.depth_1_0 - self.calc_depth_1_0(s.v_s30, s.region))
            )
            # Same as `dz_1 if dz_1 <= dz_1ref else dz_1ref`. Unknown depths
            # (NaN) have no basin term.
            F_dz1 = np.where(
                np.isnan(dz_1), 0.0, c_5 * np.where(dz_1 <= dz_1ref, dz_1, dz_1ref)
            )
        else:
            F_dz1 = 0
        # site duration
        F_S = c_4 * np.log(np.minimum(v_s30, V_1) / V_ref) + F_dz1

        # aleatory coefficients
        tau_1 = np.array([0.28, 0.25, 0.30])
        tau_2 = np.array([0.25, 0.19, 0.19])
        phi_1 = np.array([0.54, 0.43, 0.56])
        phi_2 = np.array([0.41, 0.35, 0.45])
        # aleatory variability
        tau = tau_1 + (tau_2 - tau_1) * (
            np.minimum(np.maximum(mag, 6.5), 7.0) - 6.5
        ) / (7.0 - 6.5)
        phi = phi_1 + (phi_2 - phi_1) * (
            np.minimum(np.maximum(mag, 5.5), 5.75) - 5.5
        ) / (5.75 - 5.50)

        # total druation
        self._ln_dur = np.log(F_E + F_P) + F_S

        # aleatory standard deviation
        self._std_err = np.sqrt(tau**2 + phi**2)

    @staticmethod
    def calc_depth_1_0(
        v_s30: float | np.ndarray, region: str | np.ndarray = "california"
    ) -> float | np.ndarray:
        """Calculate the average depth to 1 km/sec (:math:`Z_{1.0}`).

        Parameters
        ----------
        v_s30 : float or array_like
            time-averaged shear-wave velocity over the top 30 m of the site
            (:math:`V_{s30}`, m/s).
        region : str or array_like, optional
            basin region. The Japan relation is used for "japan", and the
            California relation otherwise. An array of regions is broadcast
            against `v_s30`.

        Returns
        -------
        depth_1_0 : float or :class:`np.ndarray`
            average depth to a shear-wave velocity of 1 km/sec [km]
        """
        if np.ndim(region) > 0:
            return np.where(
                model.equals(region, "japan"),
                AfshariStewart2016.calc_depth_1_0(v_s30, "japan"),
                AfshariStewart2016.calc_depth_1_0(v_s30, "california"),
            )
        if region in ["japan"]:
            # Japan
            power = 2
            v_ref = 412.39
            slope = -5.23 / power
        else:
            # Global
            power = 4
            v_ref = 570.94
            slope = -7.15 / power

        return (
            np.exp(
                slope
                * np.log((v_s30**power + v_ref**power) / (1360.0**power + v_ref**power))
            )
            / 1000
        )

    @property
    def duration(self):
        """Return the durations as a `np.recarray`.

        Values can be accessed with 'D_5t75', 'D_5t95', or 'D_20t80'.
        """
        return AfshariStewart2016._as_recarray(np.exp(self._ln_dur))

    @property
    def std_err(self):
        """Return the standard errors of the durations as a `np.recarray`.

        Values can be accessed with 'D_5t75', 'D_5t95', or 'D_20t80'.
        """
        std_err = self._std_err
        if np.ndim(self._ln_dur) > 1:
            std_err = np.broadcast_to(std_err, np.shape(self._ln_dur))
        return AfshariStewart2016._as_recarray(std_err)

    @staticmethod
    def _as_recarray(values):
        # The duration measures are along the last axis
        return np.rec.fromarrays(
            np.moveaxis(values, -1, 0),
            names=[
                "D_5t75",
                "D_5t95",
                "D_20t80",
            ],
        )

    def _corner_freq(self, stress_drop):
        """Corner frequency."""
        moment = model.as_column(10 ** (1.5 * self.scenario.mag + 16.05))
        return (4.9e6 * 3.2) * (stress_drop / moment) ** (1 / 3)

    def duration_model(self) -> Duration:
        """Return the :class:`~pygmm.contracts.Duration` contract.

        ``duration`` and ``std_err`` are record arrays over the three energy
        intervals; 5-95% is taken as the primary value.

        The values are floats for a scalar scenario and arrays for a vectorized
        scenario.
        """
        d, s = self.duration, self.std_err
        return Duration(
            duration=model.as_scalar_or_array(d["D_5t95"]),
            d_5_75=model.as_scalar_or_array(d["D_5t75"]),
            d_5_95=model.as_scalar_or_array(d["D_5t95"]),
            d_20_80=model.as_scalar_or_array(d["D_20t80"]),
            ln_std=model.as_scalar_or_array(s["D_5t95"]),
        )
