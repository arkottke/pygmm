"""Basic models."""

import collections
import os
import warnings
from collections.abc import Sequence
from typing import Any

import numpy as np
import numpy.typing as npt
from scipy.interpolate import interp1d

from ._types import ArrayLike, InterpKind
from .contracts import ResponseSpectrum


class Scenario(collections.UserDict):
    r"""An eathquake scenario used in all ground motion models.

    Parameters
    ----------
    depth_1_0 : float
        depth to the 1.0 km∕s shear-wave velocity horizon beneath the site,
        :math:`Z_{1.0}` in [km].
    depth_2_5 : float
        depth to the 2.5 km∕s shear-wave velocity horizon beneath the site,
        :math:`Z_{2.5}` in [km].
    depth_tor : float
        depth to the top of the rupture plane (:math:`Z_{tor}`, km).
    depth_bor : float
        depth to the bottom of the rupture plane (:math:`Z_{bor}`, km).
    depth_bot : float
        depth to bottom of seismogenic crust [km].
    depth_sed : float
        thickness of the Atlantic and Gulf coastal plain sediments beneath the
        site [km], the depth to the base of the sediments (``zSed`` in USGS
        nshmp-lib). NaN or *None* is a site that is not on the coastal plain.
    dip : float
        fault dip angle (:math:`\phi`, deg).
    dist_jb : float
        Joyner-Boore distance to the rupture plane (:math:`R_\text{JB}`, km)
    dist_crjb : float
        centroid Joyner-Boore distance, which is the shortest distance between
        the centroid of Joyner-Boore rupture surface of the potential Class 2
        earthquakes and the closest point on the edge of the Joyner-Boore
        rupture surface of the main shock (:math:`CR_\text{JB}`, km)
    dist_epi : float
        epicentral distance to the rupture plane (:math:`R_\text{epi}`, km)
    dist_hyp : float
        hypocentral distance to the rupture plane (:math:`R_\text{hyp}`, km).
    dist_rup : float
        closest distance to the rupture plane (:math:`R_\text{rup}`, km)
    dist_x : float
        site coordinate measured perpendicular to the fault strike from the
        fault line with the down-dip direction being positive (:math:`R_x`,
        km).
    dist_y0 : float
        horizontal distance off the end of the rupture measured parallel to
        strike (:math:`R_{y0}`, km).
    dpp_centered : float
        direct point parameter (DPP) for directivity effect (see Chiou and
        Spudich (2014, :cite:p:`spudich14`)) centered on the earthquake-specific
        average DPP for California.
    event_type : str
        event type. Type of event used in subduction models to distinguish
        between intraslab and interface events.
    is_aftershock : bool
        if the scenario is an aftershock.
    mag : float
        moment magnitude of the event (:math:`M_w`)
    mechanism : str
        fault mechanism. Valid options: "SS", "NS", "RS", and "U". See
        the mechanism reference in the user guide for more information.
    on_hanging_wall : bool
        If the site is located on the hanging wall of the fault. If *None*,
        then *False* is assumed.
    pga_ref : float
        peak ground accelearion in *g* at the model-specific reference
        condition.
    region : str
        region. Valid options are specified in a specific GMM.
    site_cond : str
        site condition. String description of the site condition. Valid
        options are specified in a specific GMM.
    tectonic_region : str
        tectonic region. Tectonic setting of the site typically used in
        subductin models.
    v_s30 : float
        time-averaged shear-wave velocity over the top 30 m of the site
        (:math:`V_{s30}`, m/s).
    vs_source : str
        source of the `v_s30` value.  Valid options include: "measured",
        "inferred"
    width : float
        down-dip width of the fault.

    """

    KNOWN_KEYS = [
        "depth_1_0",
        "depth_2_5",
        "depth_tor",
        "depth_bor",
        "depth_bot",
        "depth_hyp",
        "depth_sed",
        "dip",
        "dist_crjb",
        "dist_jb",
        "dist_epi",
        "dist_hyp",
        "dist_rup",
        "dist_x",
        "dist_y0",
        "dpp_centered",
        "event_type",
        "is_aftershock",
        "mag",
        "mechanism",
        "on_hanging_wall",
        "pga_ref",
        "region",
        "site_cond",
        "tectonic_region",
        "v_s30",
        "vs_source",
        "width",
    ]

    def __init__(self, **kwds):
        """Initialize the scenario."""
        super().__init__(kwds)
        self._check_keys(self.keys())

    def __getattr__(self, item):
        """Access the data with attributes."""
        return self.data[item]

    def __repr__(self):
        """Representation."""
        return "<Scenario(mag={mag}, dist_jb={dist_jb})>".format(**self.data)

    def copy_with(self, **kwds):
        """Return a copy of the scenario with parameters replaced.

        Parameters
        ----------
        **kwds
            New values of the parameters.
        """
        self._check_keys(kwds.keys())
        other = self.copy()
        other.update(**kwds)
        return other

    def _check_keys(self, keys):
        for k in keys:
            if k not in self.KNOWN_KEYS:
                raise Warning(f"{k} is not a recognized scenario key!")


class Model:
    """Base class of the ground motion models.

    A model is constructed from a :class:`Scenario`, or the keywords of one.
    Only the parameters in ``PARAMS`` are kept, and each is checked against the
    range recommended for the model.
    """

    #: Long name of the model
    NAME = ""
    #: Short name of the model
    ABBREV = ""
    #: Limits of model applicability
    LIMITS: dict = dict()
    #: Model parameters
    PARAMS: list = []

    def __init__(self, *args, **kwargs):
        """Initialize the model."""
        super().__init__()

        if len(args) == 1:
            scenario = args[0]
        else:
            scenario = Scenario(**kwargs)

        # Select the used parameters and check them against the recommended
        # values
        self._scenario = Scenario(
            **{p.name: scenario.get(p.name, None) for p in self.PARAMS}
        )
        self._check_inputs()

    def _check_inputs(self):
        for p in self.PARAMS:
            self._scenario[p.name] = p.check(self._scenario[p.name])

    @property
    def scenario(self):
        """Scenario of the model, limited to the parameters that the model uses."""
        return self._scenario


class GroundMotionModel(Model):
    """Abstract class for ground motion prediction models.

    Models that support the ``ims`` argument compute only the requested intensity
    measures, which is faster for large vectorized scenarios. With ``ims=None``
    (the default), all intensity measures are computed. The intensity measures
    are:

        +----------------------+------------------------------------------------+
        | Name                 | Description                                    |
        +======================+================================================+
        | ``"pga"``            | peak ground acceleration                       |
        +----------------------+------------------------------------------------+
        | ``"pgv"``            | peak ground velocity                           |
        +----------------------+------------------------------------------------+
        | ``"pgd"``            | peak ground displacement                       |
        +----------------------+------------------------------------------------+
        | ``"psa_1p000"``      | pseudo-spectral acceleration at one period     |
        |                      | (here 1.0 s), with "p" as the decimal point    |
        +----------------------+------------------------------------------------+
        | ``"psa_ngawest2_21"``| pseudo-spectral acceleration at the 21 periods |
        |                      | used to compare the NGA-West2 models           |
        +----------------------+------------------------------------------------+
        | ``"psa_all"``        | pseudo-spectral acceleration at all of the     |
        |                      | model's periods                                |
        +----------------------+------------------------------------------------+

    ``periods``, ``spec_accels``, and ``ln_stds`` contain only the computed
    spectral periods, in order of increasing period, and ``psa_ims`` gives their
    names.
    """

    #: Intensity measures that can be requested with the ``ims`` argument, in
    #: addition to individual periods such as "psa_1p000"
    IMS = ("pga", "pgv", "pgd", "psa_ngawest2_21", "psa_all")
    #: The 21 spectral periods (s) used to compare the NGA-West2 models
    PERIODS_NGAWEST2_21 = np.array(
        [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4]
        + [0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]
    )

    #: Indices for the spectral accelerations
    INDICES_PSA: np.ndarray = np.array([])
    #: Indices of the periods
    PERIODS: np.ndarray = np.array([])
    #: Index of the peak ground acceleration
    INDEX_PGA: int | None = None
    #: Index of the peak ground velocity
    INDEX_PGV: int | None = None
    #: Index of the peak ground displacement
    INDEX_PGD: int | None = None
    #: Scale factor to apply to get PGV in cm/sec
    PGV_SCALE = 1.0
    #: Scale factor to apply to get PGD in cm
    PGD_SCALE = 1.0

    def __init__(self, scenario: Scenario, ims: Sequence[str] | None = None):
        """Initialize the model."""
        super().__init__(scenario)

        self._ln_resp: np.ndarray | None = None
        self._ln_std: np.ndarray | None = None
        # Branches of the logic tree of the ground motion distribution (e.g., the
        # USGS epistemic branches of the median) as a list of (weight, ln_resp,
        # ln_std), or None for a single branch (_ln_resp and _ln_std)
        self._branches: list | None = None
        self._ims: tuple[str, ...] | None = None
        # Coefficient rows (periods) that are computed, or None for all rows
        self._indices: np.ndarray | None = None
        # Coefficient rows of the computed spectral accelerations
        self._psa_indices: np.ndarray
        if ims is None:
            self._psa_indices = np.asarray(self.INDICES_PSA)
        else:
            self._ims = self._check_ims(ims)
            self._indices = self._ims_indices(self._ims)
            psa = [im for im in self._ims if im.startswith("psa_")]
            self._psa_indices = (
                self._ims_indices(psa) if psa else np.array([], dtype=int)
            )

    def _check_ims(self, ims: str | Sequence[str]) -> tuple[str, ...]:
        """Check the requested intensity measures."""
        ims = (ims,) if isinstance(ims, str) else tuple(ims)
        if not ims:
            raise ValueError("ims must include at least one intensity measure")
        for im in ims:
            # Raises an error for names that are not valid or not provided
            self._im_indices(im)
        return ims

    def _im_indices(self, im: str) -> np.ndarray:
        """Coefficient rows needed for an intensity measure."""
        if not isinstance(im, str):
            raise ValueError(f"{im!r} is not a valid intensity measure")
        if im in ("pga", "pgv", "pgd"):
            index = {
                "pga": self.INDEX_PGA,
                "pgv": self.INDEX_PGV,
                "pgd": self.INDEX_PGD,
            }[im]
            if index is None:
                raise ValueError(f"{self.NAME} does not provide {im!r}")
            return np.atleast_1d(index)

        indices_psa = np.asarray(self.INDICES_PSA, dtype=int)
        periods_psa = np.asarray(self.PERIODS)[indices_psa]
        if im == "psa_all":
            if not indices_psa.size:
                raise ValueError(f"{self.NAME} does not provide {im!r}")
            return indices_psa
        if im == "psa_ngawest2_21":
            periods = self.PERIODS_NGAWEST2_21
        elif im.startswith("psa_"):
            periods = np.atleast_1d(self.psa_period(im))
        else:
            raise ValueError(
                f"{im!r} is not a valid intensity measure. Valid options are: "
                + ", ".join(self.IMS)
                + ", or a spectral period such as 'psa_1p000'"
            )

        found = np.isclose(periods_psa[:, np.newaxis], periods, rtol=1e-6, atol=0)
        missing = periods[~found.any(axis=0)]
        if missing.size:
            raise ValueError(
                f"{self.NAME} does not provide {im!r}. Missing periods (s): "
                + ", ".join(f"{p:g}" for p in missing)
                + ". Available spectral periods are: "
                + ", ".join(self.psa_name(p) for p in periods_psa)
            )
        return indices_psa[found.any(axis=1)]

    def _ims_indices(self, ims) -> np.ndarray:
        """Sorted coefficient rows needed for the intensity measures."""
        return np.unique(np.concatenate([self._im_indices(im) for im in ims]))

    @staticmethod
    def psa_name(period: float) -> str:
        """Name of the spectral acceleration intensity measure at a period.

        Parameters
        ----------
        period : float
            spectral period [s]

        Returns
        -------
        name : str
            intensity measure name, e.g., "psa_1p000" for 1.0 s

        """
        return "psa_" + f"{period:.3f}".replace(".", "p")

    @staticmethod
    def psa_period(name: str) -> float:
        """Spectral period of a spectral acceleration intensity measure name.

        Parameters
        ----------
        name : str
            intensity measure name with "p" as the decimal point, e.g.,
            "psa_1p000", "psa_1p0", or "psa_0p075". Any number of decimals can
            be used.

        Returns
        -------
        period : float
            spectral period [s]

        """
        try:
            text = name[len("psa_") :]
            if not name.startswith("psa_") or not text or text.count("p") > 1:
                raise ValueError
            period = float(text.replace("p", "."))
        except ValueError:
            raise ValueError(
                f"{name!r} is not a valid spectral acceleration name. Use, e.g., "
                "'psa_1p000' for 1.0 s"
            ) from None
        return period

    @property
    def psa_ims(self) -> list[str]:
        """Names of the computed spectral accelerations, e.g., "psa_1p000"."""
        return [self.psa_name(p) for p in self.periods]

    def _coeff_rows(self, values: Any) -> Any:
        """Select the computed coefficient rows (periods) from per-period values.

        Models that support ``ims`` use this on coefficients and periods so that
        only the requested intensity measures are computed. The result has the
        type of `values` (typed as ``Any``), so that the fields of coefficient
        record arrays (e.g., ``c.e_1``) are visible to type checkers.
        """
        if self._indices is None:
            return values
        if not isinstance(values, np.ndarray):
            values = np.asarray(values)
        # Indexing keeps the array type, so coefficient record arrays keep
        # attribute access (e.g., ``c.e_1``)
        return values[self._indices]

    def _take(self, values: np.ndarray, index, im: str) -> np.ndarray:
        """Select periods from computed values by coefficient row."""
        if self._indices is None:
            return take_periods(values, index)
        cols = np.searchsorted(self._indices, index)
        cols_clipped = np.clip(cols, 0, len(self._indices) - 1)
        if not np.all(self._indices[cols_clipped] == index):
            raise ValueError(
                f"{im} was not computed. The model was created with "
                f"ims={list(self._ims or ())}; include {im!r} in ims."
            )
        return take_periods(values, cols)

    def interp_ln_spec_accels(
        self, periods: ArrayLike, kind: InterpKind = "linear"
    ) -> np.ndarray:
        """Interpolate the spectral acceleration.

        Interpolation of the spectral acceleration is done in natural log
        space.

        Parameters
        ----------
        periods : array_like
            spectral periods to interpolate the response.
        kind : str, optional
            see :func:`scipy.interpolate.interp1d` for description of kind.
            Options include: 'linear' (default), 'nearest', 'zero', 'slinear',
            'quadratic', and 'cubic'

        Returns
        -------
        ln_spec_accels : np.ndarray
            interpolated spectral accelerations

        """
        return interp1d(
            np.log(self.periods),
            self._take_psa(self._resp_vector),
            kind=kind,
            copy=False,
            bounds_error=False,
            fill_value=np.nan,
        )(np.log(periods))

    def interp_spec_accels(
        self, periods: ArrayLike, kind: InterpKind = "linear"
    ) -> np.ndarray:
        """Interpolate the spectral acceleration.

        Interpolation of the spectral acceleration is done in natural log
        space.

        Parameters
        ----------
        periods : array_like
            spectral periods to interpolate the response.
        kind : str, optional
            see :func:`scipy.interpolate.interp1d` for description of kind.
            Options include: 'linear' (default), 'nearest', 'zero', 'slinear',
            'quadratic', and 'cubic'

        Returns
        -------
        spec_accels : np.ndarray
            interpolated spectral accelerations

        """
        return np.exp(self.interp_ln_spec_accels(periods, kind))

    def interp_ln_stds(
        self, periods: ArrayLike, kind: InterpKind = "linear"
    ) -> np.ndarray:
        r"""Interpolate the logarithmic standard deviation.

        Interpolate the logarithmic standard deviation (:math:`\sigma_{\ln}`)
        of spectral acceleration at the provided damping at specified periods.

        Parameters
        ----------
        periods : array_like
            spectral periods to interpolate the response.
        kind : str, optional
            see :func:`scipy.interpolate.interp1d` for description of kind.
            Options include: 'linear' (default), 'nearest', 'zero', 'slinear',
            'quadratic', and 'cubic'

        Returns
        -------
        ln_stds : np.ndarray
            interpolated logarithmic standard deviations

        """
        if self._ln_std is None:
            raise NotImplementedError
        else:
            return interp1d(
                np.log(self.periods),
                self._take_psa(self._std_vector),
                kind=kind,
                copy=False,
                bounds_error=False,
                fill_value=np.nan,
            )(np.log(periods))

    @property
    def periods(self) -> np.ndarray:
        """Periods of the computed spectral accelerations.

        These are all of the periods specified by the model unless ``ims``
        selects specific periods.
        """
        self._check_psa_computed()
        return self.PERIODS[self._psa_indices]

    @property
    def spec_accels(self) -> np.ndarray:
        """Pseudo-spectral accelerations computed by the model [g]."""
        return np.exp(self._take_psa(self._resp_vector))

    def response_spectrum(self, damping: float = 0.05) -> ResponseSpectrum:
        """Return the computed PSA as a response spectrum.

        The result is duck-typed for cross-package interop, matching
        :class:`~pygmm.contracts.ResponseSpectrum`, so it can be handed to
        consumers such as ``pyrvt.motions.CompatibleRvtMotion`` without those
        packages depending on pygmm.

        Parameters
        ----------
        damping : float, optional
            Fractional damping of the oscillator. Ground motion models are
            published for 5% damping, which is the default.

        Returns
        -------
        response_spectrum : :class:`~pygmm.contracts.ResponseSpectrum`
            Periods, spectral accelerations, and damping.

        """
        return ResponseSpectrum(
            periods=self.periods,
            spec_accels=self.spec_accels,
            damping=damping,
        )

    @property
    def ln_stds(self) -> np.ndarray:
        """Pseudo-spectral accelerations log-standard deviation."""
        if self._ln_std is None:
            raise NotImplementedError
        else:
            return self._take_psa(self._std_vector)

    @property
    def pga(self) -> float | np.ndarray:
        """Peak ground acceleration (PGA) computed by the model [g]."""
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGA, "pga")

    @property
    def ln_pga(self) -> float | np.ndarray:
        """Natural logarithm of the peak ground acceleration (PGA) in g.

        Equal to ``np.log(pga)``, but without the exponential and logarithm, which
        is faster for large vectorized scenarios.
        """
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._take(self._resp_vector, self.INDEX_PGA, "pga")

    @property
    def ln_std_pga(self) -> float | np.ndarray:
        """Peak ground accelaration log-standard deviation."""
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._take(self._std_vector, self.INDEX_PGA, "pga")

    @property
    def pgv(self) -> float | np.ndarray:
        """Peak ground velocity (PGV) computed by the model [cm/sec]."""
        if self.INDEX_PGV is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGV, "pgv") * self.PGV_SCALE

    @property
    def ln_std_pgv(self) -> float | np.ndarray:
        """Peak ground velocity log-standard deviation."""
        if self.INDEX_PGV is None:
            raise NotImplementedError
        else:
            return self._take(self._std_vector, self.INDEX_PGV, "pgv")

    @property
    def pgd(self) -> float | np.ndarray:
        """Peak ground displacement (PGD) computed by the model [cm]."""
        if self.INDEX_PGD is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGD, "pgd") * self.PGD_SCALE

    @property
    def ln_std_pgd(self) -> float | np.ndarray:
        """Peak ground displacement log-standard deviation."""
        if self.INDEX_PGD is None:
            raise NotImplementedError
        else:
            return self._take(self._std_vector, self.INDEX_PGD, "pgd")

    def _branch_list(self) -> list:
        """Branches as (weight, ln_resp, ln_std) of the computed coefficient rows."""
        if self._branches is None:
            return [(1.0, self._resp_vector, self._ln_std)]
        return list(self._branches)

    def ln_branches(self, im: str = "pga") -> list:
        r"""Branches of the logic tree of the ground motion distribution.

        Models with a logic tree of the median or the standard deviation, such as
        the USGS epistemic uncertainty branches of the nshmp-lib models, give the
        collapsed values in ``ln_pga`` etc. (the natural log of the weighted
        median, :math:`\ln \sum_i w_i \exp(\mu_i)`, as nshmp-lib's
        ``GroundMotions.combine``). The collapsed distribution is not the mixture
        of the branch distributions, so hazard calculations should instead sum the
        weighted exceedance probabilities of the branches, as nshmp-lib does
        (``ExceedanceModel.treeExceedanceCombined``). Models without branches
        return a single branch with a weight of 1.

        Parameters
        ----------
        im : str, optional
            intensity measure: "pga" (default), "pgv", "pgd", or "psa" for the
            computed spectral accelerations

        Returns
        -------
        branches : list of tuple
            (weight, ln_mean, ln_std) of each branch, where ln_mean and ln_std
            are the natural logs of the median (in the units of ``ln_pga`` etc.,
            i.e., without the PGV and PGD scale factors) and the standard
            deviations. Models whose branch weights depend on the period (e.g.,
            :class:`~pygmm.ground_motion.nga_east.NgaEast`) store them per computed
            period, and the weights of "psa" are then arrays with one value per
            period
        """
        if im == "psa":
            take = self._take_psa
        else:
            index = {
                "pga": self.INDEX_PGA,
                "pgv": self.INDEX_PGV,
                "pgd": self.INDEX_PGD,
            }.get(im, False)
            if index is False:
                raise ValueError(
                    f'im must be "pga", "pgv", "pgd", or "psa", not {im!r}'
                )
            if index is None:
                raise NotImplementedError

            def take(values):
                return self._take(values, index, im)

        return [
            (
                w if np.ndim(w) == 0 else take(np.asarray(w, dtype=float)),
                take(ln_resp),
                None if ln_std is None else take(ln_std),
            )
            for w, ln_resp, ln_std in self._branch_list()
        ]

    def _check_psa_computed(self) -> None:
        # Without ims, models keep their previous behavior (e.g., empty arrays
        # for models without spectral accelerations)
        if self._ims is not None and not len(self._psa_indices):
            raise ValueError(
                "Spectral accelerations were not computed. The model was created "
                f"with ims={list(self._ims)}; include 'psa_all', "
                "'psa_ngawest2_21', or a period such as 'psa_1p000' in ims."
            )

    def _take_psa(self, values: np.ndarray) -> np.ndarray:
        """Select the computed spectral accelerations from computed values."""
        self._check_psa_computed()
        return self._take(values, self._psa_indices, "psa")

    @property
    def _resp_vector(self) -> np.ndarray:
        """The computed response vector, or a clear error if there isn't one."""
        if self._ln_resp is None:
            raise RuntimeError(
                f"{type(self).__name__} has not computed a response vector; "
                "its __init__ must assign self._ln_resp."
            )
        return self._ln_resp

    @property
    def _std_vector(self) -> np.ndarray:
        """The computed standard deviations, or a clear error if absent."""
        if self._ln_std is None:
            raise RuntimeError(
                f"{type(self).__name__} has not computed standard deviations; "
                "its __init__ must assign self._ln_std."
            )
        return self._ln_std

    def _resp(self, index, im: str = "psa") -> np.ndarray:
        return np.exp(self._take(self._resp_vector, index, im))


class Parameter:
    """Model parameter.

    Parameters
    ----------
    name : str
        parameter name
    required : bool
        if the parameter is required
    default : None
        (optional) default value. Use *None* for no default value.

    """

    def __init__(self, name, required=False, default=None):
        """Initialize the parameter."""
        super().__init__()
        self._name = name
        self._required = required
        self._default = default

    def check(self, value):
        """Check the value against the limits."""
        if value is None and self.required:
            raise ValueError(self.name, "is a required parameter")

        if value is None:
            value = self.default
        return value

    @property
    def default(self):
        """Value to use as default."""
        return self._default

    @property
    def name(self):
        """Parameter name."""
        return self._name

    @property
    def required(self):
        """If the parameter is required."""
        return self._required


class NumericParameter(Parameter):
    """Numeric parameter.

    Parameters
    ----------
    name : str
        parameter name
    required : bool
        if the parameter is required
    default : float or int
        (optional) default value. Use *None* for no default value.

    """

    def __init__(
        self,
        name: str,
        required: bool = False,
        min_: float | None = None,
        max_: float | None = None,
        default: float | None = None,
    ):
        """Initialize parameter."""
        super().__init__(name, required, default)
        self._min = min_
        self._max = max_

    @property
    def min(self) -> float | None:
        """Minimum value, or None if unbounded below."""
        return self._min

    @property
    def max(self) -> float | None:
        """Maximum value, or None if unbounded above."""
        return self._max

    def check(self, value) -> float:
        """Check the value against the limits.

        The value can be a scalar or an array of values for a vectorized scenario.
        For an array, a single warning is issued for the values outside of the
        limits.
        """
        value = super().check(value)
        if value is not None and np.ndim(value) > 0:
            value = np.asarray(value)
            if self.min is not None and np.any(value < self.min):
                below = value < self.min
                warnings.warn(
                    f"{self.name} ({np.count_nonzero(below)} of {value.size} values, "
                    f"minimum of {np.min(value[below])}) "
                    f"is less than the recommended limit ({self.min}).",
                    UserWarning,
                    stacklevel=2,
                )
            if self.max is not None and np.any(self.max < value):
                above = self.max < value
                warnings.warn(
                    f"{self.name} ({np.count_nonzero(above)} of {value.size} values, "
                    f"maximum of {np.max(value[above])}) "
                    f"is greater than the recommended limit ({self.max}).",
                    UserWarning,
                    stacklevel=2,
                )
        elif value is not None:
            if self.min is not None and value < self.min:
                warnings.warn(
                    f"{self.name} ({value}) "
                    f"is less than the recommended limit ({self.min}).",
                    UserWarning,
                    stacklevel=2,
                )
            elif self.max is not None and self.max < value:
                warnings.warn(
                    f"{self.name} ({value}) "
                    f"is greater than the recommended limit ({self.max}).",
                    UserWarning,
                    stacklevel=2,
                )

        return value


class CategoricalParameter(Parameter):
    """Categorical parameter.

    Parameters
    ----------
    name : str
        parameter name
    required : bool
        if the parameter is required
    options : List[str]
        list of options
    default : str
        (optional) default option. Use *None* for no default value.

    """

    def __init__(
        self,
        name: str,
        required: bool = False,
        options: list[str] | None = None,
        default: str | None = None,
    ):
        """Initialize parameter."""
        super().__init__(name, required, default)
        self._options = options or []

    @property
    def options(self) -> list[str]:
        """Possible options."""
        return self._options

    def check(self, value) -> str:
        """Check the value against the limits.

        The value can be a single option or an array of options for a vectorized
        scenario. For an array, entries that are not one of the options are
        replaced with the default.
        """
        value = super().check(value)
        if np.ndim(value) > 0:
            value = np.asarray(value)
            # Comparing with each option is faster than np.isin for the short
            # option lists used by the models
            invalid = np.ones(value.shape, dtype=bool)
            for option in self.options:
                invalid &= ~equals(value, option)
            if np.any(invalid):
                warnings.warn(
                    f"{self.name} has {np.count_nonzero(invalid)} of {value.size} "
                    "values that are not one of the options. The following options "
                    "are possible: "
                    f"{', '.join([str(o) for o in self._options])}",
                    UserWarning,
                    stacklevel=2,
                )
                warnings.warn(
                    f"Using default value for {self.name}", UserWarning, stacklevel=2
                )
                value = np.where(invalid, self.default, value)
        elif value not in self.options:
            warnings.warn(
                f"{self.name} value of '{value}' "
                "is not one of the options. The following options are possible: "
                f"{', '.join([str(o) for o in self._options])}",
                UserWarning,
                stacklevel=2,
            )
            warnings.warn(
                f"Using default value for {self.name}", UserWarning, stacklevel=2
            )
            value = self.default

        return value


def symmetric_branches(
    ln_resp: ArrayLike, ln_std: ArrayLike, delta: ArrayLike, weights
) -> list:
    r"""Branches of a symmetric three-point logic tree of the median.

    The branches have medians of :math:`\mu - \delta`, :math:`\mu`, and
    :math:`\mu + \delta` (e.g., the USGS 5th, 50th, and 95th percentile
    epistemic branches of nshmp-lib, ``GroundMotions.createTree``), where
    :math:`\mu` is found from the collapsed median, :math:`\ln \sum_i w_i
    \exp(\mu_i)`, and all branches have the same standard deviation.

    Parameters
    ----------
    ln_resp : array_like
        natural log of the collapsed median
    ln_std : array_like
        standard deviation (natural log units)
    delta : array_like
        change of the natural log of the median of the lower and upper branches
    weights : sequence of float
        weights of the lower, central, and upper branches (the lower and upper
        weights are equal)

    Returns
    -------
    branches : list of tuple
        (weight, ln_resp, ln_std) of the lower, central, and upper branches
    """
    w_lo, w_mid, w_hi = weights
    ln_resp = np.asarray(ln_resp, dtype=float)
    delta = np.asarray(delta, dtype=float)
    ln_mid = ln_resp - np.log(w_lo * np.exp(-delta) + w_mid + w_hi * np.exp(delta))
    shape = np.broadcast_shapes(ln_mid.shape, np.shape(ln_std))
    ln_std = np.array(np.broadcast_to(ln_std, shape), dtype=float)
    return [
        (
            w,
            np.array(np.broadcast_to(ln_mid + sign * delta, shape), dtype=float),
            ln_std,
        )
        for w, sign in zip(weights, (-1.0, 0.0, 1.0))
    ]


def take_periods(values: np.ndarray, index) -> np.ndarray:
    """Select periods from a response or standard deviation array.

    A scalar scenario gives a 1-D array indexed by period. A vectorized scenario
    gives an array with periods along the last axis, e.g., shape (N, periods) for
    N scenarios.

    Parameters
    ----------
    values : array_like
        values with periods along the last axis
    index : int or array_like
        index or indices of the periods

    Returns
    -------
    values : :class:`np.ndarray`
        values at the selected periods
    """
    if np.ndim(values) <= 1:
        return values[index]
    return values[..., index]


def equals(values: npt.ArrayLike, option) -> np.ndarray:
    """Element-wise comparison of scenario values with an option.

    Equivalent to ``np.asarray(values) == option``. Arrays of fixed-width strings
    (e.g., mechanisms) are compared as integer code points, which is much faster
    than comparing strings for large vectorized scenarios.

    Parameters
    ----------
    values : array_like
        scalar or array of values
    option : str or other
        value to compare with

    Returns
    -------
    equal : :class:`np.ndarray`
        boolean array with the shape of `values`
    """
    values = np.asarray(values)
    if values.dtype.kind != "U" or values.ndim == 0 or not isinstance(option, str):
        return values == option
    width = values.dtype.itemsize // 4
    if len(option) > width:
        return np.zeros(values.shape, dtype=bool)
    codes = (
        np.ascontiguousarray(values).view(np.uint32).reshape(values.shape + (width,))
    )
    # Shorter strings are padded with null code points
    target = np.array([option], dtype=values.dtype).view(np.uint32)
    equal = codes[..., 0] == target[0]
    for j in range(1, width):
        equal &= codes[..., j] == target[j]
    return equal


def as_column(value: npt.ArrayLike) -> np.ndarray:
    """Add a trailing axis so scenario values broadcast against period coefficients.

    Coefficient arrays have one value per period. A scenario value with shape (N,)
    becomes shape (N, 1), so results have shape (N, periods). A scalar becomes shape
    (1,), which broadcasts to shape (periods,), so scalar scenarios keep their 1-D
    results.

    Parameters
    ----------
    value : array_like
        scalar or array of scenario values

    Returns
    -------
    value : :class:`np.ndarray`
        value with a trailing axis of length one
    """
    return np.asarray(value)[..., np.newaxis]


def as_scalar_or_array(value: npt.ArrayLike) -> float | np.ndarray:
    """Return a float for a scalar scenario and an array for a vectorized one.

    Parameters
    ----------
    value : array_like
        scalar or array of values

    Returns
    -------
    value : float or :class:`np.ndarray`
        a float if `value` has no dimensions, otherwise an array
    """
    if np.ndim(value) == 0:
        return float(np.asarray(value))
    return np.asarray(value)


def load_data_file(name, skip_header: int = 0) -> np.recarray:
    """Load a data file.

    Returns
    -------
    data : :class:`numpy.recarray`
       data values

    """
    fname = os.path.join(os.path.dirname(__file__), "data", name)
    return np.genfromtxt(
        fname, skip_header=skip_header, delimiter=",", names=True, case_sensitive=True
    ).view(np.recarray)


class Coefficients(collections.abc.Mapping):
    """Read-only container for model coefficients."""

    def __init__(self, **kwds):
        self._data = kwds

    def __getitem__(self, key):
        return self._data[key]

    def __getattr__(self, key):
        return self.__getitem__(key)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)
