# -*- coding: utf-8 -*-
"""Basic models."""
import collections
import logging
import os
from typing import List
from typing import Optional

import numpy as np
from scipy.interpolate import interp1d

from .types import ArrayLike


class Scenario(collections.UserDict):
    r"""An eathquake scenario used in all ground motion models.

    Parameters
    ----------
    depth_1_0 : float
        depth to the 1.0 km∕s shear-wave velocity horizon beneath the site,
        :math:`Z_{1.0}` in (km).
    depth_2_5 : float
        depth to the 2.5 km∕s shear-wave velocity horizon beneath the site,
        :math:`Z_{2.5}` in (km).
    depth_tor : float
        depth to the top of the rupture plane (:math:`Z_{tor}`, km).
    depth_bor : float
        depth to the bottom of the rupture plane (:math:`Z_{bor}`, km).
    depth_bot : float
        depth to bottom of seismogenic crust (km).
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
        Spudich (2014, :cite:`spudich14`)) centered on the earthquake-specific
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
        :ref:`Mechanism` for more information.
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
        self._check_keys(kwds.keys())
        other = self.copy()
        other.update(**kwds)
        return other

    def _check_keys(self, keys):
        for k in keys:
            if k not in self.KNOWN_KEYS:
                raise Warning("%s is not a recognized scenario key!" % k)


class Model(object):
    #: Long name of the model
    NAME = ""
    #: Short name of the model
    ABBREV = ""

    def __init__(self, *args, **kwargs):
        """Initialize the model."""
        super(Model, self).__init__()

        self._ln_resp = None
        self._ln_std = None

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
        return self._scenario


class GroundMotionModel(Model):
    """Abstract class for ground motion prediction models."""

    #: Indices for the spectral accelerations
    INDICES_PSA = np.array([])
    #: Indices of the periods
    PERIODS = np.array([])
    #: Index of the peak ground acceleration
    INDEX_PGA = None
    #: Index of the peak ground velocity
    INDEX_PGV = None
    #: Index of the peak ground displacement
    INDEX_PGD = None
    #: Limits of model applicability
    LIMITS = dict()
    #: Model parameters
    PARAMS = []
    #: Scale factor to apply to get PGV in cm/sec
    PGV_SCALE = 1.0
    #: Scale factor to apply to get PGD in cm
    PGD_SCALE = 1.0

    def __init__(self, scenario: Scenario):
        """Initialize the model."""
        super(GroundMotionModel, self).__init__(scenario)

        self._ln_resp = None
        self._ln_std = None

        # Select the used parameters and check them against the recommended
        # values
        self._scenario = Scenario(
            **{p.name: scenario.get(p.name, None) for p in self.PARAMS}
        )
        self._check_inputs()

    def interp_ln_spec_accels(
        self, periods: ArrayLike, kind: Optional[str] = "linear"
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
            self._ln_resp[self.INDICES_PSA],
            kind=kind,
            copy=False,
            bounds_error=False,
            fill_value=np.nan,
        )(np.log(periods))

    def interp_spec_accels(
        self, periods: ArrayLike, kind: Optional[str] = "linear"
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
        self, periods: ArrayLike, kind: Optional[str] = "linear"
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
                self._ln_std[self.INDICES_PSA],
                kind=kind,
                copy=False,
                bounds_error=False,
                fill_value=np.nan,
            )(np.log(periods))

    @property
    def periods(self) -> np.ndarray:
        """Periods specified by the model."""
        return self.PERIODS[self.INDICES_PSA]

    @property
    def spec_accels(self) -> np.ndarray:
        """Pseudo-spectral accelerations computed by the model (g)."""
        return self._resp(self.INDICES_PSA)

    @property
    def ln_stds(self) -> np.ndarray:
        """Pseudo-spectral accelerations log-standard deviation."""
        if self._ln_std is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDICES_PSA]

    @property
    def pga(self) -> float:
        """Peak ground acceleration (PGA) computed by the model (g)."""
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGA)

    @property
    def ln_std_pga(self) -> float:
        """Peak ground accelaration log-standard deviation."""
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDEX_PGA]

    @property
    def pgv(self) -> float:
        """Peak ground velocity (PGV) computed by the model (cm/sec)."""
        if self.INDEX_PGV is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGV) * self.PGV_SCALE

    @property
    def ln_std_pgv(self) -> float:
        """Peak ground velocity log-standard deviation."""
        if self.INDEX_PGV is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDEX_PGV]

    @property
    def pgd(self) -> float:
        """Peak ground displacement (PGD) computed by the model (cm)."""
        if self.INDEX_PGD is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGD) * self.PGD_SCALE

    @property
    def ln_std_pgd(self) -> float:
        """Peak ground displacement log-standard deviation."""
        if self.INDEX_PGD is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDEX_PGD]

    def _resp(self, index) -> np.ndarray:
        if index is not None:
            return np.exp(self._ln_resp[index])

    def _check_inputs(self) -> None:
        for p in self.PARAMS:
            self._scenario[p.name] = p.check(self._scenario[p.name])


class Parameter(object):
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
        super(Parameter, self).__init__()
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
        min_: Optional[float] = None,
        max_: Optional[float] = None,
        default: Optional[float] = None,
    ):
        """Initialize parameter."""
        super(NumericParameter, self).__init__(name, required, default)
        self._min = min_
        self._max = max_

    @property
    def min(self) -> float:
        """Minimum value."""
        return self._min

    @property
    def max(self) -> float:
        """Maximum value."""
        return self._max

    def check(self, value) -> float:
        """Check the value against the limits."""
        value = super(NumericParameter, self).check(value)
        if value is not None:
            if self.min is not None and value < self.min:
                logging.warning(
                    "%s (%g) is less than the recommended limit (%g).",
                    self.name,
                    value,
                    self.min,
                )
            elif self.max is not None and self.max < value:
                logging.warning(
                    "%s (%g) is greater than the recommended limit (%g).",
                    self.name,
                    value,
                    self.max,
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
    options : list[str]
        list of options
    default : str
        (optional) default option. Use *None* for no default value.

    """

    def __init__(
        self,
        name: str,
        required: bool = False,
        options: Optional[List[str]] = None,
        default: Optional[str] = None,
    ):
        """Initialize parameter."""
        super(CategoricalParameter, self).__init__(name, required, default)
        self._options = options or []

    @property
    def options(self) -> List[str]:
        """Possible options."""
        return self._options

    def check(self, value) -> str:
        """Check the value against the limits."""
        value = super(CategoricalParameter, self).check(value)
        if value not in self.options:
            alert = logging.error if self.required else logging.warning
            alert(
                '%s value of "%s" is not one of the options. The following'
                " options are possible: %s",
                self.name,
                value,
                ", ".join([str(o) for o in self._options]),
            )

            if not self.required:
                logging.warning("Using default value for %s", self.name)
                value = self.default

        return value


def load_data_file(name, skip_header=None) -> np.recarray:
    """Load a data file.

    Returns
    -------
    data : :class:`numpy.recarray`
       data values

    """
    fname = os.path.join(os.path.dirname(__file__), "data", name)
    return np.recfromcsv(fname, skip_header=skip_header, case_sensitive=True).view(
        np.recarray
    )
