#!/usr/bin/python
# -*- coding: utf-8 -*-

from __future__ import division

import collections
import logging
import os

import numpy as np
from scipy.interpolate import interp1d

from .types import ArrayLike


<<<<<<< HEAD
class Scenario(collections.UserDict):
    """An eathquake scenario used in all ground motion models."""

    KNOWN_KEYS = ['depth_1_0', 'depth_2_5', 'depth_tor', 'depth_bor',
                  'depth_bot', 'depth_hyp', 'dip', 'dist_jb', 'dist_epi',
                  'dist_hyp', 'dist_rup', 'dist_x', 'dist_y0', 'dpp_centered',
                  'mag', 'mechanism', 'on_hanging_wall', 'region', 'v_s30',
                  'vs_source', 'width']

    def __init__(self, **kwds):
        """Initialize the scenario.

        Keyword Args:
            depth_1_0 (float): depth to the 1.0 km∕s shear-wave velocity
                horizon beneath the site, :math:`Z_{1.0}` in (km).
            depth_2_5 (float): depth to the 2.5 km∕s shear-wave velocity
                horizon beneath the site, :math:`Z_{2.5}` in (km).
            depth_tor (float): depth to the top of the rupture plane
                (:math:`Z_{tor}`, km).
            depth_bor (float): depth to the bottom of the rupture plane
                (:math:`Z_{bor}`, km).
            depth_bot (float): depth to bottom of seismogenic crust (km).
            dip (float): fault dip angle (:math:`\phi`, deg).
            dist_jb (float): Joyner-Boore distance to the rupture plane
                (:math:`R_\\text{JB}`, km)
            dist_epi (float): Epicentral distance to the rupture plane
                (:math:`R_\\text{epi}`, km)
            dist_hyp (float): Hypocentral distance to the rupture plane
                (:math:`R_\\text{hyp}`, km).
            dist_rup (float): closest distance to the rupture plane
                (:math:`R_\\text{rup}`, km)
            dist_x (float): site coordinate measured perpendicular to the
                fault strike from the fault line with the down-dip direction
                being positive (:math:`R_x`, km).
            dist_y0 (float): the horizontal distance off the end of the
                rupture measured parallel to strike (:math:`R_{y0}`, km).
            dpp_centered (float): direct point parameter (DPP) for directivity
                effect (see Chiou and Spudich (2014, :cite:`spudich14`))
                centered on the earthquake-specific average DPP for
                California.
            mag (float): moment magnitude of the event (:math:`M_w`)
            mechanism (str): fault mechanism. Valid options: "SS", "NS", "RS",
                and "U". See :ref:`Mechanism` for more information.
            on_hanging_wall (bool): If the site is located on the hanging wall
                of the fault. If *None*, then *False* is assumed.
            region (str): region. Valid options are specified FIXME.
            v_s30 (float): time-averaged shear-wave velocity over the top 30 m
                of the site (:math:`V_{s30}`, m/s).
            vs_source (str): source of the `v_s30` value.  Valid options:
                "measured", "inferred"
            width (float): Down-dip width of the fault.
        """
        super(Scenario, self).__init__(kwds)
        for k in self.data:
            if k not in self.KNOWN_KEYS:
                raise Warning(f'{k} is not a recognized scenario key!')

    def __getitem__(self, key):
        return self.data[key]



class Model(object):
    """Abstract class for ground motion prediction models.
    """

    #: Long name of the model
    NAME = ''
    #: Short name of the model
    ABBREV = ''
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
    PGV_SCALE = 1.
    #: Scale factor to apply to get PGD in cm
    PGD_SCALE = 1.

    def __init__(self, scenario):
        """Initialize the model.
        """
=======
class Model(object):
    """Abstract class for ground motion prediction models.

    A common set of keywords is used for all ground motion models. Depending
    the model these keywords may be optional, required, or not used.

    Parameters
    ----------
    depth_1_0 : float
        depth to the 1.0 km∕s shear-wave velocity
        horizon beneath the site, :math:`Z_{1.0}` in (km).
    depth_2_5 : float
        depth to the 2.5 km∕s shear-wave velocity
        horizon beneath the site, :math:`Z_{2.5}` in (km).
    depth_tor : float
        depth to the top of the rupture plane
        (:math:`Z_{tor}`, km).
    depth_bor : float
        depth to the bottom of the rupture plane
        (:math:`Z_{bor}`, km).
    depth_bot : float
        depth to bottom of seismogenic crust (km).
    dip : float
        fault dip angle (:math:`\phi`, deg).
    dist_jb : float
        Joyner-Boore distance to the rupture plane
        (:math:`R_\\text{JB}`, km)
    dist_epi : float
        Epicentral distance to the rupture plane
        (:math:`R_\\text{epi}`, km)
    dist_hyp : float
        Hypocentral distance to the rupture plane
        (:math:`R_\\text{hyp}`, km).
    dist_rup : float
        closest distance to the rupture plane
        (:math:`R_\\text{rup}`, km)
    dist_x : float
        site coordinate measured perpendicular to the
        fault strike from the fault line with the down-dip direction
        being positive (:math:`R_x`, km).
    dist_y0 : float
        the horizontal distance off the end of the
        rupture measured parallel to strike (:math:`R_{y0}`, km).
    dpp_centered : float
        direct point parameter (DPP) for directivity
        effect (see Chiou and Spudich (2014, :cite:`spudich14`))
        centered on the earthquake-specific average DPP for
        California.
    mag : float
        moment magnitude of the event (:math:`M_w`)
    mechanism : str
        fault mechanism. Valid options: "SS", "NS", "RS",
        and "U". See :ref:`Mechanism` for more information.
    on_hanging_wall : bool
        If the site is located on the hanging wall
        of the fault. If *None*, then *False* is assumed.
    region : str
        region. Valid options are specified FIXME.
    v_s30 : float
        time-averaged shear-wave velocity over the top 30 m
        of the site (:math:`V_{s30}`, m/s).
    vs_source : str
        source of the `v_s30` value.  Valid options:
        "measured", "inferred"
    width : float
        Down-dip width of the fault.

    Attributes
    ----------
    NAME : str
        Long name of the model
    ABBREV : str
        Short name of the model
    INDICES_PSA : :class:`np.ndarray`
        Indices for the spectral accelerations
    PERIODS : :class:`np.ndarray`
        Periods of the spectral accelerations
    INDEX_PGA : int
        Index of the peak ground acceleration
    INDEX_PGV : int or None
        Index of the peak ground velocity. None if not provided by the model.
    INDEX_PGD : int or None
        Index of the peak ground displacement. None if not provided by the
        model.
    LIMITS : dict
        Limits of the model parameters. Used for checking the input.
    PARAMS : list[str]
        List of model parameters. For possible names see `Parameters`_
    PGV_SCALE : float
        Scale factor to apply to get PGV in cm/sec.
    PGD_SCALE : float
        Scale factor to apply to get PGD in cm
    """

    NAME = ''
    ABBREV = ''
    INDICES_PSA = np.array([])
    PERIODS = np.array([])
    INDEX_PGA = None
    INDEX_PGV = None
    INDEX_PGD = None
    LIMITS = dict()
    PARAMS = []
    PGV_SCALE = 1.
    PGD_SCALE = 1.

    def __init__(self, **kwds):
        """Initialize the model."""
>>>>>>> 463f156a57779d7fb9def11b795e00bc38ad0dd8
        super(Model, self).__init__()

        self._ln_resp = None
        self._ln_std = None

        # Select the used parameters and check them against the recommended
        # values
        self._scenario = {p.name: scenario.get(p.name, None) for p in self.PARAMS}
        self._check_inputs()

    def interp_spec_accels(self, periods: ArrayLike,
                           kind: str='linear') -> np.ndarray:
        """Return the pseudo-spectral acceleration at the provided damping
        at specified periods.

        Interpolation is done in natural log space.

        Parameters
        ----------
        periods : array_like
            Spectral periods to interpolate the response.
        kind : str, optional
            See :func:`scipy.interpolate.interp1d` for description of kind.
            Options include: 'linear' (default), 'nearest', 'zero', 'slinear',
            'quadratic', and 'cubic'

        Returns
        -------
        spec_accels : np.ndarray
            Interpolated spectral accelerations
        """
        return np.exp(
            interp1d(
                np.log(self.periods),
                self._ln_resp[self.INDICES_PSA],
                kind=kind,
                copy=False,
                bounds_error=False,
                fill_value=np.nan, )(np.log(periods)))

    def interp_ln_stds(self, periods: ArrayLike,
                       kind: str='linear') -> np.ndarray:
        """Return the logarithmic standard deviation
        (:math:`\\sigma_{ \\ln}`) of spectral acceleration at the provided
        damping at specified periods.

        Parameters
        ----------
        periods : array_like
            Spectral periods to interpolate the response.
        kind : str, optional
            See :func:`scipy.interpolate.interp1d` for description of kind.
            Options include: 'linear' (default), 'nearest', 'zero', 'slinear',
            'quadratic', and 'cubic'

        Returns
        -------
        ln_stds : np.ndarray
            Interpolated logarithmic standard deviations
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
                fill_value=np.nan, )(np.log(periods))

    @property
    def periods(self):
        """Periods specified by the model."""
        return self.PERIODS[self.INDICES_PSA]

    @property
    def spec_accels(self):
        """Pseudo-spectral accelerations computed by the model (g)."""
        return self._resp(self.INDICES_PSA)

    @property
    def ln_stds(self):
        """Logarithmic standard deviation of the pseudo-spectral
        accelerations.
        """
        if self._ln_std is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDICES_PSA]

    @property
    def pga(self):
        """Peak ground acceleration (PGA) computed by the model (g)."""
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGA)

    @property
    def ln_std_pga(self):
        """Logarithmic standard deviation (:math:`\\sigma_{ \\ln}`) of the
        peak ground acceleration computed by the model.
        """
        if self.INDEX_PGA is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDEX_PGA]

    @property
    def pgv(self):
        """Peak ground velocity (PGV) computed by the model (cm/sec)."""
        if self.INDEX_PGV is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGV) * self.PGV_SCALE

    @property
    def ln_std_pgv(self):
        """Logarithmic standard deviation (:math:`\\sigma_{ \\ln}`) of the
        peak ground velocity computed by the model.
        """
        if self.INDEX_PGV is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDEX_PGV]

    @property
    def pgd(self):
        """Peak ground displacement (PGD) computed by the model (cm)."""
        if self.INDEX_PGD is None:
            raise NotImplementedError
        else:
            return self._resp(self.INDEX_PGD) * self.PGD_SCALE

    @property
    def ln_std_pgd(self):
        """Logarithmic standard deviation (:math:`\\sigma_{ \\ln}`) of the
        peak ground displacement computed by the model.
        """
        if self.INDEX_PGD is None:
            raise NotImplementedError
        else:
            return self._ln_std[self.INDEX_PGD]

    def _resp(self, index):
        if index is not None:
            return np.exp(self._ln_resp[index])

    def _check_inputs(self):
        for p in self.PARAMS:
            self._scenario[p.name] = p.check(self._scenario[p.name])


class Parameter(object):
    def __init__(self, name, required=False, default=None):
        super(Parameter, self).__init__()

        self._name = name
        self._required = required
        self._default = default

    def check(self, value):
        raise NotImplementedError

    @property
    def default(self):
        return self._default

    @property
    def name(self):
        return self._name

    @property
    def required(self):
        return self._required


class NumericParameter(Parameter):
    def __init__(self,
                 name,
                 required=False,
                 min_=None,
                 max_=None,
                 default=None):
        super(NumericParameter, self).__init__(name, required, default)
        self._min = min_
        self._max = max_

    @property
    def min(self):
        return self._min

    @property
    def max(self):
        return self._max

    def check(self, value):
        if value is None and self.required:
            raise ValueError(self.name, 'is a required parameter')

        if value is None:
            value = self.default
        else:
            if self.min is not None and value < self.min:
                logging.warning(
                    '%s (%g) is less than the recommended limit (%g).',
                    self.name, value, self.min)
            elif self.max is not None and self.max < value:
                logging.warning(
                    '%s (%g) is greater than the recommended limit (%g).',
                    self.name, value, self.max)

        return value


class CategoricalParameter(Parameter):
    def __init__(self, name, required=False, options=None, default=''):
        super(CategoricalParameter, self).__init__(name, required, default)
        self._options = options or []

    @property
    def options(self):
        return self._options

    def check(self, value):
        if value is None and self.required:
            raise ValueError(self.name, 'is a required parameter')

        if value is None:
            value = self.default
        elif value not in self.options:
            alert = logging.error if self._required else logging.warning
            alert('%s value of "%s" is not one of the options. The following'
                  ' options are possible: %s', self._name, value,
                  ', '.join([str(o) for o in self._options]))

        return value


def load_data_file(name, skip_header=None):
    fname = os.path.join(os.path.dirname(__file__), 'data', name)
    return np.recfromcsv(
        fname, skip_header=skip_header, case_sensitive=True).view(np.recarray)
