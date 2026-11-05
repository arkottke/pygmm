r"""USGS nshmp-lib options of the NGA-West2 models.

The NGA-West2 models for active crustal regions,
:class:`~pygmm.ground_motion.abrahamson_silva_kamai_2014.AbrahamsonSilvaKamai2014`,
:class:`~pygmm.ground_motion.boore_stewart_seyhan_atkinson_2014.BooreStewartSeyhanAtkinson2014`,
:class:`~pygmm.ground_motion.campbell_bozorgnia_2014.CampbellBozorgnia2014`, and
:class:`~pygmm.ground_motion.chiou_youngs_2014.ChiouYoungs2014`, have options for
the variants of the models in the USGS nshmp-lib library (commit 44728a7d),
which is used for the current U.S. Geological Survey national seismic hazard
models (NSHMs). The options are off by default, which gives the published
models. ``GMM_IDS`` of each model maps the nshmp-lib ``Gmm`` ids (e.g.,
``ASK_14_BASIN``) to the options and scenario values:

- ``nshmp``: the nshmp-lib implementation of the model (``AbrahamsonEtAl_2014``,
  ``BooreEtAl_2014``, ``CampbellBozorgnia_2014``, and ``ChiouYoungs_2014``),
  which differs from the published model in details such as the treatment of
  missing basin depths and the hanging-wall conditions (see the models). In
  all of the models, a ``depth_1_0`` or ``depth_2_5`` of *None* or NaN gives
  the basin term of the :math:`V_{S30}` based reference depth (no basin term
  in ASK14, BSSA14, and CY14), and the unspecified style of faulting
  (``mechanism`` of "U", see :func:`mechanism_from_rake`) has no style of
  faulting term (as strike-slip) in ASK14, CB14, and CY14.
- ``epistemic``: the USGS additional epistemic uncertainty of the NGA-West2
  models (``NgaUncertainty.NGA_WEST2``) is represented by three branches with
  means :math:`\mu - \epsilon`, :math:`\mu`, and :math:`\mu + \epsilon`,
  where :math:`\epsilon` depends on the magnitude and the Joyner-Boore
  distance (see :func:`epistemic_epsilon`), and weights of 0.185, 0.63, and
  0.185. As in nshmp-lib (``GroundMotions.combine``), the branches are
  collapsed to a single median, :math:`\ln \sum_i w_i \exp(\mu_i)`, and the
  standard deviation, which is the same for all branches. The collapsed
  values are not suitable for hazard calculations, which should sum the
  weighted exceedance probabilities of the branches, as nshmp-lib does
  (``ExceedanceModel.treeExceedanceCombined``); the branches are given by
  :meth:`~pygmm.model.GroundMotionModel.ln_branches`.
- ``basin``: the USGS deep basin model, which keeps the basin (sediment
  depth) term of the model only for spectral periods longer than 0.5 s and
  scales it from 0 at :math:`Z_{1.0}` of 0.3 km to the full term at 0.5 km
  (ASK14, BSSA14, and CY14, see :func:`basin_scale`), or, for CB14, from the
  :math:`V_{S30}` based (reference depth) term at :math:`Z_{2.5}` of 1 km to
  the full term at 3 km. The scaled term is multiplied by 0.585 at 0.75 s.
  This is the variant used by the current CONUS NSHM (nshm-conus 6.2.0, the
  ``*_14_BASIN`` ids with weights of 0.25).
- ``cybershake`` (requires ``basin``): the 2023 NSHM CyberShake based basin
  terms for spectral periods longer than 1.9 s (alternative basin
  coefficients and a constant of 0.1 added to the basin term), which are used
  in the Los Angeles region.
- ``prvi``: the 2025 Puerto Rico and U.S. Virgin Islands adjustments of the
  constant and linear site coefficients.

The ``GMM_IDS`` of ASK14 and CY14 also have the ``vs_source`` scenario value,
which is "inferred" except for the ``*_14_VS30_MEASURED`` variants (nshmp-lib
uses the aleatory variability for an inferred :math:`V_{S30}`), so the values
are split into the model options (the arguments of the model) and the scenario
values. The ``epistemic``, ``basin``, ``cybershake``, and ``prvi`` options can
also be used without ``nshmp``, with the published model. BSSA14 uses the
published coefficients, which differ slightly from the rounded nshmp-lib
coefficients (see
:class:`~pygmm.ground_motion.boore_stewart_seyhan_atkinson_2014.BooreStewartSeyhanAtkinson2014`).

The CyberShake and PRVI coefficients are only given by nshmp-lib for PGA, PGV,
and the 21 NGA-West2 comparison periods. The results at the other periods of
the models are NaN with these options (e.g., 6 s of ASK14 with ``prvi``), so use
``ims="psa_ngawest2_21"`` for the spectral accelerations.

References
----------
USGS nshmp-lib, ``gov.usgs.earthquake.nshmp.gmm.AbrahamsonEtAl_2014``,
``BooreEtAl_2014``, ``CampbellBozorgnia_2014``, ``ChiouYoungs_2014``,
``GmmUtils``, ``NgaUncertainty``, and ``GroundMotions``, commit 44728a7d
(https://code.usgs.gov/ghsc/nshmp/nshmp-lib).
"""

from typing import Any

import numpy as np

from .. import model
from .._types import ArrayLike

#: Weights of the epistemic branches (mu - epsilon, mu, mu + epsilon)
EPI_WTS = (0.185, 0.63, 0.185)

#: NGA-West2 additional epistemic uncertainty (nshmp-lib
#: ``NgaUncertainty.NGA_WEST2``), in natural log units. Rows are for Joyner-Boore
#: distances less than 10 km, from 10 to 30 km, and 30 km and greater, and the
#: columns are for magnitudes less than 6, from 6 to 7, and 7 and greater.
EPI_NGA_WEST2 = np.array(
    [
        [0.37, 0.25, 0.40],
        [0.22, 0.23, 0.36],
        [0.22, 0.23, 0.33],
    ]
)

#: :math:`Z_{1.0}` (km) of the USGS deep basin model with a scale of 0 and 1
BASIN_Z1P0 = (0.3, 0.5)
#: :math:`Z_{2.5}` (km) of the USGS deep basin model with a scale of 0 and 1
BASIN_Z2P5 = (1.0, 3.0)
#: Spectral periods (s) longer than this have the CyberShake basin terms
CYBERSHAKE_MIN_PERIOD = 1.9
#: Constant added to the CyberShake basin terms
CYBERSHAKE_CONSTANT = 0.1


def epistemic_epsilon(mag: ArrayLike, dist_jb: ArrayLike) -> np.ndarray:
    r"""USGS NGA-West2 additional epistemic uncertainty.

    The :math:`\epsilon` of the epistemic branches with means :math:`\mu -
    \epsilon`, :math:`\mu`, and :math:`\mu + \epsilon` (nshmp-lib
    ``NgaUncertainty.NGA_WEST2``; see :data:`EPI_NGA_WEST2`).

    Parameters
    ----------
    mag : float or array_like
        moment magnitude
    dist_jb : float or array_like
        Joyner-Boore distance [km]

    Returns
    -------
    epsilon : :class:`np.ndarray`
        epistemic uncertainty in natural log units

    Examples
    --------
    >>> from pygmm.ground_motion.nshmp import epistemic_epsilon
    >>> epistemic_epsilon([5.5, 6.5, 7.5], [5.0, 20.0, 50.0]).tolist()
    [0.37, 0.23, 0.33]
    """
    mag = np.asarray(mag, dtype=float)
    dist_jb = np.asarray(dist_jb, dtype=float)
    mi = np.where(mag < 6, 0, np.where(mag < 7, 1, 2))
    ri = np.where(dist_jb < 10, 0, np.where(dist_jb < 30, 1, 2))
    return EPI_NGA_WEST2[ri, mi]


def epistemic_ln_factor(mag: ArrayLike, dist_jb: ArrayLike) -> np.ndarray:
    r"""Change of the natural log of the median from the epistemic branches.

    :math:`\ln(0.185 e^{-\epsilon} + 0.63 + 0.185 e^{\epsilon})`, which is
    added to the natural log of the median of the central branch to give the
    collapsed median, :math:`\ln \sum_i w_i \exp(\mu_i)`.

    Parameters
    ----------
    mag : float or array_like
        moment magnitude
    dist_jb : float or array_like
        Joyner-Boore distance [km]

    Returns
    -------
    ln_factor : :class:`np.ndarray`
        change of the natural log of the median
    """
    eps = epistemic_epsilon(mag, dist_jb)
    w_lo, w_mid, w_hi = EPI_WTS
    return np.log(w_lo * np.exp(-eps) + w_mid + w_hi * np.exp(eps))


def epistemic_branches(
    ln_resp: np.ndarray, ln_std: np.ndarray, mag: ArrayLike, dist_jb: ArrayLike
) -> tuple[np.ndarray, list]:
    """Collapsed median and branches of the USGS epistemic uncertainty.

    nshmp-lib ``GroundMotions.createNgaTree`` and ``combine``.

    Parameters
    ----------
    ln_resp : :class:`np.ndarray`
        natural log of the median (of the central branch), with periods along
        the last axis
    ln_std : :class:`np.ndarray`
        standard deviation (natural log units)
    mag : array_like
        moment magnitude, with a trailing axis for the periods (see
        :func:`pygmm.model.as_column`)
    dist_jb : array_like
        Joyner-Boore distance [km], with a trailing axis for the periods

    Returns
    -------
    ln_resp : :class:`np.ndarray`
        natural log of the collapsed median
    branches : list of tuple
        (weight, ln_resp, ln_std) of the branches
    """
    eps = epistemic_epsilon(mag, dist_jb)
    shape = np.broadcast_shapes(np.shape(ln_resp), np.shape(ln_std), eps.shape)
    ln_std_b = np.array(np.broadcast_to(ln_std, shape), dtype=float)
    branches = [
        (
            w,
            np.array(np.broadcast_to(ln_resp + sign * eps, shape), dtype=float),
            ln_std_b,
        )
        for w, sign in zip(EPI_WTS, (-1.0, 0.0, 1.0))
    ]
    return ln_resp + epistemic_ln_factor(mag, dist_jb), branches


def broadcast(ln_resp: np.ndarray, ln_std: np.ndarray) -> tuple:
    """Broadcast the median and standard deviation to the same shape."""
    shape = np.broadcast_shapes(np.shape(ln_resp), np.shape(ln_std))
    return (
        np.array(np.broadcast_to(ln_resp, shape), dtype=float),
        np.array(np.broadcast_to(ln_std, shape), dtype=float),
    )


def mechanism_from_rake(rake: ArrayLike) -> np.ndarray:
    """Style of faulting from the rake, as in nshmp-lib.

    nshmp-lib (``GmmUtils.rakeToFaultStyle_NSHMP``) divides the styles on 45
    degree diagonals: reverse ("RS") for rakes from 45 to 135 degrees, normal
    ("NS") for rakes from -135 to -45 degrees, strike-slip ("SS") otherwise,
    and unspecified ("U") for a rake of NaN.

    Parameters
    ----------
    rake : float or array_like
        rake angle [deg]

    Returns
    -------
    mechanism : :class:`np.ndarray`
        mechanism ("SS", "NS", "RS", or "U")

    Examples
    --------
    >>> from pygmm.ground_motion.nshmp import mechanism_from_rake
    >>> mechanism_from_rake([0.0, 90.0, -90.0, 180.0, float("nan")]).tolist()
    ['SS', 'RS', 'NS', 'SS', 'U']
    """
    rake = np.asarray(rake, dtype=float)
    return np.select(
        [np.isnan(rake), (rake >= 45) & (rake <= 135), (rake >= -135) & (rake <= -45)],
        ["U", "RS", "NS"],
        "SS",
    )


def basin_scale(
    periods: ArrayLike, depth: ArrayLike, upper: float, lower: float
) -> np.ndarray:
    """USGS deep basin scaling (``GmmUtils.deltaZ1scale`` and ``deltaZ25scale``).

    The scale is 0 for periods of 0.5 s and shorter (and PGA and PGV), depths
    of `upper` or less, or NaN depths, increases linearly to 1 at `lower`, and
    is multiplied by 0.585 at 0.75 s.

    Parameters
    ----------
    periods : array_like
        spectral periods [s], with 0 for PGA and -1 for PGV
    depth : float or array_like
        basin depth [km]
    upper : float
        depth with a scale of 0 [km]
    lower : float
        depth with a scale of 1 [km]

    Returns
    -------
    scale : :class:`np.ndarray`
        basin scale

    Examples
    --------
    >>> from pygmm.ground_motion.nshmp import basin_scale
    >>> basin_scale([0.5, 0.75, 1.0], 0.4, 0.3, 0.5).round(4).tolist()
    [0.0, 0.2925, 0.5]
    """
    periods = np.asarray(periods)
    depth = np.asarray(depth, dtype=float)
    scale = (np.clip(depth, upper, lower) - upper) / (lower - upper)
    scale = np.where(periods == 0.75, scale * 0.585, scale)
    # NaN depths compare False
    return np.where((periods > 0.5) & (depth > upper), scale, 0.0)


def load_coeffs(name: str, periods: np.ndarray) -> Any:
    """Load nshmp-lib coefficients for the periods of a model.

    The rows of the data file are matched with `periods` (with 0 for PGA and
    -1 for PGV), and the periods that are not in the file have NaN values.

    Parameters
    ----------
    name : str
        name of the data file, which has a ``period`` column
    periods : :class:`np.ndarray`
        periods of the model coefficients

    Returns
    -------
    coeff : :class:`numpy.recarray`
        coefficients with one row per period
    """
    data = model.load_data_file(name, 1)
    names = [n for n in data.dtype.names or () if n != "period"]
    coeff = np.recarray(
        len(periods), dtype=[("period", float)] + [(n, float) for n in names]
    )
    coeff["period"] = periods
    found = np.isclose(
        np.asarray(periods)[:, np.newaxis], data["period"], rtol=1e-9, atol=0
    )
    assert np.all(found.sum(axis=0) == 1), f"periods of {name} are not in the model"
    row = np.where(found.any(axis=1), found.argmax(axis=1), -1)
    for n in names:
        coeff[n] = np.where(row >= 0, data[n][row], np.nan)
    return coeff


def replace_coeffs(coeff: Any, other: Any) -> Any:
    """Coefficients with the columns of `other` (e.g., the PRVI coefficients).

    Parameters
    ----------
    coeff : :class:`numpy.recarray`
        model coefficients
    other : :class:`numpy.recarray`
        replacement coefficients with the same rows (see :func:`load_coeffs`)

    Returns
    -------
    coeff : :class:`numpy.recarray`
        copy of `coeff` with the columns of `other`
    """
    assert np.array_equal(coeff["period"], other["period"])
    coeff = coeff.copy()
    for name in other.dtype.names or ():
        coeff[name] = other[name]
    return coeff


class NshmpOptions:
    """Model options of the nshmp-lib variants of the NGA-West2 models.

    See :mod:`pygmm.ground_motion.nshmp`.
    """

    #: Model options and scenario values of the nshmp-lib ``Gmm`` ids
    GMM_IDS: dict[str, dict[str, Any]] = {}
    #: nshmp-lib PRVI coefficients
    COEFF_PRVI: Any = None
    #: nshmp-lib CyberShake basin coefficients
    COEFF_CYBERSHAKE: Any = None

    _nshmp = False
    _epistemic = False
    _basin = False
    _cybershake = False
    _prvi = False

    def _set_options(
        self,
        nshmp: bool,
        epistemic: bool,
        basin: bool,
        cybershake: bool,
        prvi: bool,
    ) -> None:
        """Check and set the options (before the scenario is checked)."""
        if cybershake and not basin:
            raise ValueError(
                "cybershake=True requires basin=True (nshmp-lib *_14_CYBERSHAKE)"
            )
        self._nshmp = bool(nshmp)
        self._epistemic = bool(epistemic)
        self._basin = bool(basin)
        self._cybershake = bool(cybershake)
        self._prvi = bool(prvi)

    @property
    def nshmp(self) -> bool:
        """If the nshmp-lib implementation of the model is used."""
        return self._nshmp

    @property
    def epistemic(self) -> bool:
        """If the USGS epistemic branches are included."""
        return self._epistemic

    @property
    def basin(self) -> bool:
        """If the USGS deep basin model is used."""
        return self._basin

    @property
    def cybershake(self) -> bool:
        """If the CyberShake basin terms are used."""
        return self._cybershake

    @property
    def prvi(self) -> bool:
        """If the PRVI coefficients are used."""
        return self._prvi

    def _options_coeff(self, coeff: Any) -> Any:
        """Model coefficients, with the PRVI coefficients if `prvi`."""
        if self._prvi:
            return replace_coeffs(coeff, self.COEFF_PRVI)
        return coeff

    def _cybershake_rows(self, periods: np.ndarray) -> np.ndarray:
        """Periods with the CyberShake basin terms (``GmmUtils.cybershakeImt``)."""
        return self._cybershake & (np.asarray(periods) > CYBERSHAKE_MIN_PERIOD)
