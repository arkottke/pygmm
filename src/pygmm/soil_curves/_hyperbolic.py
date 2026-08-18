"""Modified-hyperbolic base class for Darendeli, Menq, and Rollins models."""

from __future__ import annotations

from abc import abstractmethod

import numpy as np
import numpy.typing as npt

from ..contracts import NonlinearSoilCurves


class ModifiedHyperbolicBase:
    """Compute mod-reduc and damping using the modified-hyperbolic + Masing approach.

    Subclasses must supply ``strain_ref``, ``curvature``, and ``masing_scaling``
    as properties.  Call ``super().__init__(name, unit_wt, damping_min, strains)``
    after setting the parameters that those properties depend on.
    """

    def __init__(
        self,
        name: str,
        unit_wt: float,
        damping_min: float,
        strains: npt.ArrayLike | None = None,
    ) -> None:
        self.name = name
        self._unit_wt = unit_wt

        if strains is None:
            strains = np.logspace(-6, -1.5, num=20)
        else:
            strains = np.asarray(strains, dtype=float)

        # Modified-hyperbolic G/Gmax
        mod_reduc = 1.0 / (1.0 + (strains / self.strain_ref) ** self.curvature)

        # Masing damping
        strains_pct = strains * 100
        strain_ref_pct = self.strain_ref * 100
        masing_a1 = (100.0 / np.pi) * (
            4.0
            * (
                strains_pct
                - strain_ref_pct
                * np.log((strains_pct + strain_ref_pct) / strain_ref_pct)
            )
            / (strains_pct**2 / (strains_pct + strain_ref_pct))
            - 2.0
        )
        a = self.curvature
        c1 = -1.1143 * a**2 + 1.8618 * a + 0.2523
        c2 = 0.0805 * a**2 - 0.0710 * a - 0.0095
        c3 = -0.0005 * a**2 + 0.0002 * a + 0.0003
        damping_masing = c1 * masing_a1 + c2 * masing_a1**2 + c3 * masing_a1**3

        d_correction = self.masing_scaling * damping_masing * mod_reduc**0.1
        damping = np.maximum.accumulate(d_correction / 100.0)

        if isinstance(damping_min, np.ndarray):
            damping = damping_min + damping[:, np.newaxis]
        else:
            damping += damping_min

        self._strains = strains
        self._mod_reduc = mod_reduc
        self._damping = damping
        self._damping_min = (
            damping_min if np.ndim(damping_min) == 0 else float(damping_min.flat[0])
        )

    @property
    @abstractmethod
    def strain_ref(self) -> float:
        """Reference strain [decimal]."""
        ...

    @property
    @abstractmethod
    def curvature(self) -> float:
        """Curvature exponent."""
        ...

    @property
    @abstractmethod
    def masing_scaling(self) -> float:
        """Scaling factor for the Masing damping component."""
        ...

    def curves(self) -> NonlinearSoilCurves:
        return NonlinearSoilCurves(
            strains=self._strains,
            mod_reduc=self._mod_reduc,
            damping=self._damping,
            damping_min=self._damping_min,
            unit_wt=self._unit_wt,
            name=self.name,
        )
