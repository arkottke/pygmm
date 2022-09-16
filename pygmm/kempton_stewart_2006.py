"""Kempton and Stewart (2006, :cite:`kempton2006`) duration model."""

from __future__ import division

import numpy as np

from . import model

__author__ = ''


class KemptonStewart2006(model.Model):
    """Kempton and Stewart (2006, :cite:`kempton2006`) duration model.

    Parameters
    ----------
    scenario : :class:`pygmm.model.Scenario`
        earthquake scenario

    """

    NAME = 'Kempton Stewart (2006)'
    ABBREV = 'KS06'

    PARAMS = [
        model.NumericParameter('mag',       True,   5,   7.6),
        model.NumericParameter('dist_rup',  True,   0,   200),
        model.NumericParameter('v_s30',     True,   200, 1000),
    ]

    def __init__(self, scenario):
        super(KemptonStewart2006, self).__init__(scenario)

        #scenario
        s = self._scenario

       #stress drop indices for duration measures
        stress_drop = np.exp([6.02,                         #D5-75a
                              2.79 + 0.82 * (s.mag - 6),    #D5-95a
                              5.46,                         #D5-75v
                              1.53 + 1.34 * (s.mag - 6)])   #D5-95v
        

        #source duration
        f_0 = np.array([self.source_dur(s.mag, s_d) for s_d in stress_drop])
        
        #path duration
        f_1 = np.array([0.07, 0.15, 0.10, 0.15]) * s.dist_rup
        
        #site duration
        f_2 = (np.array([0.82, 3.00, 1.40, 3.99]) + 
               np.array([-0.0013,-0.0041,-0.0022,-0.0062]) * s.v_s30 )
        
        #total druation
        self._ln_dur = np.log( f_0 + f_1 + f_2 )
        
        #aleatory standard deviation
        self._std_err = np.array([0.57, 0.47, 0.66, 0.52])

    @property
    def acc_duration(self):
        return np.exp(self._ln_dur[:2])

    @property
    def vel_duration(self):
        return np.exp(self._ln_dur[2:])

    @property
    def acc_std_err(self):
        return self._std_err[:2]

    @property
    def vel_std_err(self):
        return self._std_err[2:]

    def source_dur(self, mag, stress_drop):
        
        #seismic moment
        moment = 10 ** (1.5 * mag + 16.05)
        
        return (stress_drop / moment) ** (-1/3) / (4.9e6 * 3.2)
