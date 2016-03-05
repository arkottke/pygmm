"""Create plots for visual comparison."""

__author__ = 'akottke'

import os

import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
import numpy as np

from ..atkinson_boore_2006 import AtkinsonBoore2006
from ..boore_stewart_seyhan_atkinson_2014 import  BooreStewartSeyhanAtkinson2014
from ..campbell_2003 import Campbell2003
from ..campbell_bozorgnia_2014 import CampbellBozorgnia2014
from ..chiou_youngs_2014 import ChiouYoungs2014
from ..idriss_2014 import Idriss2014
from ..pezeshk_zandieh_tavakoli_2011 import PezeshkZandiehTavakoli2011
from ..tavakoli_pezeshk_2005 import TavakoliPezeshk05

DEFAULT_PROPS = dict(
    depth_2_5=5,
    depth_bor=15,
    depth_hyp=9,
    depth_tor=5,
    dip=90.,
    dist=20,
    dist_jb=30.,
    dist_rup=30.,
    dist_x=30.,
    dpp_centered=0,
    flag_hw=0,
    flag_meas=0,
    mag=6,
    mechanism='SS',
    region='california',
    v_s30=500.,
    width=10,
)

# Make the figure directory if needed
if not os.path.exists('figures'):
    os.makedirs('figures')


def plot_model_with_param(model, key, values, label):
    props = dict(DEFAULT_PROPS)

    fig, ax = plt.subplots()
    for v in values:
        if isinstance(key, str):
            props[key] = v
        else:
            for k in key:
                props[k] = v
        m = model(**props)
        ax.plot(m.periods, m.spec_accels, label='%g' % v)

    ax.set_xlabel('Period (s)')
    try:
        ax.set_xscale('log')
    except ValueError:
        pass

    ax.set_ylabel('5% Damped, Spectral. Accel. (g)')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1e1)

    ax.legend(loc='upper right', title=label, fontsize='x-small')

    ax.grid()

    fig.tight_layout()

    if isinstance(key, str):
        prefix = key
    else:
        prefix = key[0]

    fig.savefig(os.path.join('figures', prefix + '-' + model.ABBREV))
    plt.close(fig)


def test_models():
    for m in [AtkinsonBoore2006, BooreStewartSeyhanAtkinson2014,
              Campbell2003, CampbellBozorgnia2014, ChiouYoungs2014,
              Idriss2014, PezeshkZandiehTavakoli2011, TavakoliPezeshk05]:
        for key, values, label in [
                ('mag', [5, 6, 7], 'Magnitude'),
                (['dist', 'dist_rup', 'dist_jb', 'dist_x'], [10, 50, 100],
                 'Distance (km)'),
                ('v_s30', [300, 650, 1000], '$V_{s30}$ (m/s)'),
                ]:
            yield plot_model_with_param, m, key, values, label

