from pygmm import Scenario
from pygmm.afshari_stewart_2016 import AfshariStewart2016


def test_run():
    """Simple test to make sure we can get outputs."""

    s = Scenario(mag=6, dist_rup=50, v_s30=300, mechanism="SS")

    m = AfshariStewart2016(s)

    m.duration

    m.std_err
