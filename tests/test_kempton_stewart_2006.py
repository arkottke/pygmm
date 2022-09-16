from pygmm import Scenario
from pygmm.kempton_stewart_2006 import KemptonStewart2006


def test_run():
    """Simple test to make sure we can get outputs."""

    s = Scenario(
        mag=6,
        dist_rup=50,
        v_s30=300,
    )

    m = KemptonStewart2006(s)

    m.duration

    m.std_err
