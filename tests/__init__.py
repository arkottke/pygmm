import gzip
import json
import pathlib


FPATH_DATA = pathlib.Path(__file__).parent / "data"


def load_tests(fname):
    fpath = FPATH_DATA / fname
    if fpath.suffix == ".gz":
        fp = gzip.open(fpath, "rt")
    else:
        fp = fpath.open()

    return json.load(fp)
