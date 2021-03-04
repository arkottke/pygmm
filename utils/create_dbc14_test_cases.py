#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create test cases for the DBC14 model."""
import gzip
import json
import os

import xlwings as xw

params = {
    "depth_hyp": [10, 0, 25],
    "dist_jb": [5, 40, 200],
    "mag": [4, 6, 7],
    "mechanism": ["NS", "RS", "SS"],
    "v_s30": [200, 400, 800],
}


def iter_parameters(parameters):
    max_len = max(len(p) for p in parameters.values())

    for i in range(max_len):
        yield {k: v[i % len(v)] for k, v in parameters.items()}


def reformat_params(params):
    """Translate to spreadsheet values."""
    params["mechanism"] = 1 + ["NS", "RS", "SS"].index(params["mechanism"])
    return params


def load_params(wb, **params):
    cell_map = {
        "dist_jb": "B08",
        "mag": "B10",
        "v_s30": "B12",
        "depth_hyp": "B14",
        "mechanism": "B16",
    }
    params = reformat_params(dict(params))
    for k, v in params.items():
        print("Setting:", k, cell_map[k], v)
        xw.Range("ANN(PGA)", cell_map[k], wkb=wb).value = v

    # Force a re-calculation. This is the same as Ctrl + Alt + F9
    wb.xl_app.CalculateFull()


def get_results(wb):
    d = {}
    # Collect PGA and PGV
    for key, row in zip(["pga", "pgv"], [13, 10]):
        d[key] = xw.Range("ANN(PGA)", "O%d" % row, wkb=wb).value

    # Collect response spectrum
    for key, rc in zip(["periods", "spec_accels"], "NO"):
        d[key] = xw.Range("ANN(PGA)", "%s17:%s78" % (rc, rc), wkb=wb).value

    return d


wb_fname = os.path.abspath("10518_2013_9481_MOESM1_ESM.xlsx")
tests = []
for p in iter_parameters(params):
    wb = xw.Workbook(wb_fname)
    load_params(wb, **p)
    tests.append({"params": p, "results": get_results(wb)})

fname = "../tests/data/dbc14_tests.json.gz"
with gzip.open(fname, "wt") as fp:
    json.dump(tests, fp, indent=4)
