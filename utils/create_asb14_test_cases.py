#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create test cases for the DBC13 model."""
import gzip
import json
import os

import xlwings as xw

params = {
    "dist": [5, 40, 200],
    "mag": [4, 6, 7],
    "mechanism": ["NS", "RS", "SS"],
    "v_s30": [200, 400, 750, 1200],
}


def iter_parameters(parameters):
    max_len = max(len(p) for p in parameters.values())

    for i in range(max_len):
        yield {k: v[i % len(v)] for k, v in parameters.items()}


def reformat_params(params):
    """Translate to spreadsheet values."""
    params["mechanism"] = ["NS", "RS", "SS"].index(params["mechanism"])
    return params


def load_params(wb, **params):
    cell_map = {
        "mag": "B1",
        "v_s30": "B6",
        "mechanism": "B5",
    }
    params = reformat_params(dict(params))
    for k, v in params.items():
        if k == "dist":
            # Update all distances
            for c in ["B2", "B3", "B4"]:
                xw.Range("Input-Output", c, wkb=wb).value = v
        else:
            xw.Range("Input-Output", cell_map[k], wkb=wb).value = v

    # Force a re-calculation. This is the same as Ctrl + Alt + F9
    wb.xl_app.CalculateFull()


def get_results(wb):
    results = []
    for col, dist in zip("GHI", ["jb", "epi", "hyp"]):
        d = {}
        # Collect PGA and PGV
        for key, row in zip(["pga", "pgv"], [3, 4]):
            d[key] = xw.Range("Input-Output", "%s%d" % (col, row), wkb=wb).value

        # Collect response spectrum
        for key, c in zip(["periods", "spec_accels"], ["F", col]):
            d[key] = xw.Range("Input-Output", "%s5:%s66" % (c, c), wkb=wb).value

        results.append(("dist_" + dist, d))

    return results


wb_fname = os.path.abspath("ASB14.xlsx")
tests = []
for p in iter_parameters(params):
    wb = xw.Workbook(wb_fname)
    load_params(wb, **p)
    tests.append({"params": p, "results": get_results(wb)})

fname = "../tests/data/asb14_tests.json.gz"
with gzip.open(fname, "wt") as fp:
    json.dump(tests, fp, indent=4)
