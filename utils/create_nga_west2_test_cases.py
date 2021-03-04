#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Create test cases for the NGA West2 models."""
import gzip
import json
import os

import xlwings as xw


def iter_parameters(parameters):
    max_len = max(len(p) for p in parameters.values())

    for i in range(max_len):
        yield {k: v[i % len(v)] for k, v in parameters.items()}


def reformat_params(params):
    """Transform from pygmm input to spreadsheet input."""
    if "mechanism" in params:
        # Replace mechanism with flags
        mech = params.pop("mechanism")
        if mech == "SS":
            params["flag_u"] = 0
            params["flag_rv"] = 0
            params["flag_nm"] = 0
        elif mech == "NS":
            params["flag_u"] = 0
            params["flag_rv"] = 0
            params["flag_nm"] = 1
        elif mech == "RS":
            params["flag_u"] = 0
            params["flag_rv"] = 1
            params["flag_nm"] = 0
        elif mech == "U":
            params["flag_u"] = 1
            params["flag_rv"] = 0
            params["flag_nm"] = 0

    try:
        on_hanging_wall = params.pop("on_hanging_wall")
        params["flag_hw"] = 1 if on_hanging_wall else 0
    except KeyError:
        pass

    try:
        params["flag_meas"] = params.pop("vs_source")
    except KeyError:
        pass

    try:
        params["region"] = params["region"].replace("_", " ").title()
    except KeyError:
        pass

    # Replace None with 999
    for key in ["depth_tor", "depth_hyp", "depth_1_0", "depth_2_5", "width", "dist_y0"]:
        if key in params and params[key] is None:
            params[key] = 999

    if "dist_rup" not in params and "dist_jb" in params:
        params["dist_rup"] = params["dist_jb"]

    try:
        # No supported to changing the DPP
        params.pop("dpp_centered")
    except KeyError:
        pass

    return params


def load_params(wb, **params):
    params = reformat_params(dict(params))
    for k, v in params.items():
        print("Setting:", k, CELLS[k], v)
        xw.Range("Main", CELLS[k], wkb=wb).value = v

    # Force a re-calculation. This is the same as Ctrl + Alt + F9
    wb.xl_app.CalculateFull()


def get_results(wb, abbrev):
    sheetname = abbrev
    result_columns = RESULT_COLUMNS[abbrev]
    # Collect response spectrum
    d = {}
    for key, rc in zip(["periods", "spec_accels", "ln_stds"], result_columns):
        d[key] = xw.Range(sheetname, "{col}6:{col}26".format(col=rc), wkb=wb).value

    # Collect the PGA and PGV
    for key, row in zip(["pga", "pgv"], [28, 29]):
        values = xw.Range(
            sheetname,
            "{col0:}{row:d}:{col1:}{row:d}".format(
                col0=result_columns[1], col1=result_columns[2], row=row
            ),
            wkb=wb,
        ).value
        for subkey, v in zip(["", "_ln_std"], values):
            try:
                float(v)
            except ValueError:
                v = None
            d[key + subkey] = v

    return d


wb_fname = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__), "NGAW2_GMPE_Spreadsheets_v5.7_041415_Protected.xlsm"
    )
)

# Parameters
params_site = {
    "v_s30": [300, 450, 760, 1000],
    "depth_1_0": [0.4, 0.3, 0.04, 0.005],
    "depth_2_5": [1.8, 1.1, 0.4, 0.3],
    "vs_source": ["measured", "inferred"],
    "region": [
        "global",
        "california",
        "china",
        "japan",
        "italy",
        "new_zealand",
        "turkey",
        "taiwan",
    ],
}

params_fault = [
    {
        "mag": 6.0,
        "dist_rup": 3.16,
        "dist_jb": 1.0,
        "dist_x": 1.0,
        "dist_y0": None,
        "mechanism": "SS",
        "on_hanging_wall": False,
        "dip": 90,
        "depth_tor": 0,
        "depth_hyp": 8,
        "width": 10,
    },
    {
        "mag": 7.1,
        "dist_rup": 10,
        "dist_jb": 2,
        "dist_x": 16,
        "dist_y0": None,
        "mechanism": "RS",
        "on_hanging_wall": True,
        "dip": 50,
        "depth_tor": 1,
        "depth_hyp": 8,
        "width": 15,
    },
    {
        "mag": 4.5,
        "dist_rup": 10,
        "dist_jb": 8.0,
        "dist_x": -8,
        "dist_y0": None,
        "mechanism": "RS",
        "on_hanging_wall": False,
        "dip": 60,
        "depth_tor": 4,
        "depth_hyp": 5,
        "width": 2,
    },
    {
        "mag": 6.0,
        "dist_rup": 5,
        "dist_jb": 0.0,
        "dist_x": 0.1,
        "dist_y0": None,
        "mechanism": "NS",
        "on_hanging_wall": True,
        "dip": 60,
        "depth_tor": 4,
        "depth_hyp": 5,
        "width": 8,
    },
    {
        "mag": 6.5,
        "dist_rup": 100,
        "dist_jb": 100,
        "dist_x": 60,
        "dist_y0": 70,
        "mechanism": "U",
        "on_hanging_wall": False,
        "dip": 90,
        "depth_tor": 2,
        "depth_hyp": 7,
        "width": 15,
    },
]

CELLS = dict(
    mag="B24",
    dist_rup="B27",
    dist_jb="B30",
    dist_x="B33",
    dist_y0="B36",
    v_s30="B39",
    flag_u="B42",
    flag_rv="B45",
    flag_nm="B48",
    flag_hw="B51",
    dip="B54",
    depth_tor="B57",
    depth_hyp="B60",
    depth_1_0="B63",
    depth_2_5="B66",
    width="B69",
    flag_meas="B72",
    flag_as="B75",
    region="B78",
    dpp_centered="B83",
    depth_bor="B89",
)

RESULT_COLUMNS = dict(
    ASK14="DEF",
    BSSA14="DST",
    CB14="DYZ",
    CY14="DUV",
    I14="DLM",
)

abbreviations = ["ASK14", "BSSA14", "CB14", "CY14", "I14"]

# Join the parameters of the site and fault
parameters = dict(params_site)
for pf in params_fault:
    for k, v in pf.items():
        if k not in parameters:
            parameters[k] = list()
        parameters[k].append(v)

tests = []
for p in iter_parameters(parameters):
    wb = xw.Workbook(wb_fname)
    load_params(wb, **p)

    results = {abb: get_results(wb, abb) for abb in abbreviations}

    if p["region"] == "taiwan":
        # Load ASK14 specific values.
        load_params(wb, region="TWA")
        results["ASK14"] = get_results(wb, "ASK14")

    tests.append({"params": p, "results": results})

fname = "../tests/data/ngaw2_tests.json.gz"
with gzip.open(fname, "wt") as fp:
    json.dump(tests, fp, indent=4)
