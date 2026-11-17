#!/usr/bin/env python3
"""bench_IW_allele Part B, step 1: RCK's simulated-input degradations of every case's multi-clone input.

For each bench_IW_intpsi case and variant (chil, chil_avg, chil_err; rck.ecdna.rck_inputs), writes a degraded
copy of the case's labelled DIAG folder (subgraphs/bench_IW_intpsi_cnp/<case>/union/rck_input/, only read) and
checks it loads with the multi-clone route's reader (rck.ecdna.multiclone.load_rck_input: RCK's own readers
and alignment checks; no Gurobi). Seeds: [case seed, variant index], so every input is reproducible.

    python make_rck_inputs.py [--jobs 8]      # existing outputs are rebuilt (cheap, deterministic)
    python make_rck_inputs.py --summary       # rck_inputs_summary.tsv here, from the transforms.json files

Outputs: /n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/<case>/<variant>/rck_input/ (+ transforms.json).
"""
import argparse
import json
import os
import sys
import traceback
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allele_bench as ab  # noqa: E402  (configures bench_iw for iw_intpsi)

bi = ab.bi
OUT = os.environ.get("BENCH_IW_RCKSTYLE_OUT", "/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle")
SUMMARY = os.path.join(HERE, "rck_inputs_summary.tsv")


def input_dir(case, variant):
    return os.path.join(OUT, case["name"], variant, "rck_input")


def tasks():
    from rck.ecdna.rck_inputs import VARIANTS
    return [(case, variant) for case in bi.cases() for variant in VARIANTS]


def run_one(task):
    case, variant = task
    from rck.ecdna import multiclone, rck_inputs
    dst = input_dir(case, variant)
    try:
        rec = rck_inputs.degrade(case["rck_input"], dst, variant, seed=[int(case["seed"]), rck_inputs.VARIANTS.index(variant)])
        inp = multiclone.load_rck_input(dst)
        rec.update(loaded_segments=len(inp.segments), loaded_novel=len(inp.input_adjacencies),
                   loaded_groups=len(inp.adjacency_groups), loaded_fragments=len(set(inp.segments_to_fragments.values())))
        with open(os.path.join(dst, "transforms.json"), "w") as handle:
            json.dump(rec, handle, indent=2)
        return case["name"], variant, "ok"
    except Exception:  # noqa: BLE001 - recorded, the other tasks go on
        os.makedirs(dst, exist_ok=True)
        with open(os.path.join(dst, "error.txt"), "w") as handle:
            handle.write(traceback.format_exc())
        return case["name"], variant, "error"


def summary():
    rows = []
    for case, variant in tasks():
        path = os.path.join(input_dir(case, variant), "transforms.json")
        row = {"name": case["name"], "variant": variant}
        if os.path.exists(path):
            with open(path) as handle:
                rec = json.load(handle)
            for key in ("shifted", "spurious"):
                if isinstance(rec.get(key), list):
                    rec["n_" + key] = len(rec.pop(key))
            rec.pop("source", None)
            rec.pop("seed", None)
            row.update(rec)
        else:
            row["status"] = "missing"
        rows.append(row)
    pd.DataFrame(rows).to_csv(SUMMARY, sep="\t", index=False)
    print("{n} rows -> {s}".format(n=len(rows), s=SUMMARY))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--summary", action="store_true")
    a = p.parse_args()
    if a.summary:
        return summary()
    import logging
    logging.disable(logging.INFO)
    states = {}
    with Pool(a.jobs) as pool:
        for name, variant, state in pool.imap_unordered(run_one, tasks()):
            states[state] = states.get(state, 0) + 1
            if state != "ok":
                print(name, variant, state, flush=True)
    print(states)
    summary()


if __name__ == "__main__":
    main()
