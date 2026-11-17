#!/usr/bin/env python3
"""bench_IW_allele Part B, step 4: score the multi-clone route on RCK's degraded inputs (PLAN.md "Part B").

One row per (case, variant), variants ``labelled`` (the existing bench_IW_intpsi rck_mc run on the labelled
input: the reference) and ``chil`` / ``chil_avg`` / ``chil_err`` (rck_partB.py runs, exported with
``rck_partB.py export``). Every variant is scored by the same code:

  rck_eval.score      novel-adjacency FNR / FDR (unlabelled; labelled within phase blocks), bulk Eq. 11,
                      generalized-IS violations, linear walks
  input use           input adjacencies used, and RCK's FPR over the spurious ones (chil_err)
  clone_cn            Eq. 11 summed over clones and copy-number state agreement, in two frames:
                      ``pb`` (phase blocks, the labelled input's fragments: the common frame) and
                      ``frag`` (the variant's own input fragments: RCK's Eq. 11 as written)
  allele_metrics      allele-specific ecDNA copy-number recall / BAF error within phase blocks (items 2-3)

Why phase blocks as the common frame: CHIL flips A/B per fragment (per segment in ``chil``); swapping the
answer per fragment to match the truth would forgive exactly the phasing errors CHIL tests. Phase across
phase blocks is not identifiable from any data; within a block the model has to reconstruct it from
copy-number balance and junctions, so that is what the main scores measure.

    python score_partB.py run [--jobs 16]     # finished tasks are skipped
    python score_partB.py collect             # partB_table.tsv here

Outputs: <run>/score/partB.json per run (the labelled reference's under
/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/<case>/labelled/score/partB.json). Local, no Gurobi.
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
import make_rck_inputs as mi  # noqa: E402
import rck_partB as rb  # noqa: E402

bi = ab.bi
VARIANTS = ["labelled", "chil", "chil_avg", "chil_err"]
TABLE = os.path.join(HERE, "partB_table.tsv")


def paths(case, variant):
    """(run dir, rck_input dir, answer cycles, output json) of one (case, variant)."""
    if variant == "labelled":
        run = bi.run_dir(case, "rck_mc")
        cycles = os.path.join(bi.score_folder(case, "rck_mc"), "answer", "{c}_rck_mc_cycles.txt".format(c=case["name"]))
        out = os.path.join(mi.OUT, case["name"], "labelled", "score", "partB.json")
        return run, case["rck_input"], cycles, out
    run = rb.run_dir(case, variant)
    cycles = os.path.join(run, "score", "answer", "{c}_rck_mc_cycles.txt".format(c=case["name"]))
    return run, mi.input_dir(case, variant), cycles, os.path.join(run, "score", "partB.json")


def _count(value):
    return len(value) if isinstance(value, (list, dict)) else value


def tasks():
    return [(case, variant) for case in bi.cases() for variant in VARIANTS]


def score_one(task):
    case, variant = task
    run, inp, cycles, out = paths(case, variant)
    if os.path.exists(out):
        return case["name"], variant, "done already"
    if not os.path.exists(cycles):
        return case["name"], variant, "no answer"
    try:
        from rck.ecdna import allele_metrics as am, bench_metrics as bm, evaluate as ev, rck_eval
        truth = ev.load_truth(case["truth"])
        props = bm.read_proportions(os.path.join(inp, "clone_proportions.tsv"))
        weights = bm.multiclone_bulk_weights(run, props)
        strict, _ = am.strict_cycles(cycles, weights)
        pb_path = os.path.join(case["rck_input"], "rck.fragments.tsv")  # phase blocks: the common frame
        blocks = ev.read_phase_blocks(pb_path)
        result = rck_eval.score(strict, truth, blocks)
        result["input"] = rck_eval.multiclone_input_use(run, inp, truth)
        truth_scnt = os.path.join(case["truth"], "truth", "rck", "rck.scnt.tsv")
        regions = rck_eval.ecdna_regions(truth)
        for frame, fragments in (("pb", pb_path), ("frag", os.path.join(inp, "rck.fragments.tsv"))):
            scores = rck_eval.clone_cn_scores(os.path.join(run, "rck.scnt.tsv"), truth_scnt, fragments, regions)
            scores.pop("clones", None)
            result["clone_cn_" + frame] = scores
        allele, _ = am.score(cycles, truth, weights, blocks)
        result["allele"] = {k: allele.get(k) for k in ("as_cn_recall_block", "as_cn_precision_block", "baf_mae_block",
                                                        "as_cn_mae_block", "blind_cn_recall", "mixed_bp_frac")}
        with open(os.path.join(run, "summary.json")) as handle:
            s = json.load(handle)
        result["run"] = {"status": s.get("status"), "mip_gap": s.get("mip_gap"), "runtime_rck_s": s.get("runtime_rck_s"),
                         "runtime_s": s.get("runtime_s"), "fragments": s.get("fragments"),
                         "fragments_flipped": len(s.get("fragments_flipped", [])),
                         "n_cycles": _count(s["decomposition"].get("cycles")), "n_paths": _count(s["decomposition"].get("paths"))}
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w") as handle:
            json.dump(result, handle, indent=2)
        return case["name"], variant, "ok"
    except Exception:  # noqa: BLE001 - recorded, the other tasks go on
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(os.path.splitext(out)[0] + "_error.txt", "w") as handle:
            handle.write(traceback.format_exc())
        return case["name"], variant, "error"


def run(jobs):
    import logging
    logging.disable(logging.INFO)
    states = {}
    with Pool(jobs) as pool:
        for name, variant, state in pool.imap_unordered(score_one, tasks()):
            states[state] = states.get(state, 0) + 1
            if state == "error":
                print(name, variant, state, flush=True)
    print(states)


def collect():
    from rck.ecdna import rck_eval
    rows = []
    for case, variant in tasks():
        out = paths(case, variant)[3]
        row = {k: case[k] for k in ("name", "class", "S", "replicate")}
        row["variant"] = variant
        if os.path.exists(out):
            with open(out) as handle:
                row.update(rck_eval.flatten(json.load(handle)))
            row["status"] = "ok"
        else:
            row["status"] = "missing"
        rows.append(row)
    frame = pd.DataFrame(rows)
    allele = pd.read_csv(ab.TABLE, sep="\t")
    mixed = allele[allele["phasing"] == "native"].groupby("name")["mixed_bp_frac"].first()
    frame["group"] = frame["name"].map(lambda n: "mixed" if mixed.get(n, 0) > 0 else "single")
    frame.to_csv(TABLE, sep="\t", index=False)
    print("{n} rows ({ok} ok) -> {t}".format(n=len(frame), ok=int((frame["status"] == "ok").sum()), t=TABLE))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=16)
    sub.add_parser("collect")
    a = p.parse_args()
    if a.cmd == "run":
        run(a.jobs)
    else:
        collect()


if __name__ == "__main__":
    main()
