"""Summary of the step 2 noise-tolerance diagnostic (NEXT_STEPS.md), one row per case and setting.

    python noise_summary.py LIST [--archive DIR]

LIST is the rck_diag.sbatch list (NAME SUBGRAPH TRUTH OUTDIR [RCK OPTIONS]); every run is a
--keep-pool run at the λ of DIAG_POOL_LAMBDA (0.2), so its rank 1 is that λ's answer. Runs without
pool/summary.json are reported as pending. Columns: rank-1 significant species, weights, steps,
residual and mean composite; the best mean composite in the pool, its pool rank and significant
species, and whether rank 1 is within 0.005 of it (``picked``); pool size and how many pool
solutions do not fit the integer flow exactly (``inexact``, residual > 0); paths loaded and the
read share rank 1 explains; runtime. Prints a TSV, and with --archive copies each finished run's
small outputs to DIR/<name>/ (the pool's full decompositions.tsv, several MB, stays behind).
"""
import argparse
import gzip
import json
import os
import shutil
import sys

import pandas as pd

KEEP = ["summary.json", "evaluation.tsv", "block_flips.tsv", "flow_segments.tsv", "flow_adjacencies.tsv",
        "pool/pool.tsv", "pool/pool.json", "pool/evaluation.tsv"]
GZIP = ["decompositions.tsv"]


def mean_composites(evaluation):
    return evaluation.groupby("decomposition_rank")["composite"].mean()


def significant(decompositions, rank):
    d = decompositions[decompositions["decomposition_rank"] == rank].drop_duplicates("species")
    return d[d["status"] == "species"]


def row(name, options, run):
    pool = os.path.join(run, "pool")
    out = {"name": name, "options": " ".join(options) or "(defaults)"}
    if not os.path.exists(os.path.join(pool, "summary.json")):
        return dict(out, status="pending")
    summary = json.load(open(os.path.join(pool, "summary.json")))
    dec = pd.read_csv(os.path.join(pool, "decompositions.tsv"), sep="\t")
    ev = pd.read_csv(os.path.join(pool, "evaluation.tsv"), sep="\t")
    pev = pd.read_csv(os.path.join(pool, "pool", "evaluation.tsv"), sep="\t")
    ptab = pd.read_csv(os.path.join(pool, "pool", "pool.tsv"), sep="\t")
    sig = significant(dec, 1)
    rank1 = ev[ev["decomposition_rank"] == 1]
    means = mean_composites(pev)
    best_rank = int(means.idxmax())
    best_sig = ptab.loc[ptab["pool_index"] == best_rank - 1, "n_significant"]
    mean1 = float(rank1["composite"].mean())
    return dict(out, status="done",
                n_significant=len(sig), n_truth=int(rank1["n_species_truth"].iloc[0]),
                weights=",".join("{w:g}".format(w=w) for w in sig["weight"]), steps=int(sig["n_steps"].sum()),
                residual=round(float(dec.loc[dec["decomposition_rank"] == 1, "residual"].iloc[0]), 3),
                composites=",".join("{t}={c:.3f}".format(t=t, c=c) for t, c in zip(rank1["truth"], rank1["composite"])),
                mean_composite=round(mean1, 4),
                best_pool_composite=round(float(means.max()), 4), best_pool_rank=best_rank,
                best_pool_significant=int(best_sig.iloc[0]) if len(best_sig) else None,
                picked=mean1 >= float(means.max()) - 0.005,
                pool_size=len(ptab), inexact=int((ptab["residual"] > 1e-6).sum()),
                decided_by=(summary.get("ambiguity") or {}).get("decided_by", ""),
                paths_loaded=summary.get("paths_loaded"),
                reads_explained=dec.loc[dec["decomposition_rank"] == 1, "reads_explained"].iloc[0],
                runtime_s=summary.get("runtime_s"))


def archive(run, target):
    pool = os.path.join(run, "pool")
    for rel in KEEP:
        src = os.path.join(pool, rel)
        if os.path.exists(src):
            dst = os.path.join(target, rel)
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copyfile(src, dst)
    for rel in GZIP:
        src = os.path.join(pool, rel)
        if os.path.exists(src):
            with open(src, "rb") as fin, gzip.open(os.path.join(target, rel + ".gz"), "wb") as fout:
                shutil.copyfileobj(fin, fout, length=4 << 20)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("list")
    parser.add_argument("--archive", default=None, help="Copy each finished run's small outputs to ARCHIVE/<name>/")
    args = parser.parse_args()
    rows = []
    for line in open(args.list):
        if not line.strip():
            continue
        name, _sub, _truth, run, *options = line.split()
        rows.append(row(name, options, run))
        if args.archive and rows[-1]["status"] == "done":
            archive(run, os.path.join(args.archive, name))
    pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
