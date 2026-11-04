"""Re-score the bench through the method-agnostic AA-format path and compare with today's scores.

    python diag/aa_regression.py [--out DIR] [--lambdas 0 0.2] [--jobs 16] [--cases demo3 ...]

For every case of run_case.py and each λ (0: the lambda_0 run's rank 1; λ > 0: the pool's rank 1
at λ, as collect.py re-ranks it):

1. the truth mixture graph is read from truth/aa/mixture_* when it is there (the untracked
   simulations on scratch or /n/fs/yubi-iw/runs, and the tracked demo5-demo8 runs, whose
   mixtures are in amplicon-allele-analysis);
   otherwise it is built with amplicon-allele-analysis simulations/truth_mixture.py, in place when
   the truth folder is on scratch or /n/fs/yubi-iw/runs, else under OUT/mixtures/ (demo4);
   runs without truth/aa (demo3) get no graph-level metrics;
2. the answer and its first 2 same-fit alternatives are exported to AA format
   (rck.ecdna.aa_export) under OUT/<case>/lambda_<λ>/ and scored by rck.ecdna.aa_eval;
3. the strict scores are compared with rck-ecdna-eval's (lambda_0/recovery_matrix.tsv, or
   pool/pool/evaluation.tsv at the pool entry), field by field at full precision (read with
   float_precision="round_trip": pandas' default parser is off by one ulp on some values), and with
   bench_sweep.tsv (composites, mean_composite, n_species, n_truth, merged).

Writes OUT/aa_regression.tsv (one row per case and λ: every compared field and whether it
matches), OUT/aa_metrics.tsv (the new metrics, strict and best-of) and OUT/aa_species.tsv.
"""
import argparse
import os
import sys
import traceback
from multiprocessing import Pool

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_case  # noqa: E402  (puts the RCK root on sys.path)
from rck.ecdna import aa_eval, aa_export  # noqa: E402
from rck.ecdna import evaluate as ev  # noqa: E402
from rck.ecdna.rerank import rerank  # noqa: E402

sys.path.insert(0, run_case.SIM)
import truth_mixture  # noqa: E402  (amplicon-allele-analysis simulations/truth_mixture.py)

DEFAULT_OUT = "/n/fs/yubi-iw/runs/bench_nested_aa"
PARTS = ["inferred", "segment_jaccard", "haplotype_jaccard", "junction_f1", "phased_junction_f1", "phasing_accuracy",
         "order_match", "cn_rel_err", "composite"]
OLD_NAMES = {"segment_jaccard": "jaccard_total", "haplotype_jaccard": "jaccard_hap"}
JUNCTION_MATCH = 0.5  # collect.py's merged rule


def truth_sample(case):
    opts = case["options"]
    return opts[opts.index("--truth-sample") + 1] if "--truth-sample" in opts else None


def mixture_dir(case, out):
    """Build (once) and return the case's truth mixture folder, or None without truth/aa."""
    truth = case["truth"]
    if not os.path.isdir(os.path.join(truth, "truth", "aa")):
        return None
    in_place = os.path.join(truth, "truth", "aa")
    if os.path.exists(os.path.join(in_place, "mixture_amplicons.tsv")):
        return in_place  # written there already: the untracked simulations, and the tracked demo5-demo8 runs
    untracked = os.path.realpath(truth).startswith(("/n/fs/scratch/", "/n/fs/yubi-iw/runs/"))
    target = in_place if untracked else os.path.join(out, "mixtures", case["name"])
    if not os.path.exists(os.path.join(target, "mixture_amplicons.tsv")):
        try:
            truth_mixture.write_all(truth, target, truth_sample(case))
        except (ValueError, OSError) as error:
            print("no truth mixture for {c}: {e}".format(c=case["name"], e=error), flush=True)
            return None
    return target


def merged(truth, strict, rows):
    """collect.py's merged flag on the AA-parsed cycles."""
    carries = {c.name: [t.name for t in truth if ev.junction_scores(t.steps, c.steps, 10_000, False, t.circular, c.circular)[0] >= JUNCTION_MATCH]
               for c in strict}
    notes = ["{c} carries {t}".format(c=c, t="+".join(ts)) for c, ts in carries.items() if len(ts) >= 2]
    for t in rows.loc[rows["inferred"].isna(), "truth"]:
        holders = [c for c, ts in carries.items() if t in ts]
        if holders:
            notes.append("{t} missed, carried by {c}".format(t=t, c="+".join(holders)))
    return bool(notes)


def same(a, b):
    if (a is None or (isinstance(a, float) and np.isnan(a))) and (b is None or (isinstance(b, float) and np.isnan(b))):
        return True
    if isinstance(a, (bool, np.bool_)) or isinstance(b, (bool, np.bool_)):
        return bool(a) == bool(b)
    try:
        return float(a) == float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)


def one(task):
    case, lam, out = task
    base = os.path.join(run_case.OUT, case["name"])
    try:
        if lam == 0:
            run_dir, index = os.path.join(base, "lambda_0"), None
            if not run_case.done(run_dir):
                return None
            old = pd.read_csv(os.path.join(run_dir, "recovery_matrix.tsv"), sep="\t", float_precision="round_trip")
            old = old[old["in_scope"].astype(bool)]
        else:
            run_dir = os.path.join(base, "pool")
            if not run_case.done(run_dir, pool=True):
                return None
            index = int(rerank(run_dir, [lam])["pool_index"].iloc[0])
            old = pd.read_csv(os.path.join(run_dir, "pool", "evaluation.tsv"), sep="\t", float_precision="round_trip")
            old = old[old["decomposition_rank"] == index + 1].rename(columns={v: k for k, v in OLD_NAMES.items()})
        folder = os.path.join(out, case["name"], "lambda_%g" % lam)
        paths = aa_export.export(run_dir, folder, index, alternatives=2)
        truth = ev.load_truth(case["truth"], None, truth_sample(case))
        mdir = mixture_dir(case, out)
        mixture = aa_eval.load_mixture(mdir) if mdir else None
        subgraph_paths = os.path.join(case["subgraph"], "paths.tsv")
        result = aa_eval.score(paths["cycles"], truth, paths.get("graph"), paths["alt_cycles"], mixture, paths.get("segment_cn"),
                               subgraph_paths if os.path.exists(subgraph_paths) else None)
        aa_eval.write(folder, result)
        strict = result.species[result.species["mode"] == "strict"].reset_index(drop=True)
        row = {"case": case["name"], "kind": case["kind"], "scenario": case["scenario"], "step_penalty": lam,
               "pool_index": index, "n_alternatives_exported": len(paths["alt_cycles"])}
        mismatches = []
        if list(old["truth"]) != list(strict["truth"]):
            mismatches.append("truth species {a} vs {b}".format(a=list(old["truth"]), b=list(strict["truth"])))
        else:
            for col in PARTS:
                for a, b, t in zip(old[col], strict[col], strict["truth"]):
                    if not same(a, b):
                        mismatches.append("{t}.{c}: {a} vs {b}".format(t=t, c=col, a=a, b=b))
        row["composites"] = ",".join("%s=%.3f" % (r.truth, r.composite) for r in strict.itertuples())
        row["mean_composite"] = round(float(strict["composite"].mean()), 4) if len(strict) else None
        row["n_species"] = len(result.strict)
        row["n_truth"] = len(strict)
        row["merged"] = merged(result.truth, result.strict, strict)
        row["parts_mismatches"] = "; ".join(mismatches)
        new = result.sample.set_index("mode")
        for mode in ("strict", "best_of"):
            for col, v in new.loc[mode].items():
                if col != "mode":
                    row["%s.%s" % (mode, col)] = v
        species = result.species.assign(case=case["name"], step_penalty=lam)
        return row, species, result.graph.assign(case=case["name"], step_penalty=lam)
    except Exception:  # noqa: BLE001 - report every failing case, keep going
        return {"case": case["name"], "step_penalty": lam, "error": traceback.format_exc()}, None, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=DEFAULT_OUT)
    ap.add_argument("--lambdas", type=float, nargs="+", default=[0.0, 0.2])
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--cases", nargs="*", default=None)
    a = ap.parse_args()
    cases = [c for c in run_case.cases() if a.cases is None or c["name"] in a.cases]
    tasks = [(c, lam, a.out) for c in cases for lam in a.lambdas]
    os.makedirs(a.out, exist_ok=True)
    # build mixtures first, one process per case, so parallel tasks never write the same files
    with Pool(min(a.jobs, len(cases))) as pool:
        pool.starmap(mixture_dir, [(c, a.out) for c in cases])
        results = [r for r in pool.map(one, tasks) if r is not None]
    rows = [r[0] for r in results]
    table = pd.DataFrame(rows)
    sweep = pd.read_csv(os.path.join(os.path.dirname(HERE), "bench_sweep.tsv"), sep="\t")
    checks = []
    for r in table.to_dict("records"):
        if r.get("error"):
            checks.append({"sweep_match": False, "sweep_mismatches": "error"})
            continue
        s = sweep[(sweep["case"] == r["case"]) & np.isclose(sweep["step_penalty"], r["step_penalty"])]
        if not len(s):
            checks.append({"sweep_match": None, "sweep_mismatches": "not in bench_sweep.tsv"})
            continue
        s = s.iloc[0]
        bad = [c for c in ("composites", "mean_composite", "n_species", "n_truth", "merged") if not same(r[c], s[c])]
        checks.append({"sweep_match": not bad, "sweep_mismatches": ",".join("%s: %s vs %s" % (c, r[c], s[c]) for c in bad)})
    table = pd.concat([table, pd.DataFrame(checks)], axis=1)
    table.to_csv(os.path.join(a.out, "aa_regression.tsv"), sep="\t", index=False)
    species = pd.concat([r[1] for r in results if r[1] is not None], ignore_index=True)
    species.to_csv(os.path.join(a.out, "aa_species.tsv"), sep="\t", index=False)
    graphs = [r[2] for r in results if r[2] is not None and len(r[2])]
    if graphs:
        pd.concat(graphs, ignore_index=True).to_csv(os.path.join(a.out, "aa_graph.tsv"), sep="\t", index=False)
    ok = table["sweep_match"].fillna(False).astype(bool) & (table.get("parts_mismatches", "").fillna("") == "")
    print("{n} (case, λ) scored; {k} identical to bench_sweep.tsv and to rck-ecdna-eval's parts; {e} errors".format(
        n=len(table), k=int(ok.sum()), e=int(table.get("error", pd.Series(dtype=str)).notna().sum())))
    bad = table[~ok]
    if len(bad):
        cols = [c for c in ("case", "step_penalty", "parts_mismatches", "sweep_mismatches", "error") if c in bad]
        with pd.option_context("display.max_colwidth", 300, "display.width", 250):
            print(bad[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
