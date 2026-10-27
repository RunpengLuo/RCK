"""bench_nested RCK runs: one array task per case, two rck-ecdna runs per case.

    python run_case.py list                 # index, case, subgraph
    python run_case.py run 7                # case 7 (SLURM_ARRAY_TASK_ID)
    python run_case.py run demo8            # or by name

Per case, under runs/<case>/:
  lambda_0/   --step-penalty 0, the current default (a λ = 0 run solves a narrower pool, so it is
              run, not re-ranked);
  pool/       --step-penalty 0.1 --keep-pool: the widened pool, which every λ > 0 shares; λ > 0
              is re-ranked from pool/pool/ offline (python -m rck.ecdna.rerank).
A run is skipped when its outputs exist, so the array can be resubmitted. The academic WLS
licence allows 2 sessions, and a finished process's session stays counted for a few minutes, so
a failed run is retried up to 4 times, 2 min apart, and each task releases its session at the end.

Cases: the regression cases of results/step_sweep/ (same subgraphs and options), then every
accepted bench_nested seed of amplicon-allele-analysis configs/bench_nested/seeds.tsv, in
scenario order, with --block-flip unoriented. Every run has --truth-dir, --no-plots,
--threads 4 and --time-limit 1800, and the subgraph folder as its only input.

BENCH_NESTED_ROOT, as for amplicon-allele-analysis bench_nested.py, takes the bench cases from a
benchmark built elsewhere ($BENCH_NESTED_ROOT/configs/seeds.tsv, subgraphs/<run>/union and
runs/<run>); BENCH_NESTED_OUT puts the runs elsewhere.
"""
import json
import os
import sys
import time
import logging

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RCK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, RCK)
AAA = os.environ.get("AAA", os.path.expanduser("~/amplicon-allele-analysis"))
SIM = os.path.join(AAA, "simulations")
DATA = os.path.join(RCK, "tests", "ecdna", "data")
OUT = os.environ.get("BENCH_NESTED_OUT", os.path.join(HERE, "runs"))
BENCH_ROOT = os.environ.get("BENCH_NESTED_ROOT")
SEEDS = os.path.join(BENCH_ROOT, "configs", "seeds.tsv") if BENCH_ROOT else os.path.join(SIM, "configs", "bench_nested", "seeds.tsv")
SUBGRAPHS = os.path.join(BENCH_ROOT, "subgraphs") if BENCH_ROOT else os.path.join(SIM, "subgraphs", "bench_nested")
SIM_RUNS = os.path.join(BENCH_ROOT, "runs") if BENCH_ROOT else os.path.join(SIM, "runs", "bench_nested")
POOL_LAMBDA = "0.1"
LAMBDAS = [0.0, 0.02, 0.05, 0.1, 0.2, 0.5]
RETRIES = 5  # licence sessions of finished tasks can take a few minutes to free up
COMMON = ["--no-plots", "--threads", os.environ.get("SLURM_CPUS_PER_TASK", "4"), "--time-limit", "1800"]

# the step sweep's cases (results/step_sweep/sweep.py), unchanged
REGRESSION = [
    ("demo3", f"{DATA}/amplicon_0", f"{DATA}/demo3_truth", ["--depth-per-copy", "15", "--junction-reads-per-copy", "6"]),
    ("demo4_amp0", f"{DATA}/demo4/amplicon_0", f"{DATA}/demo4_truth",
     ["--depth-per-copy", "13.5", "--truth-sample", "tumor1", "--min-species-weight", "2"]),
    ("demo4_amp1", f"{DATA}/demo4/amplicon_1", f"{DATA}/demo4_truth",
     ["--depth-per-copy", "13.5", "--truth-sample", "tumor1", "--min-species-weight", "2"]),
    ("demo5", f"{SIM}/runs/demo5_cycle1_ecdna/subgraph/amplicon_0", f"{SIM}/runs/demo5_cycle1_ecdna", []),
    ("demo6_union", f"{SIM}/runs/demo6_cycle2_ecdna/subgraph_union/union", f"{SIM}/runs/demo6_cycle2_ecdna", []),
    ("demo7_blockflip", f"{SIM}/runs/demo7_nestedcycles_ecdna/subgraph_union_blocks/union", f"{SIM}/runs/demo7_nestedcycles_ecdna",
     ["--block-flip", "unoriented"]),
    ("demo7_totalonly", f"{SIM}/runs/demo7_nestedcycles_ecdna/subgraph_union/union", f"{SIM}/runs/demo7_nestedcycles_ecdna",
     ["--min-snps-allele", "100000000"]),
    ("demo8", f"{SIM}/runs/demo8_nestedcycles_ecdna/subgraph_union/union", f"{SIM}/runs/demo8_nestedcycles_ecdna", []),
]
VALIDATION = ["demo3", "demo7_totalonly", "demo8"]


def bench_cases():
    """(case, subgraph, truth dir, options, scenario) of every accepted bench_nested seed."""
    path = SEEDS
    if not os.path.exists(path):
        return []
    seeds = pd.read_csv(path, sep="\t")
    seeds = seeds[seeds["accepted"].astype(str).str.lower() == "true"]
    sys.path.insert(0, SIM)
    from make_scenarios import NESTED_SCENARIOS

    order = {s: i for i, s in enumerate(NESTED_SCENARIOS)}
    seeds = seeds.assign(o=seeds["scenario"].map(order)).sort_values(["o", "slot"])
    return [(r.run, f"{SUBGRAPHS}/{r.run}/union", f"{SIM_RUNS}/{r.run}",
             ["--block-flip", "unoriented"], r.scenario) for r in seeds.itertuples()]


def cases():
    """Every case as a dict: name, subgraph, truth, options, kind (regression / bench), scenario."""
    out = [dict(name=n, subgraph=s, truth=t, options=o, kind="regression", scenario=n) for n, s, t, o in REGRESSION]
    out += [dict(name=n, subgraph=s, truth=t, options=o, kind="bench", scenario=sc) for n, s, t, o, sc in bench_cases()]
    return out


def done(outdir, pool=False):
    need = ["summary.json", "evaluation.json"] + (["pool/pool.json", "pool/evaluation.tsv"] if pool else [])
    return all(os.path.exists(os.path.join(outdir, f)) for f in need)


def rck_run(case, outdir, extra, pool=False):
    from rck.ecdna_run import build_parser, run

    if done(outdir, pool):
        print("done already: " + outdir, flush=True)
        return True
    argv = ["--subgraph-dir", case["subgraph"], "--truth-dir", case["truth"], "--outdir", outdir] + COMMON + case["options"] + extra
    for attempt in range(1, RETRIES + 1):
        started = time.time()
        try:
            run(build_parser().parse_args(argv))
            print("{c}: {o} in {t:.1f} s".format(c=case["name"], o=os.path.basename(outdir), t=time.time() - started), flush=True)
            return True
        except Exception as error:  # noqa: BLE001 - licence hiccups: retry, then report
            print("attempt {a} failed: {c} {o}: {e!r}".format(a=attempt, c=case["name"], o=outdir, e=error), flush=True)
            if attempt < RETRIES:
                time.sleep(120)
    return False


def release_licence():
    """End this process's WLS session now; an exited process's session otherwise stays counted
    for a while, and short tasks then exceed the academic licence's 2 sessions."""
    try:
        import gurobipy

        gurobipy.disposeDefaultEnv()
    except Exception:  # noqa: BLE001 - nothing to release
        pass


def run_case(key):
    table = cases()
    case = table[int(key)] if str(key).isdigit() else next((c for c in table if c["name"] == key), None)
    if case is None:
        print("no case {k} ({n} cases; see run_case.py list)".format(k=key, n=len(table)))
        return 0
    if not os.path.exists(os.path.join(case["subgraph"], "edges.tsv")):
        print("SKIP {c}: no subgraph at {s}".format(c=case["name"], s=case["subgraph"]))
        return 1
    base = os.path.join(OUT, case["name"])
    ok = rck_run(case, os.path.join(base, "lambda_0"), ["--step-penalty", "0"])
    ok &= rck_run(case, os.path.join(base, "pool"), ["--step-penalty", POOL_LAMBDA, "--keep-pool",
                                                     "--pool-lambdas"] + [str(x) for x in LAMBDAS], pool=True)
    release_licence()
    with open(os.path.join(base, "case.json"), "w") as handle:
        json.dump(case, handle, indent=2)
    return 0 if ok else 1


def main():
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    if len(sys.argv) < 2 or sys.argv[1] not in ("list", "run"):
        sys.exit(__doc__)
    if sys.argv[1] == "list":
        for i, c in enumerate(cases()):
            print(i, c["name"], c["kind"], c["subgraph"], sep="\t")
        return 0
    return run_case(sys.argv[2])


if __name__ == "__main__":
    sys.exit(main())
