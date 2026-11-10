"""Sweep --step-penalty over the demo runs; write sweep.tsv (one row per case x lambda)."""
import json, os, sys, time, logging
import pandas as pd
sys.path.insert(0, "/n/fs/ugrad/ug28/ym3292/RCK")
from rck.ecdna_run import build_parser, run

AAA = "/n/fs/ugrad/ug28/ym3292/amplicon-allele-analysis/simulations/runs"
DATA = "/n/fs/ugrad/ug28/ym3292/RCK/tests/ecdna/data"
OUT = sys.argv[1]
CASES = {
    "demo3": (f"{DATA}/amplicon_0", f"{DATA}/demo3_truth", ["--depth-per-copy", "15", "--junction-reads-per-copy", "6"]),
    "demo4_amp0": (f"{DATA}/demo4/amplicon_0", f"{DATA}/demo4_truth", ["--depth-per-copy", "13.5", "--truth-sample", "tumor1", "--min-species-weight", "2"]),
    "demo4_amp1": (f"{DATA}/demo4/amplicon_1", f"{DATA}/demo4_truth", ["--depth-per-copy", "13.5", "--truth-sample", "tumor1", "--min-species-weight", "2"]),
    "demo5": (f"{AAA}/demo5_cycle1_ecdna/subgraph/amplicon_0", f"{AAA}/demo5_cycle1_ecdna", []),
    "demo6_union": (f"{AAA}/demo6_cycle2_ecdna/subgraph_union/union", f"{AAA}/demo6_cycle2_ecdna", []),
    "demo7_blockflip": (f"{AAA}/demo7_nestedcycles_ecdna/subgraph_union_blocks/union", f"{AAA}/demo7_nestedcycles_ecdna", ["--block-flip", "unoriented"]),
    "demo7_totalonly": (f"{AAA}/demo7_nestedcycles_ecdna/subgraph_union/union", f"{AAA}/demo7_nestedcycles_ecdna", ["--min-snps-allele", "100000000"]),
    "demo8": (f"{AAA}/demo8_nestedcycles_ecdna/subgraph_union/union", f"{AAA}/demo8_nestedcycles_ecdna", []),
}
LAMBDAS = ["0", "0.02", "0.05", "0.1", "0.2"]
logging.basicConfig(level=logging.WARNING)
rows = []
for case, (sub, truth, extra) in CASES.items():
    for lam in LAMBDAS:
        out = os.path.join(OUT, case, "lambda_" + lam)
        t = time.time()
        argv = ["--subgraph-dir", sub, "--truth-dir", truth, "--outdir", out, "--no-plots", "--num-alternatives", "5",
                "--threads", os.environ.get("SLURM_CPUS_PER_TASK", "4"), "--step-penalty", lam] + extra
        row = {"case": case, "step_penalty": float(lam)}
        try:
            res = run(build_parser().parse_args(argv))
            top = res["decompositions"][0]
            amb = top.ambiguity or {}
            m = pd.read_csv(os.path.join(out, "recovery_matrix.tsv"), sep="\t")
            m = m[m["in_scope"].astype(bool)]
            row.update({"n_species": sum(1 for s in top.species if not s.low_weight),
                        "weights": ",".join("%g" % s.weight for s in top.species if not s.low_weight),
                        "steps": top.steps, "residual": top.residual,
                        "composites": ",".join("%s=%.3f" % (r.truth, r.composite) for r in m.itertuples()),
                        "mean_composite": round(float(m["composite"].mean()), 4),
                        "ambiguous": amb.get("ambiguous"), "decided_by": amb.get("decided_by"), "gap": amb.get("gap"),
                        "n_alternatives": amb.get("n_alternatives"), "runtime_s": round(time.time() - t, 1)})
        except Exception as e:  # keep sweeping
            row["error"] = repr(e)
        rows.append(row)
        print(row, flush=True)
        pd.DataFrame(rows).to_csv(os.path.join(OUT, "sweep.tsv"), sep="\t", index=False)
