#!/bin/bash
#SBATCH --job-name=rck-valid
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --output=logs/%x_%A_%a.out
# Validation run of a bench_nested case with extra rck-ecdna options, using the case's exact bench
# configuration (subgraph, truth and options from run_case.py, so regression cases keep theirs; unlike
# rck_diag.sbatch, which adds --block-flip unoriented to every case). One --keep-pool run at
# VALID_POOL_LAMBDA (default 0.2), so rank 1 is that λ's answer; no λ = 0 run.
# As an array, task N takes line N+1 of $LIST: NAME OUTDIR [RCK OPTIONS], NAME as in `run_case.py list`.
# Submit from a folder with a logs/ subfolder. The licence allows 2 sessions: use --array=...%2.
set -euo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
export AAA=${AAA:-$HOME/amplicon-allele-analysis}
cd ~/RCK/results/bench_nested
set -- $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST")
python - "$@" <<'PY'
import os, sys
import run_case as rc
name, out = sys.argv[1:3]
case = [c for c in rc.cases() if c["name"] == name]
if len(case) != 1:
    sys.exit("unknown case " + name)
case = dict(case[0], options=case[0]["options"] + sys.argv[3:])
lam = os.environ.get("VALID_POOL_LAMBDA", "0.2")
ok = rc.rck_run(case, os.path.join(out, "pool"), ["--step-penalty", lam, "--keep-pool", "--pool-lambdas"] + [str(x) for x in rc.LAMBDAS], pool=True)
rc.release_licence()
sys.exit(0 if ok else 1)
PY
