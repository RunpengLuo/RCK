#!/bin/bash
#SBATCH --job-name=rck-diag
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%j.out
# Submit from a folder with a logs/ subfolder, e.g. /n/fs/scratch/$USER/bench_nested_v2/diag.
# Diagnostic RCK run (λ = 0 and a --keep-pool run at 0.1) on one rebuilt subgraph, with the
# options of RCK results/bench_nested/run_case.py. Usage: sbatch rck_diag.sbatch NAME SUBGRAPH TRUTH OUTDIR [RCK OPTIONS]
# As an array (--array=0-N), each task takes line N+1 of $LIST instead of the arguments.
# DIAG_POOL_LAMBDA sets the pool run's --step-penalty (default run_case.py's 0.1; any λ > 0 gives the
# same pool, and rank 1 is that λ's), and DIAG_LAMBDA0=0 skips the λ = 0 run.
set -euo pipefail
export PATH=$HOME/.conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd ~/RCK/results/bench_nested
if [ -n "${LIST:-}" ]; then set -- $(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$LIST"); fi
python - "$@" <<'PY'
import os, sys
import run_case as rc
name, sub, truth, out = sys.argv[1:5]
case = dict(name=name, subgraph=sub, truth=truth, options=["--block-flip", "unoriented"] + sys.argv[5:])
ok = True
if os.environ.get("DIAG_LAMBDA0", "1") != "0":
    ok &= rc.rck_run(case, os.path.join(out, "lambda_0"), ["--step-penalty", "0"])
pool_lambda = os.environ.get("DIAG_POOL_LAMBDA", rc.POOL_LAMBDA)
ok &= rc.rck_run(case, os.path.join(out, "pool"), ["--step-penalty", pool_lambda, "--keep-pool", "--pool-lambdas"] + [str(x) for x in rc.LAMBDAS], pool=True)
rc.release_licence()
sys.exit(0 if ok else 1)
PY
