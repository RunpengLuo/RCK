#!/bin/bash
#SBATCH --job-name=bench-iw
#SBATCH --array=0-95%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%A_%a.out
#
# bench_IW: one (case, method) run per array task, in the order of `python bench_iw.py list`:
# 24 accepted cases x {rck_bulk, rck_mc, coral, ce} = 96 tasks. Submit from here:
#
#   cd ~/RCK-cnp/results/bench_IW
#   L=/n/fs/yubi-iw/runs/bench_IW/logs/slurm; mkdir -p $L; sbatch --output=$L/%x_%A_%a.out bench_iw.sbatch
#   python bench_iw.py status          # offline, afterwards; then score and collect
#
# %2: every method uses the academic Gurobi WLS licence (2 sessions). Finished runs are skipped, so
# resubmitting redoes only what is missing. Memory: rck-ecdna peaked at 366 MB on bench_nested
# graphs, CoRAL and CE are as small; 4G is the graph route's request. Time: rck_bulk --time-limit
# 1800 per MILP (bench_nested pool runs: max 1896 s), CoRAL 1800 s per MILP, CE 7200 s per MILP
# (hard-coded), rck_mc 60 s.

set -euo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
echo "== task $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas; print('gurobi', gurobipy.gurobi.version())"
python bench_iw.py run "$SLURM_ARRAY_TASK_ID"
