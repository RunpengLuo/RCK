#!/bin/bash
#SBATCH --job-name=iwa-partB
#SBATCH --cpus-per-task=4
#
# bench_IW_allele Part B: one multi-clone run per array task (rck_partB.py run IDX; 72 tasks, variant-major).
# Submit with submit_partB.sh, which sets --array (with %1: Gurobi licence), --time, --mem and --output.
# A finished run is skipped, so the array can be resubmitted.
set -euo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
echo "== Part B task $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas; print('gurobi', gurobipy.gurobi.version())"
python rck_partB.py run "$SLURM_ARRAY_TASK_ID"
