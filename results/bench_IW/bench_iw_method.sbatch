#!/bin/bash
#SBATCH --job-name=bench-iw-m
#SBATCH --cpus-per-task=4
#
# One method of one bench_IW set, one case per array task (bench_iw.py run --method M IDX, IDX = case):
# the one-array-per-method layout, so the decomposition-only baselines (coral, ce) can depend task by
# task on rck_bulk (--dependency=aftercorr:<rck_bulk array>). Submit with submit_intpsi.sh, which sets
# --array, --time, --mem, --output and the dependency, and exports SET and METHOD:
#
#   SET=iw_intpsi METHOD=coral sbatch --array=0-23%1 --dependency=aftercorr:<bulk job> ... bench_iw_method.sbatch
#
# A finished run is skipped, so an array can be resubmitted. A coral / ce task whose case has no finished
# rck_bulk run exits 2 without running (bench_iw.py bulk_graph). The older bench_iw.sbatch (one 96-task
# array, array 32267060) is kept as it ran.

set -euo pipefail
: "${SET:?export SET (iw or iw_intpsi)}"
: "${METHOD:?export METHOD (rck_bulk, rck_mc, coral or ce)}"
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
echo "== $SET $METHOD case $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas; print('gurobi', gurobipy.gurobi.version())"
python bench_iw.py --set "$SET" run --method "$METHOD" "$SLURM_ARRAY_TASK_ID"
