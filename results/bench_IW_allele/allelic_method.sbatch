#!/bin/bash
#SBATCH --job-name=iwa-allelic
#SBATCH --cpus-per-task=4
#
# PLAN.md item 1: one method on the allelic-sibling cases, one case per array task
# (allelic_bench.py run --method $METHOD IDX). Submit with submit_allelic.sh (--array, --time, --mem, --output,
# the aftercorr dependency of coral / ce on rck_bulk). A finished run is skipped, so an array can be resubmitted.
set -euo pipefail
: "${METHOD:?export METHOD (rck_bulk, rck_mc, coral or ce)}"
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
echo "== iw_allelic $METHOD case $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas; print('gurobi', gurobipy.gurobi.version())"
python allelic_bench.py run --method "$METHOD" "$SLURM_ARRAY_TASK_ID"
