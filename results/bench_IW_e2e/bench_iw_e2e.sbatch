#!/bin/bash
#SBATCH --job-name=iw-e2e
#
# One case of one bench_IW_e2e step per array task (bench_iw_e2e.py; IDX = case index):
#   STEP=reads                  simulated reads + BAMs + CNVkit (no Gurobi)
#   STEP=run METHOD=coral_e2e   coral seed + reconstruct (Gurobi)
#   STEP=run METHOD=ce_e2e      CycleExtractor on coral_e2e's graphs (Gurobi)
# Submit with submit_e2e.sh (sets --array, --time, --mem, --cpus-per-task, --output, dependencies).
# A finished task is skipped, so an array can be resubmitted.
set -euo pipefail
: "${STEP:?export STEP (reads or run)}"
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
echo "== $STEP ${METHOD:-} case $SLURM_ARRAY_TASK_ID on $(hostname)"
if [ "$STEP" = reads ]; then
  python bench_iw_e2e.py reads "$SLURM_ARRAY_TASK_ID"
else
  : "${METHOD:?export METHOD (coral_e2e or ce_e2e)}"
  python bench_iw_e2e.py run --method "$METHOD" "$SLURM_ARRAY_TASK_ID"
fi
