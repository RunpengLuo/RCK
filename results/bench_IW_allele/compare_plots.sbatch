#!/bin/bash
#SBATCH --job-name=iwa-compare
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#
# RCK MC vs CoRAL topology comparison (make_compare_plots.py): one case per array task, no Gurobi.
# Submit from this folder:
#   sbatch --array=0-7 --output=/n/fs/yubi-iw/runs/bench_IW_allelic/logs/slurm/%x_%A_%a.out compare_plots.sbatch
set -euo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
cd "$SLURM_SUBMIT_DIR"
echo "== compare case $SLURM_ARRAY_TASK_ID on $(hostname)"
python make_compare_plots.py --index "$SLURM_ARRAY_TASK_ID"
