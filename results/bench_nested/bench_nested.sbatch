#!/bin/bash
#SBATCH --job-name=bench-nested
#SBATCH --array=0-21%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:30:00
#SBATCH --output=logs/%x_%A_%a.out
#
# RCK step of bench_nested: one array task per case (python run_case.py list), each a λ = 0 run
# and a --keep-pool run at λ = 0.1 (see README.md). Submit from this folder, after the
# simulation array (amplicon-allele-analysis simulations/bench_nested.sbatch) has finished:
#
#   cd ~/RCK/results/bench_nested && mkdir -p logs && sbatch bench_nested.sbatch
#   sbatch --array=0,6,7 bench_nested.sbatch          # demo3, demo7_totalonly, demo8 only: the
#                                                      # re-ranking check (python collect.py --validate)
#
# %2: the academic Gurobi WLS licence allows 2 concurrent sessions. Finished runs are skipped,
# so resubmitting redoes only what is missing. 0-21 = 8 regression cases + 7 scenarios x 2 seeds;
# tasks past the last accepted seed exit with "no case". Memory: the 40-run step sweep peaked at
# 187 MB, demo7 subgraph + RCK at 0.8 GB.

set -euo pipefail
# the conda env's bin on PATH, as results/step_sweep/sweep.sbatch: `conda activate` needs conda's
# shell setup, which ~/.bashrc does not give a batch job on a compute node
export PATH=$HOME/.conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
export AAA=${AAA:-$HOME/amplicon-allele-analysis}
cd "$SLURM_SUBMIT_DIR"
echo "== task $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas, scipy; print('gurobi', gurobipy.gurobi.version())"
python run_case.py run "$SLURM_ARRAY_TASK_ID"
