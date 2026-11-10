#!/bin/bash
#SBATCH --job-name=graph-route-rt
#SBATCH --array=0-143%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=01:15:00
#SBATCH --output=logs/%x_%A_%a.out
#
# Round-trip test of the graph route (Task 2 step 3), Gurobi part: one RCK-ecDNA total-only run
# per line of roundtrip.list (python roundtrip.py list): 36 cases (demo5-demo8 + 32 bench) x 4 arms
# (original, rt_comment, rt_none_profile, rt_none). Prepare first (offline), then submit from here:
#
#   cd ~/RCK/results/graph_route && python roundtrip.py prepare      # done 2026-11-03; rerun only if the code changes
#   L=/n/fs/scratch/$USER/graph_route/roundtrip/logs; mkdir -p $L; sbatch --output=$L/%x_%A_%a.out roundtrip.sbatch   # job 32096382, 2026-11-03
#   python roundtrip.py compare                                       # after the array, offline
#
# %2: the academic Gurobi WLS licence allows 2 concurrent sessions. Finished runs are skipped, so
# resubmitting redoes only what is missing. Memory as bench_nested.sbatch (bench peak 187 MB,
# demo7 0.8 GB); time: two MILPs at --time-limit 1800 each, bench runs took 1-358 s.

set -euo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
unset BENCH_NESTED_ROOT BENCH_NESTED_OUT
echo "== task $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas; print('gurobi', gurobipy.gurobi.version())"
python roundtrip.py run "$SLURM_ARRAY_TASK_ID"
