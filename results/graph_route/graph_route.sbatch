#!/bin/bash
#SBATCH --job-name=graph-route
#SBATCH --array=0-224%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%A_%a.out
#
# Graph route (Task 2) head-to-head, priority 1 a: one run per line of graph_route.list
# (python graph_route.py list): 36 cases (demo5-demo8 + 32 bench) x graphs {truth_rp: 39 truth
# amplicon graphs with reads and path constraints, inferred: 36 exported subgraphs} x methods
# {CoRAL 3.1.1 cycle, RCK-ecDNA on aa_import --imputed none, CycleExtractor 1.0.0}.
# Lines 0-149: CoRAL and RCK-ecDNA, per graph; lines 150-224: CycleExtractor. Prepare first
# (offline), then submit from here:
#
#   cd ~/RCK/results/graph_route && python graph_route.py prepare     # done 2026-11-04; rerun only if inputs change
#   L=/n/fs/scratch/$USER/graph_route/h2h/logs; mkdir -p $L; sbatch --output=$L/%x_%A_%a.out graph_route.sbatch   # job 32104590, 2026-11-04
#   python graph_route.py status                                       # offline, afterwards
#
# %2: CoRAL, CycleExtractor (--solver gurobi) and RCK-ecDNA all use the academic Gurobi WLS licence
# (2 sessions). Finished runs are skipped, so resubmitting redoes only what is missing. Memory:
# RCK-ecDNA peaked at 366 MB in the round trip; CoRAL and CE graphs are as small. Time: RCK-ecDNA
# --time-limit 1800 per MILP (bench pool runs: max 1896 s), CoRAL --solver-time-limit 1800 per MILP,
# CE's hard-coded 7200 s per MILP (HiGHS took < 4 s on these graphs).

set -euo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
unset BENCH_NESTED_ROOT BENCH_NESTED_OUT
echo "== task $SLURM_ARRAY_TASK_ID on $(hostname)"
python -c "import gurobipy, pandas; print('gurobi', gurobipy.gurobi.version())"
python graph_route.py run "$SLURM_ARRAY_TASK_ID"
