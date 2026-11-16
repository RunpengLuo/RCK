#!/bin/bash
# bench_IW_e2e arrays (needs Yubi's go-ahead; prints the commands unless --submit):
# NOT SIZED YET: the limits below predate the sentinel windows (reads ~5x longer) and there are no arrays yet
# for inputs / rck_bulk_e2e / rck_mc_e2e. Resize from the pilot array 32318772 first (NEXT_STEPS.md, "Resume here" 3).
#
#   ./submit_e2e.sh            # print
#   ./submit_e2e.sh --submit   # submit
#
#   R reads      24 tasks, %8, 8 cpus, 16G, 1 h    no Gurobi (NanoSim ~4 min, minimap2 + CNVkit; sized from the pilots)
#   C coral_e2e  24 tasks, %1, 8 cpus, 8G, 3 h     aftercorr:R; Gurobi, --solver-time-limit 1800 per MILP
#   E ce_e2e     24 tasks, %1, 4 cpus, 4G, 2 h     aftercorr:C; Gurobi
#
# The Gurobi WLS licence has 2 sessions: C and E hold at most one each (%1). A task waits for the
# licence up to 15 x 2 min (bench_iw.LICENCE_RETRIES).
set -euo pipefail
cd "$(dirname "$0")"
L=/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/logs/slurm
N=$(python bench_iw_e2e.py list | wc -l)
[ "$N" -eq 24 ] || { echo "expected 24 cases, found $N"; exit 1; }
common() { echo "--parsable --array=0-$((N - 1))%$1 --output=$L/%x_%A_%a.out --job-name=iwe-$2"; }
R_ARGS="$(common 8 reads) --cpus-per-task=8 --mem=16G --time=01:00:00 --export=ALL,STEP=reads"
C_ARGS="$(common 1 coral) --cpus-per-task=8 --mem=8G --time=03:00:00 --export=ALL,STEP=run,METHOD=coral_e2e"
E_ARGS="$(common 1 ce) --cpus-per-task=4 --mem=4G --time=02:00:00 --export=ALL,STEP=run,METHOD=ce_e2e"
if [ "${1:-}" != "--submit" ]; then
  echo "mkdir -p $L"
  echo "R=\$(sbatch $R_ARGS bench_iw_e2e.sbatch)"
  echo "C=\$(sbatch $C_ARGS --dependency=aftercorr:\$R bench_iw_e2e.sbatch)"
  echo "E=\$(sbatch $E_ARGS --dependency=aftercorr:\$C bench_iw_e2e.sbatch)"
  exit 0
fi
mkdir -p "$L"
R=$(sbatch $R_ARGS bench_iw_e2e.sbatch)
C=$(sbatch $C_ARGS --dependency=aftercorr:"$R" bench_iw_e2e.sbatch)
E=$(sbatch $E_ARGS --dependency=aftercorr:"$C" bench_iw_e2e.sbatch)
echo "reads $R  coral_e2e $C  ce_e2e $E" | tee -a "$L/submitted.txt"
