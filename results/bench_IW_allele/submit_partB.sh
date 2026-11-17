#!/bin/bash
# bench_IW_allele Part B (PLAN.md): the multi-clone route on RCK's degraded inputs. NEEDS YUBI'S GO-AHEAD.
#
#   ./submit_partB.sh                    # print the sbatch command only
#   ./submit_partB.sh --pilot            # print the 3-task pilot (one case, each variant)
#   ./submit_partB.sh --submit [--pilot] # submit
#
# 72 tasks (3 variants x 24 cases, variant-major: 0-23 chil, 24-47 chil_avg, 48-71 chil_err), %1 for the
# Gurobi WLS licence (2 sessions). Reference rck_mc runs: solves <= 5 s, <= 367 s per task with licence
# waits; Part B allows 600 s per MILP, so 30 min per task. Memory 2G (as bench_IW_intpsi's rck_mc).
# Pilot: tasks 21, 45, 69 = hard_S4_s16003 (the largest DIAG, 214 segments) in each variant.
set -euo pipefail
cd "$(dirname "$0")"
L=/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/logs/slurm
N=$(python rck_partB.py list | wc -l)
[ "$N" -eq 72 ] || { echo "expected 72 tasks, found $N"; exit 1; }
PILOT=$(python rck_partB.py list | awk -F'\t' '$3=="hard_S4_s16003"{printf "%s%s", s, $1; s=","}')
ARR="0-$((N - 1))%1"
NAME=iwa-partB
for a in "$@"; do [ "$a" = "--pilot" ] && ARR="$PILOT%1" && NAME=iwa-partB-pilot; done
CMD="sbatch --parsable --array=$ARR --mem=2G --time=00:30:00 --output=$L/%x_%A_%a.out --job-name=$NAME rck_partB.sbatch"
if [ "${1:-}" != "--submit" ]; then
  echo "mkdir -p $L"
  echo "$CMD"
  exit 0
fi
mkdir -p "$L"
J=$($CMD)
echo "$NAME $J array $ARR" | tee -a "$L/submitted.txt"
