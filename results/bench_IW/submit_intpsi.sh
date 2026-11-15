#!/bin/bash
# bench_IW_intpsi methods: four 24-task arrays, one per method (NEXT_STEPS.md step 5; needs Yubi's go-ahead).
#
#   ./submit_intpsi.sh            # print the sbatch commands only
#   ./submit_intpsi.sh --submit   # submit them
#
#   A rck_bulk  %1, 90 min                      max 1908 s in array 32267060, plus licence waits
#   B rck_mc    %1, 15 min                      max 367 s (licence waits included), solves <= 4 s
#   C coral     %1, 90 min, aftercorr:A         max 2002 s, plus licence waits
#   D ce        %1, 60 min, aftercorr:A         max 374 s (CE's own limit is 7200 s per MILP, hard-coded)
#
# aftercorr: case i's coral / ce task starts only after case i's rck_bulk task succeeded (a failed bulk
# task leaves them pending with DependencyNeverSatisfied; cancel those). Memory 2G (peak 754 MB in
# 32267060). The Gurobi WLS licence has 2 sessions; A and B hold at most one each, C and D join as bulk
# tasks finish, so up to 4 may ask at once: a task waits for the licence up to 15 x 2 min
# (bench_iw.py LICENCE_RETRIES) and fails after that; resubmitting redoes only unfinished runs.
set -euo pipefail
cd "$(dirname "$0")"
L=/n/fs/yubi-iw/runs/bench_IW_intpsi/logs/slurm
N=$(python bench_iw.py --set iw_intpsi list --method rck_bulk | wc -l)
[ "$N" -eq 24 ] || { echo "expected 24 cases, found $N"; exit 1; }
ARR="0-$((N - 1))%1"
common() { echo "--parsable --array=$ARR --mem=2G --output=$L/%x_%A_%a.out --export=ALL,SET=iw_intpsi,METHOD=$1 --job-name=iwi-$1"; }

if [ "${1:-}" != "--submit" ]; then
  echo "mkdir -p $L"
  echo "A=\$(sbatch $(common rck_bulk) --time=01:30:00 bench_iw_method.sbatch)"
  echo "B=\$(sbatch $(common rck_mc) --time=00:15:00 bench_iw_method.sbatch)"
  echo "C=\$(sbatch $(common coral) --time=01:30:00 --dependency=aftercorr:\$A bench_iw_method.sbatch)"
  echo "D=\$(sbatch $(common ce) --time=01:00:00 --dependency=aftercorr:\$A bench_iw_method.sbatch)"
  exit 0
fi
mkdir -p "$L"
A=$(sbatch $(common rck_bulk) --time=01:30:00 bench_iw_method.sbatch)
B=$(sbatch $(common rck_mc) --time=00:15:00 bench_iw_method.sbatch)
C=$(sbatch $(common coral) --time=01:30:00 --dependency=aftercorr:"$A" bench_iw_method.sbatch)
D=$(sbatch $(common ce) --time=01:00:00 --dependency=aftercorr:"$A" bench_iw_method.sbatch)
echo "rck_bulk $A  rck_mc $B  coral $C  ce $D" | tee -a "$L/submitted.txt"
