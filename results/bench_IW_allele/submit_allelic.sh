#!/bin/bash
# PLAN.md item 1: four 12-task arrays on the allelic-sibling cases, one per method (each needs Yubi's go-ahead).
#
#   ./submit_allelic.sh            # print the sbatch commands only
#   ./submit_allelic.sh --submit   # submit them
#
# Sizes from bench_IW_intpsi's arrays (results/bench_IW/submit_intpsi.sh, same settings and similar cases):
#   A rck_bulk  %1, 90 min                     (max 1908 s there, plus licence waits)
#   B rck_mc    %1, 15 min                     (solves <= 4 s there)
#   C coral     %1, 90 min, aftercorr:A        (max 2002 s)
#   D ce        %1, 60 min, aftercorr:A        (max 374 s)
# Memory 2G (peak 796 MB there). The Gurobi WLS licence has 2 sessions; A and B hold at most one each, C and D
# join as bulk tasks finish, so a task may wait for the licence (bench_iw LICENCE_RETRIES x 2 min) and fails after
# that; resubmitting redoes only unfinished runs. Only submit when no other array (e.g. Part B) holds sessions.
set -euo pipefail
cd "$(dirname "$0")"
L=/n/fs/yubi-iw/runs/bench_IW_allelic/logs/slurm
N=$(python allelic_bench.py list --method rck_bulk | wc -l)
[ "$N" -eq 12 ] || { echo "expected 12 cases, found $N"; exit 1; }
ARR="0-$((N - 1))%1"
common() { echo "--parsable --array=$ARR --mem=2G --output=$L/%x_%A_%a.out --export=ALL,METHOD=$1 --job-name=iwa-$1"; }
if [ "${1:-}" != "--submit" ]; then
  echo "mkdir -p $L"
  echo "A=\$(sbatch $(common rck_bulk) --time=01:30:00 allelic_method.sbatch)"
  echo "B=\$(sbatch $(common rck_mc) --time=00:15:00 allelic_method.sbatch)"
  echo "C=\$(sbatch $(common coral) --time=01:30:00 --dependency=aftercorr:\$A allelic_method.sbatch)"
  echo "D=\$(sbatch $(common ce) --time=01:00:00 --dependency=aftercorr:\$A allelic_method.sbatch)"
  exit 0
fi
mkdir -p "$L"
A=$(sbatch $(common rck_bulk) --time=01:30:00 allelic_method.sbatch)
B=$(sbatch $(common rck_mc) --time=00:15:00 allelic_method.sbatch)
C=$(sbatch $(common coral) --time=01:30:00 --dependency=aftercorr:"$A" allelic_method.sbatch)
D=$(sbatch $(common ce) --time=01:00:00 --dependency=aftercorr:"$A" allelic_method.sbatch)
echo "rck_bulk $A  rck_mc $B  coral $C  ce $D" | tee -a "$L/submitted.txt"
