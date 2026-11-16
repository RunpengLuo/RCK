#!/bin/bash
#SBATCH --job-name=iwe-pilot
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=06:00:00
#
# Pilot of the whole bench_IW_e2e chain for one case per array task (IDX = case index of
# `bench_iw_e2e.py list`; the pilots are 12 = hard_S1_s14001 and 21 = hard_S4_s16003):
#   reads (emit_reads + align_cnv; a finished step is skipped) -> check_reads -> inputs (obs_from_bam +
#   build_subgraphs_e2e) -> coral_e2e -> ce_e2e -> rck_bulk_e2e -> rck_mc_e2e -> score --case IDX
# A method that fails does not stop the others. Timings per step go to the Slurm log (sizes the arrays).
#   sbatch --array=12,21 --output=/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/logs/slurm/%x_%A_%a.out pilot_e2e.sbatch
set -uo pipefail
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd "$SLURM_SUBMIT_DIR"
IDX=$SLURM_ARRAY_TASK_ID
CASE=$(python bench_iw_e2e.py list | awk -v i="$IDX" '$1 == i {print $2}')
SIM=${AAA:-$HOME/amplicon-allele-analysis-cnp}/simulations
E2E=/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e
step() { local t0=$(date +%s); echo "== $(date +%T) $*"; "$@"; local rc=$?; echo "== rc $rc, $(( $(date +%s) - t0 )) s: $*"; return $rc; }
echo "== pilot $CASE (case $IDX) on $(hostname), $SLURM_CPUS_PER_TASK cpus"
step python bench_iw_e2e.py reads "$IDX" || exit 1
step /n/fs/yubi-iw/conda/envs/e2e-align/bin/python "$SIM/e2e/check_reads.py" --run "/n/fs/yubi-iw/runs/bench_IW_intpsi/runs/$CASE" \
  --config "$SIM/configs/bench_IW_intpsi/$CASE.yaml" --reads "$E2E/reads/$CASE"
step python bench_iw_e2e.py inputs "$IDX"
for m in coral_e2e ce_e2e rck_bulk_e2e rck_mc_e2e; do
  step python bench_iw_e2e.py run --method "$m" "$IDX"
done
step python bench_iw_e2e.py score --case "$IDX" --jobs 4
python bench_iw_e2e.py status
