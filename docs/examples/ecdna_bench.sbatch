#!/bin/bash
#SBATCH --job-name=ecbench
#SBATCH --array=0-32%2
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x_%A_%a.out
#
# RCK step of the ecDNA benchmark: one array task per config in configs/bench.
# Simulations (runs/bench) and subgraphs (subgraphs/bench) must already exist.
# %2: the academic Gurobi WLS licence allows 2 concurrent sessions.
# Finished runs (summary.json present) are skipped, so resubmitting only redoes failures.

source ~/.bashrc
conda activate gurobi-env
shopt -s nullglob

export GRB_LICENSE_FILE=$HOME/gurobi.lic
export TMPDIR=$HOME/.cache/rck-ecdna-tests
mkdir -p $TMPDIR

SIM=~/amplicon-allele-analysis/simulations
CFGS=($SIM/configs/bench/*.yaml)
CFG=${CFGS[$SLURM_ARRAY_TASK_ID]}
if [ -z "$CFG" ]; then echo "no config for task $SLURM_ARRAY_TASK_ID"; exit 1; fi
NAME=$(basename $CFG .yaml)
RUN=$SIM/runs/bench/$NAME
SUB=$SIM/subgraphs/bench/$NAME
RES=~/RCK/results/bench/$NAME
echo "== $NAME (task $SLURM_ARRAY_TASK_ID) on $(hostname)"

# inputs present?
if [ ! -d $RUN/obs/windows ]; then echo "SKIP $NAME: no $RUN/obs/windows (simulation missing or not copied)"; exit 0; fi
if [ ! -d $RUN/truth ]; then echo "WARNING $NAME: no truth/ folder; RCK will run but evaluation will fail"; fi

AMPS=($SUB/amplicon_*/ $SUB/union*/)
echo "$NAME: ${#AMPS[@]} subgraph(s)"
if [ ${#AMPS[@]} -eq 0 ]; then echo "SKIP $NAME: no subgraphs in $SUB"; exit 0; fi

cd ~/RCK
for A in "${AMPS[@]}"; do
  AMP=$(basename $A)
  for S in sim_normal sim_tumor1 sim_tumor2; do
    OUTD=$RES/$S/$AMP
    if [ -s $OUTD/summary.json ]; then echo "done already: $NAME $S $AMP"; continue; fi
    for TRY in 1 2 3 4 5; do
      python -m rck.ecdna_run --subgraph-dir $A --obs-dir $RUN/obs --sample $S --outdir $OUTD \
          --truth-dir $RUN --time-limit 1800 && break
      echo "attempt $TRY failed: $NAME $S $AMP; retrying in 60 s"
      sleep 60
    done
    [ -s $OUTD/summary.json ] || echo "RCK failed: $NAME $S $AMP"
  done
done
echo "== $NAME done"