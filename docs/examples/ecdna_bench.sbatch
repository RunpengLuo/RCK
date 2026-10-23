#!/bin/bash
#SBATCH --job-name=ecbench
#SBATCH --array=0-32
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=%x_%A_%a.out
source ~/.bashrc
conda activate gurobi-env
export GRB_LICENSE_FILE=$HOME/gurobi.lic TMPDIR=$HOME/.cache/rck-ecdna-tests
mkdir -p $TMPDIR
SIM=~/amplicon-allele-analysis/simulations
REF=$SIM/../data/COLO320/COLO320DM_0911/GRCh38_reference
CFGS=($SIM/configs/bench/*.yaml)
CFG=${CFGS[$SLURM_ARRAY_TASK_ID]}
NAME=$(basename $CFG .yaml)
RUN=$SIM/runs/bench/$NAME
SUB=$SIM/subgraphs/bench/$NAME
RES=~/RCK/results/bench/$NAME
cd ~/RCK
for A in $SUB/amplicon_*/ $SUB/union*/; do
  for S in sim_normal sim_tumor1 sim_tumor2; do
    OUTD=$RES/$S/$(basename $A)
    [ -s $OUTD/summary.json ] && continue          # already done
    for TRY in 1 2 3 4 5; do
      python -m rck.ecdna_run --subgraph-dir $A --obs-dir $RUN/obs --sample $S --outdir $OUTD --truth-dir $RUN --time-limit 1800 && break
      echo "attempt $TRY failed: $NAME $S $(basename $A); retrying in 60 s"; sleep 60
    done
    [ -s $OUTD/summary.json ] || echo "RCK failed: $NAME $S $(basename $A)"
  done
done