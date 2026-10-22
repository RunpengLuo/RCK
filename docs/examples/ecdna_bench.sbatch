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
cd $SIM
python src/simulation.py -c $CFG -g $REF/hg38.sizes -b $REF/hg38.regions.bed -o $RUN --no-plots
python subgraph_construct.py --sim-dir $RUN --out-dir $SUB --no-plot
python verify_subgraph.py --run-dir $RUN --out-dir $SUB
cd ~/RCK
for A in $SUB/amplicon_*/; do
  for S in sim_normal sim_tumor1 sim_tumor2; do
    python -m rck.ecdna_run --subgraph-dir $A --obs-dir $RUN/obs --sample $S --outdir $RES/$S/$(basename $A) --truth-dir $RUN --time-limit 1800
  done
done