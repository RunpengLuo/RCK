#!/bin/bash
# Run rck-ecdna on every amplicon of one subgraph set, then the ecDNA test suite.
#
#   sbatch docs/examples/ecdna_ionic.sbatch
#
# Edit the paths below. Requires the conda environment from docs/ecDNA.md section 7 and a
# Gurobi license usable on compute nodes (academic WLS gurobi.lic, or a cluster module).
#SBATCH --job-name=rck-ecdna
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=rck-ecdna-%j.out

set -euo pipefail

AAA=$HOME/amplicon-allele-analysis          # amplicon-allele-analysis checkout
RCK=$HOME/RCK                               # this repository (yubi_dev)
RUN=$AAA/simulations/runs/demo3_clone1_ecdna
SUBGRAPHS=$AAA/subgraphs/demo3
OUT=$HOME/rck_ecdna_out/demo3

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rck-ecdna
# export GRB_LICENSE_FILE=$HOME/gurobi.lic   # if the license is not in a default location

python -c "import gurobipy; m = gurobipy.Model(); print('gurobi', gurobipy.gurobi.version())"

for amp in "$SUBGRAPHS"/amplicon_*; do
    name=$(basename "$amp")
    rck-ecdna \
        --subgraph-dir "$amp" \
        --obs-dir "$RUN/obs" \
        --truth-dir "$RUN" \
        --outdir "$OUT/$name" \
        --threads "${SLURM_CPUS_PER_TASK:-4}" \
        --time-limit 1800
done

cd "$RCK"
RCK_ECDNA_DEMO3_OBS=$RUN/obs python -m unittest discover -s tests/ecdna -t . -v
