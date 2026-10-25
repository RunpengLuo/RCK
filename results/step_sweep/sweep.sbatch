#!/bin/bash
#SBATCH --job-name=step-sweep
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=01:30:00
#SBATCH --output=/n/fs/ugrad/ug28/ym3292/tmp/step_sweep/%j.out
set -euo pipefail
export PATH=/u/ym3292/.conda/envs/gurobi-env/bin:$PATH
export GRB_LICENSE_FILE=$HOME/gurobi.lic
cd /n/fs/ugrad/ug28/ym3292/tmp/step_sweep
python sweep.py /n/fs/ugrad/ug28/ym3292/tmp/step_sweep/out2
