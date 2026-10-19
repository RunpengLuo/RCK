#!/bin/bash
# Run the rck-ecdna test suite, including both demo3 tests, on ionic.
#
#   cd ~/RCK && sbatch docs/examples/ecdna_demo3_tests.sbatch
#
# Needs: conda env gurobi-env (numpy, pandas, gurobi), ~/gurobi.lic (academic WLS), and the
# amplicon-allele-analysis checkout below. Test outputs go to ~/.cache/rck-ecdna-tests rather
# than the node's /tmp, which is shared and often full.
#SBATCH --job-name=rck-ecdna-demo3
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=rck-ecdna-demo3-%j.out

set -euo pipefail
AAA=${AAA:-$HOME/amplicon-allele-analysis}
RCK=${RCK:-$HOME/RCK}

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gurobi-env
export GRB_LICENSE_FILE="$HOME/gurobi.lic"
export RCK_ECDNA_TEST_TMP="$HOME/.cache/rck-ecdna-tests"
export TMPDIR="$RCK_ECDNA_TEST_TMP"
export RCK_ECDNA_DEMO3_OBS="$AAA/simulations/runs/demo3_clone1_ecdna/obs"
mkdir -p "$RCK_ECDNA_TEST_TMP"

cd "$RCK"
echo "commit: $(git log --oneline -1)"
python -c "import gurobipy as gp; m = gp.Model(); print('gurobi', gp.gurobi.version()); m.dispose()"
python -m unittest discover -s tests/ecdna -t . -v
rm -rf "$RCK_ECDNA_TEST_TMP"/tmp*
