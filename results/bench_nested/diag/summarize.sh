#!/bin/bash
# Rank-1 answers of a diagnostic run: λ = 0 run, then the pool re-ranked over the bench λ grid.
cd ~/RCK
for r in "$@"; do
  echo "== $(basename $r)"
  ~/.conda/envs/gurobi-env/bin/python -m rck.ecdna.rerank $r/pool --lambdas 0 0.02 0.05 0.1 0.11 0.2 0.5 2>/dev/null \
    | cut -f1,3,4,5,8,9,11 | column -t -s$'\t'
  ~/.conda/envs/gurobi-env/bin/python - "$r/lambda_0" <<'PY'
import sys, pandas as pd
e = pd.read_csv(sys.argv[1] + "/evaluation.tsv", sep="\t"); e = e[e.decomposition_rank == 1]
d = pd.read_csv(sys.argv[1] + "/decompositions.tsv", sep="\t"); d = d[(d.decomposition_rank == 1)].drop_duplicates("species")
print("λ = 0 run: species", list(zip(d.species, d.weight, d.status)), "composites", dict(zip(e.truth, e.composite.round(3))), "mean %.4f" % e.composite.mean())
PY
done
