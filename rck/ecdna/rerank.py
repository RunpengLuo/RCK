"""Re-rank a kept decomposition pool at other ``--step-penalty`` values, offline.

``rck-ecdna --step-penalty λ --keep-pool`` (any λ > 0) writes ``pool/``: every solution of the
widened decomposition pool, before the λ-dependent pruning, annotated, refitted and evaluated
like a ranked decomposition, plus the same-fit alternatives of each solution that is rank 1 at
some λ of ``--pool-lambdas``. The decomposition MILP and its pool do not depend on λ (the
penalty acts in the ranking only, and every λ > 0 widens the pool to the same
``--ambiguity-window``), so for any λ > 0 this gives exactly what a run at that λ ranks first:

* the pruning (:func:`rck.ecdna.decompose.pool_keeps`) and the ranking
  (:func:`rck.ecdna.decompose.ranking_key`) are the functions the run itself uses;
* the composites come from ``pool/evaluation.tsv``; a decomposition's scores depend on that
  decomposition alone;
* ``ambiguous``, ``decided_by`` and ``gap`` come from :func:`rck.ecdna.decompose.ambiguity_summary`
  on the recorded same-fit alternatives. They are missing for a λ whose rank-1 solution was not
  rank 1 at any λ of ``--pool-lambdas``.

λ = 0 is emulated: a λ = 0 run solves a narrower pool (``--pool-gap`` relative, fewer
solutions), here approximated by the solutions within ``--pool-gap`` of the best. Solutions that
tie in the ranking may then come out in another order (``source`` says ``pool, λ=0 emulated``).

    python -m rck.ecdna.rerank RUN_DIR --lambdas 0 0.02 0.05 0.1 0.2 0.5 -o rerank.tsv
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Sequence

import pandas as pd

from rck.ecdna.decompose import DecompositionParams, ambiguity_summary, rank_pool

COLUMNS = ["step_penalty", "pool_index", "n_species", "weights", "steps", "residual", "objective", "composites",
           "mean_composite", "ambiguous", "decided_by", "gap", "n_alternatives", "source"]


def load_pool(run_dir: str) -> Dict[str, object]:
    path = os.path.join(run_dir, "pool", "pool.json")
    if not os.path.exists(path):
        raise FileNotFoundError("{p} not found: run rck-ecdna with --step-penalty > 0 --keep-pool".format(p=path))
    with open(path) as handle:
        return json.load(handle)


def pool_params(pool: Dict[str, object], step_penalty: float) -> DecompositionParams:
    """The ranking settings of the pool's run, at ``step_penalty``."""
    p = pool["params"]
    return DecompositionParams(species_penalty=float(p["species_penalty"]), step_penalty=float(step_penalty),
                               pool_gap=float(p["pool_gap"]), ambiguity_window=float(p["ambiguity_window"]),
                               rank_residual=str(p.get("rank_residual", "first")))


def _species_rows(run_dir: str) -> Optional[pd.DataFrame]:
    path = os.path.join(run_dir, "pool", "decompositions.tsv")
    if not os.path.exists(path):
        return None
    table = pd.read_csv(path, sep="\t")
    return table[table["ordering_rank"] == 1]


def _evaluation(run_dir: str) -> Optional[pd.DataFrame]:
    path = os.path.join(run_dir, "pool", "evaluation.tsv")
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        return pd.read_csv(path, sep="\t")
    except pd.errors.EmptyDataError:
        return None


def rerank(run_dir: str, lambdas: Sequence[float]) -> pd.DataFrame:
    """One row per λ: the rank-1 pool solution, its refitted significant species, its composites
    per truth species and the ambiguity fields, as a run at that λ reports them."""
    pool = load_pool(run_dir)
    entries = pool["entries"]
    species = _species_rows(run_dir)
    frame = _evaluation(run_dir)
    rows = []
    for lam in lambdas:
        params = pool_params(pool, lam)
        order = rank_pool(entries, params)
        row = {"step_penalty": float(lam), "source": "pool" if lam > 0 else "pool, λ=0 emulated"}
        if not order:
            rows.append(row)
            continue
        top = entries[order[0]]
        row.update({"pool_index": top["index"], "steps": top["steps"], "residual": top["residual"],
                    "objective": round(top["objective"], 6)})
        if species is not None:
            mine = species[(species["decomposition_rank"] == top["index"] + 1) & (species["status"] != "low_weight")]
            row["n_species"] = len(mine)
            row["weights"] = ",".join("%g" % w for w in mine["weight"])
        if frame is not None and len(frame):
            scored = frame[frame["decomposition_rank"] == top["index"] + 1]
            row["composites"] = ",".join("%s=%.3f" % (r.truth, r.composite) for r in scored.itertuples())
            row["mean_composite"] = round(float(scored["composite"].mean()), 4) if len(scored) else None
        alts = pool["same_fit"].get(str(top["index"]))
        if params.ambiguity_window > 0 and alts is not None:
            report = ambiguity_summary(top["objective"], top["steps"], top["n_significant"], alts, params, warn=False)
            row.update({"ambiguous": report["ambiguous"], "decided_by": report.get("decided_by"), "gap": report.get("gap"),
                        "n_alternatives": report["n_alternatives"]})
        rows.append(row)
    return pd.DataFrame(rows, columns=COLUMNS)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir", help="rck-ecdna --keep-pool output folder (with pool/)")
    parser.add_argument("--lambdas", type=float, nargs="+", default=None, help="Step penalties (default: the pool's --pool-lambdas)")
    parser.add_argument("-o", "--out", default=None, help="Write the table here (TSV) instead of printing it")
    args = parser.parse_args(argv)
    lambdas = args.lambdas if args.lambdas is not None else load_pool(args.run_dir)["lambdas"]
    table = rerank(args.run_dir, lambdas)
    if args.out:
        table.to_csv(args.out, sep="\t", index=False)
    else:
        table.to_csv(sys.stdout, sep="\t", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
