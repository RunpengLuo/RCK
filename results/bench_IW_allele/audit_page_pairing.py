#!/usr/bin/env python3
"""Which rck_mc_rank1_ec<k>.png pages pair in the absolute frame differently from the within-block pairing (2026-11-21).

For every RCK MC page in topology_plots/compare/<k>_<case>/ and topology_plots/<k>_<case>/, compares the pair the
rck-ecdna-plot page drew from the run's evaluation.tsv (absolute frame: run's own pairing and composite) with the
within-block pair of species_matrix/<set>/<case>/rck_mc_recovery_matrix_block.tsv (the benchmark table's pairing and
composite). Uses the *_absolute.png pages where the folder has been regenerated within blocks (all 16 folders since
2026-11-21). Writes topology_plots/rck_mc_pages_frame_check.tsv. Reads only; no Gurobi.

    python audit_page_pairing.py
"""
import glob
import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allele_bench as alb  # noqa: E402
import make_block_figures as mbf  # noqa: E402
import make_compare_plots as mcp  # noqa: E402
import make_topology_plots as mtp  # noqa: E402

TOPO = os.path.join(HERE, "topology_plots")


def folders():
    for k, (s, n, _w) in enumerate(mcp.CASES):
        yield "compare", s, n, mcp.folder_of(k, n)
    for k, (n, _d) in enumerate(mtp.CASES, start=1):
        yield "numbered", "iw_intpsi", n, os.path.join(TOPO, "{k}_{n}".format(k=k, n=n))


def main():
    import logging
    logging.disable(logging.INFO)
    rows = []
    for where, set_name, name, folder in folders():
        alb.use_set(set_name)
        case = next(c for c in alb.bi.cases() if c["name"] == name)
        ev = pd.read_csv(os.path.join(alb.bi.run_dir(case, "rck_mc"), "evaluation.tsv"), sep="\t")
        ev = ev[ev["decomposition_rank"] == ev["decomposition_rank"].min()]
        old = {r.inferred: (r.truth, float(r.composite)) for r in ev.itertuples() if isinstance(r.inferred, str)}
        m = pd.read_csv(os.path.join(mbf.OUT, set_name, name, "rck_mc_recovery_matrix_block.tsv"), sep="\t")
        new = {r.inferred: (r.truth, float(r.composite)) for r in m.itertuples() if isinstance(r.inferred, str)}
        pages = sorted(p for p in glob.glob(os.path.join(folder, "rck_mc_rank*_*.png")) if not p.endswith("_absolute.png"))
        absolute = sorted(glob.glob(os.path.join(folder, "rck_mc_rank*_absolute.png")))
        for p in (absolute or pages):
            sp = re.match(r"rck_mc_rank\d+_(.+?)(_absolute)?\.png$", os.path.basename(p)).group(1)
            o, n = old.get(sp, (None, None)), new.get(sp, (None, None))
            rows.append({"where": where, "set": set_name, "case": name, "folder": os.path.relpath(folder, TOPO),
                         "page": os.path.basename(p), "inferred": sp, "absolute_truth": o[0] or "",
                         "absolute_composite": o[1], "block_truth": n[0] or "", "block_composite": n[1],
                         "pair_differs": (o[0] or "") != (n[0] or ""),
                         "regenerated_within_blocks": bool(absolute)})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(TOPO, "rck_mc_pages_frame_check.tsv"), sep="\t", index=False, float_format="%.3f")
    print(out.to_string(index=False, float_format=lambda v: "%.3f" % v))


if __name__ == "__main__":
    main()
