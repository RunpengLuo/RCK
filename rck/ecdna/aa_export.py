"""Export an ``rck-ecdna`` answer to AA format for the method-agnostic scorer (``aa_eval``).

    python -m rck.ecdna.aa_export RUN_DIR OUT [--pool-index I | --lambda L] [--alternatives 2]

``rck-ecdna`` already writes ``aa/<name>_graph.txt`` and ``aa/<name>_decomp<k>_cycles.txt``,
but those files lack what the scorer needs to reproduce ``rck-ecdna-eval`` exactly: the
``low_weight`` status, a species' other orderings, and the solutions of a kept pool
(``--keep-pool``, which the λ > 0 results come from). This module writes, from
``decompositions.tsv`` (or ``pool/decompositions.tsv``):

* ``<name>_cycles.txt``: the answer, one cycle per species with ``Copy_count`` = the refitted
  weight, ``Haplotypes=``, ``Status=``, ``Species=`` and the other orderings
  (``Alt_segments=`` / ``Alt_haplotypes=``); ``Interval`` lines are the run's regions;
* ``<name>_alt<k>_cycles.txt``: the first k same-fit alternatives of the answer (the set
  metric of results/bench_nested/diag/grow_rules.py), for best-of scoring;
* ``<name>_graph.txt``: the run's AA graph;
* ``<name>_segment_cn.tsv``: allele copy numbers per segment (ecDNA + background, from
  ``flow_segments.tsv``).

Coordinates: the scorer reads AA segments as ``start - 1`` and ``end - 1`` (0-based inclusive,
how the simulator's 1-based AA truth is read). ``rck-ecdna``'s walks and ``aa/`` files carry the
subgraph's segment coordinates unchanged, so every position is written here as the walk's
coordinate + 1: the scorer then reads back exactly the steps ``rck-ecdna-eval`` reads from
``decompositions.tsv``, and its scores are identical.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from typing import Dict, List, Optional, Sequence

import pandas as pd

SHIFT = 1  # walk coordinate -> AA coordinate (see the module docstring)


def _token(text: str):
    chrom, span, hap, strand = text.rsplit(":", 3)
    start, end = span.split("-")
    return chrom, int(start), int(end), hap, strand


def write_cycles(path: str, rows: pd.DataFrame, regions: Sequence[str]) -> None:
    """One decomposition's species (``rows``: its decompositions.tsv rows, every ordering)."""
    walks = {}
    for (species, _), group in rows.groupby(["species", "ordering_rank"], sort=False):
        walks.setdefault(species, []).append([_token(t) for t in str(group["walk"].iloc[0]).split(",")])
    segs = sorted({(c, s, e) for ws in walks.values() for w in ws for c, s, e, _, __ in w})
    number = {seg: i for i, seg in enumerate(segs, start=1)}
    lines = ["# rck-ecdna answer for the method-agnostic scorer (rck.ecdna.aa_export); coordinates = walk + {s}".format(s=SHIFT)]
    for i, text in enumerate(regions, start=1):
        chrom, span = text.rsplit(":", 1)
        a, b = span.split("-")
        lines.append("Interval\t{i}\t{c}\t{a}\t{b}".format(i=i, c=chrom, a=int(a) + SHIFT, b=int(b) + SHIFT))
    lines.append("List of cycle segments")
    for seg in segs:
        lines.append("Segment\t{n}\t{c}\t{s}\t{e}".format(n=number[seg], c=seg[0], s=seg[1] + SHIFT, e=seg[2] + SHIFT))
    first = rows.drop_duplicates("species")
    for k, r in enumerate(first.to_dict("records"), start=1):
        linear = str(r.get("topology", "circular")) == "linear"

        def encode(walk):
            segtext = ",".join("{n}{st}".format(n=number[(c, s, e)], st=st) for c, s, e, _, st in walk)
            haps = ",".join(h for _, __, ___, h, ____ in walk)
            return ("0+," + segtext + ",0-", ".," + haps + ",.") if linear else (segtext, haps)

        main_segs, main_haps = encode(walks[r["species"]][0])
        fields = ["Cycle={k}".format(k=k), "Copy_count={w!r}".format(w=float(r["weight"])), "Segments=" + main_segs,
                  "Haplotypes=" + main_haps, "Status=" + str(r.get("status", "species")), "Species=" + str(r["species"])]
        alts = [encode(w) for w in walks[r["species"]][1:]]
        if alts:
            fields += ["Alt_segments=" + "|".join(a for a, _ in alts), "Alt_haplotypes=" + "|".join(h for _, h in alts)]
        lines.append(";".join(fields))
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def shift_graph(src: str, dst: str) -> None:
    """Copy an ``rck-ecdna`` AA graph with every position + SHIFT."""
    from rck.ecdna.aa_format import parse_breakend

    def be(text):
        c, p, s = parse_breakend(text)
        return "{c}:{p}{s}".format(c=c, p=p + SHIFT, s=s)

    out = []
    with open(src) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] == "sequence":
                f[1], f[2] = be(f[1]), be(f[2])
            elif f[0] in ("concordant", "discordant", "source"):
                a, b = f[1].split("->")
                f[1] = be(a) + "->" + be(b)
            out.append("\t".join(f))
    with open(dst, "w") as handle:
        handle.write("\n".join(out) + "\n")


def write_segment_cn(src: str, dst: str) -> None:
    t = pd.read_csv(src, sep="\t")
    a = t["ecdna_A"] + t["bg_A"]
    b = t["ecdna_B"] + t["bg_B"]
    pd.DataFrame({"chrom": t["chrom"], "start": t["start"] + SHIFT, "end": t["end"] + SHIFT, "cn": a + b, "cn_a": a,
                  "cn_b": b}).to_csv(dst, sep="\t", index=False)


def run_alternatives(run_dir: str, k: int) -> List[int]:
    """Decomposition ranks of rank 1's first k same-fit alternatives in a run (summary.json
    ``ambiguity``); alternatives the run did not report (rank null) are skipped, as unmatched
    ones are in grow_rules.py's set metric."""
    with open(os.path.join(run_dir, "summary.json")) as handle:
        amb = json.load(handle).get("ambiguity") or {}
    return [a["rank"] for a in (amb.get("alternatives") or [])[:k] if a.get("rank") is not None]


def pool_alternatives(pool: Dict[str, object], top: int, k: int) -> List[int]:
    """Pool indices of pool entry ``top``'s first k same-fit alternatives (pool.json
    ``same_fit``), matched to entries by objective, steps and significant weights as
    results/bench_nested/diag/grow_rules.py's set metric does."""
    entries = pool["entries"]
    sig = lambda e: sorted(w for w, x in zip(e["weights_int"], e["low_weight"]) if not x)  # noqa: E731
    out = []
    for alt in (pool["same_fit"].get(str(top)) or [])[:k]:
        hit = [e["index"] for e in entries if abs(e["objective"] - alt["objective"]) < 1e-3
               and e["steps"] == alt["steps"] and sig(e) == sorted(alt["weights"])]
        if hit:
            out.append(hit[0])
    return out


def export(run_dir: str, outdir: str, pool_index: Optional[int] = None, alternatives: int = 2, name: Optional[str] = None) -> Dict[str, object]:
    """Write the answer of ``run_dir`` (rank 1; or pool entry ``pool_index`` of its kept pool)
    and its alternatives to ``outdir``. Returns the paths for :func:`rck.ecdna.aa_eval.score`."""
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(run_dir, "summary.json")) as handle:
        summary = json.load(handle)
    regions = summary["graph"]["regions"]
    name = name or summary["graph"]["name"]
    if pool_index is None:
        table = pd.read_csv(os.path.join(run_dir, "decompositions.tsv"), sep="\t")
        ranks = [int(table["decomposition_rank"].min())] if len(table) else []
        ranks += run_alternatives(run_dir, alternatives)
    else:
        table = pd.read_csv(os.path.join(run_dir, "pool", "decompositions.tsv"), sep="\t")
        with open(os.path.join(run_dir, "pool", "pool.json")) as handle:
            pool = json.load(handle)
        ranks = [pool_index + 1] + [i + 1 for i in pool_alternatives(pool, pool_index, alternatives)]
    paths = {"cycles": os.path.join(outdir, name + "_cycles.txt"), "alt_cycles": [], "ranks": ranks}
    for k, rank in enumerate(ranks):
        path = paths["cycles"] if k == 0 else os.path.join(outdir, "{n}_alt{k}_cycles.txt".format(n=name, k=k))
        write_cycles(path, table[table["decomposition_rank"] == rank], regions)
        if k:
            paths["alt_cycles"].append(path)
    graphs = glob.glob(os.path.join(run_dir, "aa", "*_graph.txt"))
    if graphs:
        paths["graph"] = os.path.join(outdir, name + "_graph.txt")
        shift_graph(graphs[0], paths["graph"])
    if os.path.exists(os.path.join(run_dir, "flow_segments.tsv")):
        paths["segment_cn"] = os.path.join(outdir, name + "_segment_cn.tsv")
        write_segment_cn(os.path.join(run_dir, "flow_segments.tsv"), paths["segment_cn"])
    return paths


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("outdir")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--pool-index", type=int, default=None, help="export this entry of the kept pool (0-based)")
    group.add_argument("--lambda", dest="lam", type=float, default=None, help="export the pool's rank 1 at this --step-penalty")
    parser.add_argument("--alternatives", type=int, default=2, help="same-fit alternatives to export for best-of (default 2)")
    args = parser.parse_args(argv)
    index = args.pool_index
    if args.lam is not None:
        from rck.ecdna.rerank import rerank

        index = int(rerank(args.run_dir, [args.lam])["pool_index"].iloc[0])
    paths = export(args.run_dir, args.outdir, index, args.alternatives)
    json.dump(paths, sys.stdout, indent=2)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
