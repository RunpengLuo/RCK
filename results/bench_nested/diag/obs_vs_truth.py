"""Which observations are closer to the truth: segment depth or junction reads? One row per run.

    python obs_vs_truth.py RUN_DIR [RUN_DIR ...] > obs_vs_truth.tsv

RUN_DIR is an rck-ecdna output folder with --truth-dir (truth_cycles.tsv, flow_segments.tsv,
flow_adjacencies.tsv). Truth copy numbers: per flow segment, the sum over truth species of weight x
traversals at the segment midpoint (both haplotypes, as flow_vs_truth.py); per novel junction, the sum
of weight x traversals of the junction in the truth walks (consecutive truth segments that are not
reference neighbours on the same strand; the walk is closed). Columns: median and length-weighted mean
relative error |obs / truth - 1| of segments (obs_total, segments ≥ 50 kb with truth > 0) and of novel
junctions (obs_cn); the same for the flow's ecDNA + background copies; and how many junctions the truth
walks use that the graph has not got (unmatched).
"""
import os
import sys

import pandas as pd


def extremities(row, entering):
    """(chrom, pos, side) of a truth segment's extremity: the one a walk leaves by (entering False) or
    enters by (True), in the flow tables' notation (start '-', end '+')."""
    plus = row.strand == "+"
    if entering:
        return (row.chrom, row.start, "-") if plus else (row.chrom, row.end, "+")
    return (row.chrom, row.end, "+") if plus else (row.chrom, row.start, "-")


def parse_extremity(text):
    chrom, rest = text.rsplit(":", 1)
    pos, side = rest.split("_")[0][:-1], rest.split("_")[0][-1]
    return chrom, int(pos), side


def truth_junctions(truth):
    cn = {}
    for _species, walk in truth.groupby("species", sort=False):
        walk = walk.sort_values("order")
        rows = list(walk.itertuples())
        w = rows[0].weight
        for a, b in zip(rows, rows[1:] + rows[:1]):
            if a.chrom == b.chrom and a.strand == b.strand and (
                    (a.strand == "+" and b.start == a.end + 1) or (a.strand == "-" and a.start == b.end + 1)):
                continue
            key = frozenset([extremities(a, False), extremities(b, True)])
            cn[key] = cn.get(key, 0.0) + w
    return cn


def row(run):
    truth = pd.read_csv(os.path.join(run, "truth_cycles.tsv"), sep="\t")
    seg = pd.read_csv(os.path.join(run, "flow_segments.tsv"), sep="\t")
    adj = pd.read_csv(os.path.join(run, "flow_adjacencies.tsv"), sep="\t")

    def truth_cn(r):
        mid = (r.start + r.end) / 2
        return truth[(truth.chrom == r.chrom) & (truth.start <= mid) & (truth.end >= mid)].weight.sum()

    seg["truth"] = [truth_cn(r) for r in seg.itertuples()]
    seg = seg[(seg.truth > 0) & (seg.length >= 50_000) & seg.obs_total.notna()].copy()
    seg["flow"] = seg.ecdna_A + seg.ecdna_B
    seg["err_obs"] = ((seg.obs_total - seg.bg_A - seg.bg_B) / seg.truth - 1).abs()
    seg["err_flow"] = (seg.flow / seg.truth - 1).abs()

    nov = adj[adj.kind == "novel"].groupby("aid").agg(u=("u", "first"), v=("v", "first"), obs=("obs_cn", "first"),
                                                       flow=("ecdna", "sum"), bg=("bg", "sum")).reset_index()
    tj = truth_junctions(truth)
    nov["truth"] = [tj.get(frozenset([parse_extremity(r.u), parse_extremity(r.v)]), 0.0) for r in nov.itertuples()]
    used = nov[(nov.truth > 0) & nov.obs.notna()].copy()
    used["err_obs"] = ((used.obs - used.bg) / used.truth - 1).abs()
    used["err_flow"] = (used.flow / used.truth - 1).abs()
    matched = {frozenset([parse_extremity(r.u), parse_extremity(r.v)]) for r in nov.itertuples()}
    L = seg.length.sum()
    return {"run": os.path.basename(os.path.dirname(os.path.normpath(run))),
            "segments": len(seg), "seg_obs_median": round(seg.err_obs.median(), 3),
            "seg_obs_lenmean": round((seg.err_obs * seg.length).sum() / L, 3),
            "seg_flow_lenmean": round((seg.err_flow * seg.length).sum() / L, 3),
            "junctions": len(used), "junc_obs_median": round(used.err_obs.median(), 3),
            "junc_obs_mean": round(used.err_obs.mean(), 3), "junc_flow_mean": round(used.err_flow.mean(), 3),
            "unmatched": sum(1 for k in tj if k not in matched)}


def main():
    pd.DataFrame([row(r) for r in sys.argv[1:]]).to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()
