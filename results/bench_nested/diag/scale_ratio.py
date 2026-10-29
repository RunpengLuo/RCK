"""Observed / truth ecDNA copy number, length-weighted over segments (depth scale) and over junction
traversals (junction-read scale), per rck-ecdna run folder.

    python scale_ratio.py <bench_sweep.tsv> <run folder>...

``scales(run)`` is also used by ../collect.py for the ``scale`` failure stage."""
import sys
import pandas as pd


def scales(run):
    """(segment scale, junction scale, junctions compared) of the run folder ``run``: observed over
    truth ecDNA copy number, length-weighted over ecDNA segments >= 20 kb, and summed over the novel
    junctions that carry ecDNA flow and match a truth junction within 2 kb."""
    t = pd.read_csv(f"{run}/truth_cycles.tsv", sep="\t")
    f = pd.read_csv(f"{run}/flow_segments.tsv", sep="\t")
    f = f[((f.ecdna_A + f.ecdna_B) > 0) & f.obs_total.notna() & (f.length >= 20000)]
    def tcn(r):
        mid = (r.start + r.end) / 2
        return t[(t.chrom == r.chrom) & (t.start <= mid) & (t.end >= mid)].weight.sum()
    f = f.assign(truth=[tcn(r) for r in f.itertuples()])
    f = f[f.truth > 0]
    seg = (f.obs_total * f.length).sum() / (f.truth * f.length).sum() if len(f) else float("nan")
    # junctions: novel adjacencies carrying ecDNA flow; truth copies = flow-independent truth from the walk
    a = pd.read_csv(f"{run}/flow_adjacencies.tsv", sep="\t")
    a = a[(a.kind == "novel")].groupby("aid").agg(obs=("obs_cn", "first"), ec=("ecdna", "sum"), u=("u", "first"), v=("v", "first"))
    a = a[a.ec > 0]
    # truth junction copies: count consecutive-step pairs in truth walks whose ends match u / v within 2 kb
    def ends(sp):
        rows = list(sp.itertuples())
        out = []
        for i, r in enumerate(rows):
            n = rows[(i + 1) % len(rows)]
            e1 = (r.chrom, r.end if r.strand == "+" else r.start)
            e2 = (n.chrom, n.start if n.strand == "+" else n.end)
            out.append((e1, e2, r.weight))
        return out
    pairs = [p for _, sp in t.sort_values(["species", "order"]).groupby("species") for p in ends(sp)]
    def parse(x):
        c, rest = x.split(":"); return c, int(rest.rstrip("+-_").rstrip("+-_").rstrip("+-"))
    tj = []
    for r in a.itertuples():
        (c1, p1), (c2, p2) = parse(r.u), parse(r.v)
        w = sum(wt for e1, e2, wt in pairs if {(e1[0], round(e1[1] / 2000)), (e2[0], round(e2[1] / 2000))} ==
                {(c1, round(p1 / 2000)), (c2, round(p2 / 2000))})
        tj.append(w)
    a = a.assign(truth=tj)
    a = a[a.truth > 0]
    junc = a.obs.sum() / a.truth.sum() if len(a) else float("nan")
    return seg, junc, len(a)


if __name__ == "__main__":
    print("case\tseg_obs/truth\tjunc_obs/truth\tratio(seg/junc)\tmean_composite_l0.1")
    sweep = pd.read_csv(sys.argv[1], sep="\t")
    for run in sys.argv[2:]:
        seg, junc, n = scales(run)
        name = run.rstrip("/").split("/")[-2]
        mc = sweep[(sweep.case == name) & (sweep.step_penalty == 0.1)].mean_composite
        print(f"{name}\t{seg:.2f}\t{junc:.2f} ({n} junctions)\t{seg / junc:.2f}\t{mc.iloc[0] if len(mc) else ''}")
