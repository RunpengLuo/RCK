"""Per flow segment: truth ecDNA copies per allele (sum of species weight x traversals) against the
flow's ecdna_A / ecdna_B. Usage: flow_vs_truth.py RUN_DIR [RUN_DIR ...] (an rck-ecdna lambda_0 folder)."""
import sys
import pandas as pd
for run in sys.argv[1:]:
    t = pd.read_csv(f"{run}/truth_cycles.tsv", sep="\t")
    f = pd.read_csv(f"{run}/flow_segments.tsv", sep="\t")
    f = f[(f.ecdna_A + f.ecdna_B) > 0].copy()
    def truth_cn(r, hap):
        mid = (r.start + r.end) / 2
        hit = t[(t.chrom == r.chrom) & (t.start <= mid) & (t.end >= mid) & (t.hap == hap)]
        return hit.weight.sum()
    f["truth_A"] = [truth_cn(r, "A") for r in f.itertuples()]
    f["truth_B"] = [truth_cn(r, "B") for r in f.itertuples()]
    f["err"] = (f.ecdna_A - f.truth_A).abs() + (f.ecdna_B - f.truth_B).abs()
    L = f.length.sum()
    print("==", run.split("/")[-2], "| species weights:", t.groupby("species").weight.first().round(2).to_dict())
    print("  length-weighted |flow - truth| per segment: %.2f copies; segments off by > 1.5 copies: %d of %d (%.2f Mb)"
          % ((f.err * f.length).sum() / L, (f.err > 1.5).sum(), len(f), f[f.err > 1.5].length.sum() / 1e6))
    print(f[f.err > 1.5][["chrom", "start", "end", "length", "obs_A", "obs_B", "ecdna_A", "ecdna_B", "truth_A", "truth_B"]]
          .round(2).head(8).to_string(index=False))
