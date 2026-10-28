"""One line per pool solution: residual, significant species, weights, steps, composites."""
import sys
import pandas as pd
for run in sys.argv[1:]:
    d = pd.read_csv(f"{run}/pool/pool/decompositions.tsv", sep="\t")
    e = pd.read_csv(f"{run}/pool/pool/evaluation.tsv", sep="\t")
    print("==", run.rstrip("/").split("/")[-1])
    rows = []
    for k, g in d.groupby("decomposition_rank"):
        g = g.drop_duplicates("species")
        sig = g[g["status"] == "species"]
        ev = e[e["decomposition_rank"] == k]
        comp = ",".join(f"{t}={c:.3f}" for t, c in zip(ev["truth"], ev["composite"]))
        rows.append((k, round(g["residual"].iloc[0], 3), len(sig), ",".join(f"{w:g}" for w in sig["weight"]),
                     int(sig["n_steps"].sum()), len(g) - len(sig), comp, round(ev["composite"].mean(), 4)))
    t = pd.DataFrame(rows, columns=["rank", "residual", "n_sig", "weights", "steps", "n_low", "composites", "mean"])
    print(t.sort_values(["residual", "n_sig", "steps"]).head(12).to_string(index=False))
    print("BEST by mean:", t.sort_values("mean", ascending=False).head(3)[["rank","n_sig","weights","steps","composites","mean"]].to_string(index=False))
    print("solutions:", len(t), "| by n_sig:", t.groupby("n_sig").size().to_dict(), "| best mean among 2-species:",
          t[t.n_sig == 2]["mean"].max() if (t.n_sig == 2).any() else None)
