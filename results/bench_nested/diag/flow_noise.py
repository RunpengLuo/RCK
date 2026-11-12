"""Step 0 of flow-stage noise handling: what the depth noise looks like, what the flow did with it, and
whether option 1 (a per-window bias variable) or option 2 (relative-error segment weights) would make
the truth's copy numbers cheaper than the flow RCK chose. Offline: reads saved run folders only;
gurobipy is blocked at import.

    python diag/flow_noise.py --bench RCK_RUNS --sigma0 SIGMA0_DIR [-o PREFIX]

RCK_RUNS holds <case>/pool (the bench runs at sigma 0.1); SIGMA0_DIR holds <case>_sigma0/pool. Every
bench case (not demo) is read at sigma 0.1; the cases with a sigma 0 run are read at sigma 0 too.

Per run (``PREFIX_runs.tsv``) and per segment (``PREFIX_segments.tsv``):

* Observation noise. Signed relative error of each ecDNA segment's depth against the truth,
  e_i = (obs_total - flow background) / truth_ec - 1 (segments >= 50 kb carrying truth). Its spread,
  and how strongly neighbouring segments share it: the length-weighted share of the error explained
  by a common shift per 200 kb window (``window_share``, 1 = the error is one shift per window, 0 = it
  is independent per segment), and the lag-1 correlation of e between neighbouring segments (pooled
  over runs in the printout). Junction reads: the same relative error for novel junctions the truth
  uses; if the depth bias did not reach junction reads, their error does not change with sigma.
* What the flow did. Length-weighted copies by which the flow's ecDNA copy number misses the truth's
  integer copy number (truth species weights rounded, x traversals), and by which it misses the
  observation; ``follows_obs`` is the length share of missed segments where the flow sits on the
  observed side of the truth.
* Pricing. Data cost of the flow RCK chose and of the truth's integer flow (same background), under
  - ``current``: the flow MILP's segment term w x |x - obs| / sqrt(obs) with w = min(length/50 kb, 1),
  - ``opt2``: the same with 1 / sqrt(obs + (s obs)^2), s = 0.1 (relative-error weights),
  - ``opt1 k``: current weights against obs x (1 - beta_k) with one beta per 200 kb window,
    |beta| <= 0.3, at cost k x |beta| (k = 5, 10, 20; 10 makes a 10 % window shift cost 1),
  plus the novel-junction read term (0.5 x |flow - obs| / sqrt(obs)) and 0.05 per novel junction used,
  which no option changes. ``d_<scheme>`` = cost(truth) - cost(chosen flow); <= 0 means that scheme
  would not prefer the chosen flow over the truth on these terms.

Approximations, stated: totals only (the MILP fits allele-specific and phase-block pieces where it has
them, which only adds cost to a haplotype-wrong flow); reference-adjacency read terms, the background
prior and linear-end and imputed penalties are left out; windows are independent (no smoothness term).
So a negative d says the option removes the data-cost advantage of the noisy flow, not that a solve
would return the truth.
"""
import sys

sys.modules["gurobipy"] = None  # no licence use possible

import argparse  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from obs_vs_truth import parse_extremity, truth_junctions  # noqa: E402

SEG_CAP, MIN_W = 50_000, 0.02  # ilp_circular.CircularParams.seg_length_cap, min_seg_weight
ALPHA_NOVEL, NOVEL_PENALTY = 0.5, 0.05  # CircularParams.alpha_novel, novel_penalty
WINDOW, BETA_MAX, S = 200_000, 0.3, 0.1
KAPPAS = [1.0, 2.0, 3.0, 5.0, 10.0]
JUNC_MULTS = [2.0, 4.0]  # junction read term x this: trust the reads the depth bias does not reach


def window_beta_cost(x, o, c, kappa):
    """min over beta in [-BETA_MAX, BETA_MAX] of sum c |x - o (1 - beta)| + kappa |beta| (convex,
    piecewise linear: the minimum is at a breakpoint or a bound)."""
    cands = {0.0, BETA_MAX, -BETA_MAX}
    for xi, oi in zip(x, o):
        if oi > 0:
            b = 1.0 - xi / oi
            if abs(b) <= BETA_MAX:
                cands.add(b)
    return min(float(np.sum(c * np.abs(x - o * (1.0 - b)))) + kappa * abs(b) for b in cands)


def run_tables(pool):
    truth = pd.read_csv(os.path.join(pool, "truth_cycles.tsv"), sep="\t")
    seg = pd.read_csv(os.path.join(pool, "flow_segments.tsv"), sep="\t")
    adj = pd.read_csv(os.path.join(pool, "flow_adjacencies.tsv"), sep="\t")
    tint = truth.assign(weight=truth.weight.round())

    def cover(r, t):
        mid = (r.start + r.end) / 2
        return t[(t.chrom == r.chrom) & (t.start <= mid) & (t.end >= mid)].weight.sum()

    seg["truth_ec"] = [cover(r, truth) for r in seg.itertuples()]
    seg["truth_int"] = [cover(r, tint) for r in seg.itertuples()]
    seg["flow_ec"] = seg.ecdna_A + seg.ecdna_B
    seg["bg"] = seg.bg_A + seg.bg_B
    seg["w"] = np.maximum(np.minimum(seg.length / SEG_CAP, 1.0), MIN_W)
    seg["window"] = seg.chrom + ":" + (((seg.start + seg.end) // 2) // WINDOW).astype(str)

    nov = adj[adj.kind == "novel"].groupby("aid").agg(u=("u", "first"), v=("v", "first"), obs=("obs_cn", "first"),
                                                       ec=("ecdna", "sum"), bg=("bg", "sum")).reset_index()
    tj, tji = truth_junctions(truth), truth_junctions(tint)
    key = [frozenset([parse_extremity(r.u), parse_extremity(r.v)]) for r in nov.itertuples()]
    nov["truth"] = [tj.get(k, 0.0) for k in key]
    nov["truth_int"] = [tji.get(k, 0.0) for k in key]
    return seg, nov


def map_flow(seg, nov, seg0, nov0):
    """The sigma 0 run's flow carried onto the sigma 0.1 graph: segment copies (ecDNA, background)
    at each sigma 0.1 segment's midpoint, novel-junction copies by junction (both extremities)."""
    seg = seg.copy()
    nov = nov.copy()

    def at(r, col):
        mid = (r.start + r.end) / 2
        hit = seg0[(seg0.chrom == r.chrom) & (seg0.start <= mid) & (seg0.end >= mid)]
        return float(hit[col].iloc[0]) if len(hit) else 0.0

    seg["truth_int"] = [at(r, "flow_ec") for r in seg.itertuples()]
    seg["bg0"] = [at(r, "bg") for r in seg.itertuples()]
    k0 = {frozenset([parse_extremity(r.u), parse_extremity(r.v)]): (r.ec, r.bg) for r in nov0.itertuples()}
    key = [frozenset([parse_extremity(r.u), parse_extremity(r.v)]) for r in nov.itertuples()]
    nov["truth_int"] = [k0.get(k, (0.0, 0.0))[0] for k in key]
    nov["bg0"] = [k0.get(k, (0.0, 0.0))[1] for k in key]
    matched = sum(1 for k in key if k in k0)
    return seg, nov, matched, len(k0)


def costs(seg, nov, ref_bg=False):
    """{scheme: cost(chosen flow), cost(reference flow)} on the terms named in the docstring. The
    reference is ``truth_int`` (+ ``bg0`` when ref_bg, else the chosen flow's background)."""
    s = seg[seg.obs_total.notna()]
    o = s.obs_total.to_numpy(float)
    xf = (s.flow_ec + s.bg).to_numpy(float)
    xt = (s.truth_int + (s.bg0 if ref_bg else s.bg)).to_numpy(float)
    w = s.w.to_numpy(float)
    cur = w / np.sqrt(np.maximum(o, 1.0))
    rel = w / np.sqrt(np.maximum(o, 1.0) + (S * o) ** 2)

    j = nov[nov.obs.notna()]
    jo = j.obs.to_numpy(float)
    jw = ALPHA_NOVEL / np.sqrt(np.maximum(jo, 1.0))
    jb, nb = (j.bg0, nov.bg0) if ref_bg else (j.bg, nov.bg)
    jdev_f = float(np.sum(jw * np.abs((j.ec + j.bg).to_numpy(float) - jo)))
    jdev_t = float(np.sum(jw * np.abs((j.truth_int + jb).to_numpy(float) - jo)))
    pen_f = NOVEL_PENALTY * int(((nov.ec + nov.bg) > 0).sum())
    pen_t = NOVEL_PENALTY * int(((nov.truth_int + nb) > 0).sum())
    junc_f, junc_t = jdev_f + pen_f, jdev_t + pen_t
    seg_f, seg_t = float(np.sum(cur * np.abs(xf - o))), float(np.sum(cur * np.abs(xt - o)))

    out = {"seg": (seg_f, seg_t), "junc": (junc_f, junc_t), "current": (seg_f + junc_f, seg_t + junc_t),
           "opt2": (float(np.sum(rel * np.abs(xf - o))) + junc_f, float(np.sum(rel * np.abs(xt - o))) + junc_t)}
    win = s.window.to_numpy()
    for k in KAPPAS:
        cf = ct = 0.0
        for wid in np.unique(win):
            m = win == wid
            cf += window_beta_cost(xf[m], o[m], cur[m], k)
            ct += window_beta_cost(xt[m], o[m], cur[m], k)
        out["opt1 %g" % k] = (cf + junc_f, ct + junc_t)
    for mult in JUNC_MULTS:
        out["junc x%g" % mult] = (seg_f + mult * jdev_f + pen_f, seg_t + mult * jdev_t + pen_t)
    return out


def noise_stats(seg, nov):
    s = seg[(seg.length >= 50_000) & (seg.truth_ec > 0) & seg.obs_total.notna()].copy()
    s["e"] = (s.obs_total - s.bg) / s.truth_ec - 1.0
    L = s.length.sum()
    mean = (s.e * s.length).sum() / L
    tot = ((s.e - mean) ** 2 * s.length).sum()
    wmean = s.groupby("window").apply(lambda g: (g.e * g.length).sum() / g.length.sum())
    within = ((s.e - s.window.map(wmean)) ** 2 * s.length).sum()
    s = s.sort_values(["chrom", "start"])
    pairs = [(a.e, b.e) for a, b in zip(s.itertuples(), s.iloc[1:].itertuples())
             if a.chrom == b.chrom and b.start - a.end <= WINDOW]
    used = nov[(nov.truth > 0) & nov.obs.notna()]
    je = ((used.obs - used.bg) / used.truth - 1.0)

    f = seg[(seg.truth_ec > 0) | (seg.flow_ec > 0)]
    miss = f[(f.flow_ec - f.truth_int).abs() >= 1]
    toward = miss[np.sign(miss.flow_ec - miss.truth_int) == np.sign(miss.obs_total - miss.bg - miss.truth_int)]
    Lf = f.length.sum()
    return {"segments": len(s), "seg_err_sd": round(math.sqrt(tot / L), 4) if L else None,
            "seg_err_absmean": round((s.e.abs() * s.length).sum() / L, 4) if L else None,
            "window_share": round(1 - within / tot, 3) if tot > 0 else None,
            "n_pairs": len(pairs), "pairs": pairs,
            "junc_err_absmean": round(float(je.abs().mean()), 4) if len(je) else None,
            "flow_miss_copies": round(float(((f.flow_ec - f.truth_int).abs() * f.length).sum() / Lf), 3) if Lf else None,
            "flow_miss_mb": round(miss.length.sum() / 1e6, 3),
            "follows_obs": round(toward.length.sum() / miss.length.sum(), 3) if miss.length.sum() else None,
            "junc_miss": int(((nov.ec - nov.truth_int).abs() >= 1).sum())}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bench", required=True)
    ap.add_argument("--sigma0", required=True)
    ap.add_argument("-o", "--out", help="prefix for PREFIX_runs.tsv and PREFIX_segments.tsv")
    a = ap.parse_args()
    fails = set(pd.read_csv(os.path.join(HERE, "..", "pipeline.tsv"), sep="\t").query("failure_stage == 'RCK'").case)
    runs = [(c, "0.1", os.path.join(a.bench, c, "pool")) for c in sorted(os.listdir(a.bench)) if not c.startswith("demo")]
    runs += [(c[:-len("_sigma0")], "0", os.path.join(a.sigma0, c, "pool")) for c in sorted(os.listdir(a.sigma0))
             if c.endswith("_sigma0")]
    rows, segs, pairs, tables = [], [], {"0.1": [], "0": []}, {}
    for case, sigma, pool in runs:
        if not os.path.exists(os.path.join(pool, "flow_segments.tsv")):
            continue
        seg, nov = run_tables(pool)
        tables[case, sigma] = (seg, nov)
        st = noise_stats(seg, nov)
        pairs[sigma].extend(st.pop("pairs"))
        row = {"case": case, "sigma": sigma, "fails": case in fails, "n_ec_segments": int((seg.flow_ec > 0).sum())}
        row.update(st)
        for scheme, (cf, ct) in costs(seg, nov).items():
            row["d_" + scheme] = round(ct - cf, 3)
        rows.append(row)
        segs.append(seg.assign(case=case, sigma=sigma))
    # the sigma 0 flow (which decomposed right in 5 of 8) priced on the sigma 0.1 observations
    ref = []
    for (case, sigma), (seg0, nov0) in tables.items():
        if sigma != "0" or (case, "0.1") not in tables:
            continue
        seg, nov = tables[case, "0.1"]
        mseg, mnov, matched, n0 = map_flow(seg, nov, seg0, nov0)
        r = {"case": case, "junctions_matched": "%d/%d" % (matched, n0),
             "flow_diff_copies": round(float(((mseg.flow_ec - mseg.truth_int).abs() * mseg.length).sum() / mseg.length.sum()), 3)}
        for scheme, (cf, ct) in costs(mseg, mnov, ref_bg=True).items():
            r["d_" + scheme] = round(ct - cf, 3)
        ref.append(r)
    ref = pd.DataFrame(ref)
    t = pd.DataFrame(rows)
    if a.out:
        ref.to_csv(a.out + "_sigma0flow.tsv", sep="\t", index=False)
        t.to_csv(a.out + "_runs.tsv", sep="\t", index=False)
        pd.concat(segs).to_csv(a.out + "_segments.tsv", sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 200, "display.max_columns", 40):
        cols = ["case", "sigma", "fails", "n_ec_segments", "seg_err_absmean", "window_share", "junc_err_absmean",
                "flow_miss_copies", "flow_miss_mb", "follows_obs", "junc_miss"] + [c for c in t.columns if c.startswith("d_")]
        print(t[cols].to_string(index=False))
    print()
    for sigma, p in pairs.items():
        if len(p) > 2:
            x, y = np.array(p).T
            print("sigma %s: lag-1 correlation of segment depth error between neighbours (<= 200 kb apart): %.3f over %d pairs"
                  % (sigma, np.corrcoef(x, y)[0, 1], len(p)))
    for sigma, g in t.groupby("sigma"):
        print("sigma %s, %d runs: mean |depth error| %.4f, window share %.3f (median), mean |junction error| %.4f"
              % (sigma, len(g), g.seg_err_absmean.mean(), g.window_share.median(), g.junc_err_absmean.mean()))
    b = t[t.sigma == "0.1"]
    for grp, g in (("failing", b[b.fails]), ("passing", b[~b.fails])):
        print("\nsigma 0.1, %s bench cases (%d): flow misses truth by %.3f copies (length-weighted mean), follows obs %.2f"
              % (grp, len(g), g.flow_miss_copies.mean(), g.follows_obs.mean()))
        for c in [c for c in t.columns if c.startswith("d_")]:
            print("   %-10s truth cheaper or equal in %2d of %2d; median d %.3f" % (c, int((g[c] <= 1e-9).sum()), len(g), g[c].median()))
    print("\nsigma 0 flow priced on the sigma 0.1 observations (d = cost(sigma 0 flow) - cost(sigma 0.1 flow)):")
    print(ref.to_string(index=False))
    print("gurobipy loaded: %s" % (sys.modules.get("gurobipy") is not None))


if __name__ == "__main__":
    main()
