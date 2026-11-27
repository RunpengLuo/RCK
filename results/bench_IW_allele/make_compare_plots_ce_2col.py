#!/usr/bin/env python3
"""Two-column, print-size truth / RCK MC / CoRAL / CE topology comparison (Yubi, 2026-11-26).

The same figure as make_compare_plots_ce.py (commit 3414b84), laid out for a paper page: drawn at the final printed
size (6.5 in wide = full text width, at most 8.5 in tall) so font sizes are real, no text under 6.5 pt.

  left column : Truth (top), RCK MC        right column: CoRAL (top), CE

Both columns share one genome x-axis (same regions, same scale) and one row height; the two columns start at the same
height so truth sits beside CoRAL. Each walk has a one- or two-line label above it ("ec1 → ecdna1, 16.9 copies; psi
..."; → = the table scorer's pairing). No figure-wide title (it goes in the LaTeX caption); one shared legend at the
bottom. Everything else is drawn as in make_compare_plots_ce.py: colours (A blue, B orange, grey for CoRAL/CE), arrows
for strand, connectors in traversal order, dotted closing join, RCK MC haplotypes after the scorer's one A/B swap per
phase block, steps over LONG_STEP (background in linear walks) clipped to the axis.

Data loading, pairing and the step merging are those of make_compare_plots_ce / make_compare_plots (imported, not
edited); draw_topology and position_ticks (rck.ecdna.plot) are copied here with print-size fonts. If the walks do
not fit in 8.5 in at the minimum row height (MIN_ROW), the script says so and reports the height needed instead of
shrinking the fonts.

Writes topology_plots/compare_ce/<case>/: topology_compare_ce_2col.pdf (vector), topology_compare_ce_2col.png
(300 dpi), topology_compare_ce_2col_README.txt. Existing files there are not touched. Run folders are only read; no
Gurobi; well under a minute per case.

    python make_compare_plots_ce_2col.py                       # the cases in make_compare_plots_ce.CASES
    python make_compare_plots_ce_2col.py --case SET CASE
"""
import argparse
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_compare_plots as mcp  # noqa: E402
import make_compare_plots_ce as mce  # noqa: E402
import make_topology_plots as mtp  # noqa: E402
import allele_bench as alb  # noqa: E402

OUT = mce.OUT
STEM = "topology_compare_ce_2col"

# Print geometry (inches) and fonts (pt), all at the final printed size.
FIG_W, MAX_H = 6.5, 8.5
LEFT, RIGHT, GUTTER = 0.20, 0.06, 0.30  # LEFT/GUTTER leave room for the step numbers left of each column
TOP = 0.04
TITLE_H, LABEL_LINE_H, WALK_GAP, PANEL_GAP = 0.20, 0.115, 0.07, 0.14
BOTTOM = 0.95  # tick labels, "position (Mb)", chromosome names, legend
BASE = 0.35  # rows of clearance under the last step, so the closing join clears the axis line
MAX_ROW, MIN_ROW = 0.16, 6.5 / 72 + 0.02  # a row must hold its 6.5 pt step number
FS_TITLE, FS_LABEL, FS_STEP, FS_TICK, FS_AXIS, FS_LEGEND = 8.0, 7.0, 6.5, 6.5, 7.0, 7.0
COLUMNS = [["truth", "RCK MC"], ["CoRAL", "CE"]]
TITLES = {"truth": "Truth", "RCK MC": "RCK MC", "CoRAL": "CoRAL", "CE": "CE"}


def steps_of(c):
    """make_compare_plots_ce's walk for drawing (reference-contiguous steps merged; the walk is unchanged)."""
    from rck.ecdna import plot as rp

    out = []
    for x in c.steps:
        if out:
            p = out[-1]
            same = p.chrom == x.chrom and p.hap == x.hap and p.strand == x.strand
            if same and x.strand == "+" and x.start == p.end + 1:
                out[-1] = rp.PStep(p.chrom, p.start, x.end, p.hap, p.strand)
                continue
            if same and x.strand == "-" and x.end == p.start - 1:
                out[-1] = rp.PStep(p.chrom, x.start, p.end, p.hap, p.strand)
                continue
        out.append(rp.PStep(x.chrom, x.start, x.end, x.hap, x.strand))
    return out


def load(case):
    """The walks of make_compare_plots_ce.topology_compare_ce, same loading and pairing: returns
    (rows [(kind, cycle, [label parts])], pairs {kind: {inferred: truth}})."""
    import pandas as pd
    from rck.ecdna import bench_metrics as bm, evaluate as ev

    bi = alb.bi
    truth = ev.load_truth(case["truth"])
    psi_t, bulk_t = mtp.clone_psi(case["truth"])
    run = bi.run_dir(case, "rck_mc")
    weights = bm.multiclone_bulk_weights(run, bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv")))
    rck = mcp.strict_cycles(mtp.answer_files(case, "rck_mc")[0], weights)
    blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
    by_t, by_c = {t.name: t for t in truth}, {c.name: c for c in rck}
    swaps = ev.block_swaps([(by_t[a], by_c[b]) for a, b in ev._match(truth, rck, False).items()], blocks)
    rck = [ev.relabel(c, blocks, swaps) for c in rck]
    pairs = {"RCK MC": mcp.scorer_pairs(case, "rck_mc")}
    check = {v: k for k, v in ev._match(truth, rck, True).items()}
    if check != pairs["RCK MC"]:
        print("{n}: RCK MC pairing: scorer {s} vs re-matched {m}; drawing the scorer's".format(n=case["name"], s=pairs["RCK MC"], m=check))
    dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t")
    dec = dec[dec["decomposition_rank"] == dec["decomposition_rank"].min()].drop_duplicates("species").set_index("species")
    clones = [c[len("psi_"):] for c in dec.columns if c.startswith("psi_")]
    blind = {}
    for label, method, variant in mce.BLIND:
        blind[label] = mcp.strict_cycles(alb.answer_files(case, method, variant)[0])
        pairs[label] = mcp.scorer_pairs(case, method, variant)
        check = {v: k for k, v in ev._match(truth, blind[label], False).items()}
        if check != pairs[label]:
            print("{n}: {l} pairing: scorer {s} vs re-matched {m}; drawing the scorer's".format(
                n=case["name"], l=label, s=pairs[label], m=check))

    rows = []
    for t in truth:
        base = t.name.split(".")[0]
        psi = ", ".join("{c} {v:g}".format(c=c, v=v) for c, v in sorted(psi_t.get(base, {}).items()))
        rows.append(("truth", t, ["{n}, {b:.1f} copies".format(n=t.name, b=bulk_t.get(base, t.weight)), "psi " + psi]))
    for c in rck:
        psi = ", ".join("{k} {v:g}".format(k=k, v=float(dec.loc[c.name, "psi_" + k])) for k in clones) if c.name in dec.index else ""
        rows.append(("RCK MC", c, ["{n} → {p}, {b:.1f} copies".format(n=c.name, p=pairs["RCK MC"].get(c.name, "no pair"), b=c.weight),
                                   "psi " + psi]))
    for label, _, _ in mce.BLIND:
        for c in blind[label]:
            rows.append((label, c, ["{n} → {p}, {b:.1f} copies".format(n=c.name, p=pairs[label].get(c.name, "no pair"), b=c.weight),
                                    "cycle" if c.circular else "linear walk"]))
    return rows, pairs


def text_width_in(text, size):
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties

    return TextPath((0, 0), text, size=size, prop=FontProperties()).get_extents().width / 72.0


def draw_topology(ax, axis, steps, circular):
    """rck.ecdna.plot.draw_topology with the step numbers at FS_STEP and no y label (the label goes above the walk)."""
    import numpy as np
    from matplotlib.patches import Rectangle
    from rck.ecdna import plot as rp

    lows, highs = rp._intervals(axis, steps)
    count = len(steps)
    for i, s in enumerate(steps):
        low, high = lows[i], highs[i]
        if not (np.isfinite(low) and np.isfinite(high)):
            continue
        y = count - i - 0.5
        ax.add_patch(Rectangle((low, y - 0.28), max(high - low, 1.0), 0.56, facecolor=rp.HAP_COLORS.get(s.hap, "#9a9a9a"),
                               edgecolor="white", linewidth=0.4, zorder=3))
        forward = s.strand == "+"
        ax.annotate("", xy=(high if forward else low, y), xytext=(low if forward else high, y),
                    arrowprops={"arrowstyle": "-|>", "color": "white", "lw": 0.6, "shrinkA": 0, "shrinkB": 0,
                                "mutation_scale": 6}, zorder=4)
        ax.text(axis.segments[0].axis_start, y, str(i + 1) + " ", ha="right", va="center", fontsize=FS_STEP,
                color=rp.INK_MUTED, transform=ax.transData, clip_on=False)
    last = count if circular else count - 1
    for i in range(last):
        nxt = (i + 1) % count
        here, there = steps[i], steps[nxt]
        exit_x = highs[i] if here.strand == "+" else lows[i]
        entry_x = lows[nxt] if there.strand == "+" else highs[nxt]
        if not (np.isfinite(exit_x) and np.isfinite(entry_x)):
            continue
        y0, y1 = count - i - 0.5, count - nxt - 0.5
        mid = y0 - 0.5
        ax.plot([exit_x, exit_x, entry_x, entry_x], [y0, mid, mid, y1], color=rp.INK_MUTED, linewidth=0.6,
                linestyle=":" if nxt == 0 else "-", zorder=2)
    ax.set_ylim(-BASE, count + 0.05)
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)


def position_ticks(ax, axis):
    """rck.ecdna.plot.position_ticks with print-size fonts and the chromosome names a fixed distance under the axis."""
    import math
    from matplotlib.transforms import offset_copy
    from rck.ecdna import plot as rp

    positions, labels = [], []
    for seg in axis.segments:
        shift = seg.axis_start - seg.raw_start
        step = rp._nice_step(seg.raw_end - seg.raw_start, target=3)
        first = int(math.ceil(seg.raw_start / step)) * step
        for p in range(first, int(seg.raw_end), step):
            positions.append(p + shift)
            labels.append("{v:g}".format(v=round(p / 1e6, 3)))
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=FS_TICK, rotation=45, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=2, width=0.5, pad=1.5)
    ax.set_xlabel("position (Mb)", fontsize=FS_AXIS, color=rp.INK_MUTED, labelpad=1)
    below = offset_copy(ax.get_xaxis_transform(), fig=ax.figure, y=-36, units="points")
    for x, name in zip(axis.xtick_chrs, axis.xlab_chrs):
        ax.text(x, 0, str(name), ha="center", va="top", fontsize=FS_AXIS, color=rp.INK, transform=below)


def build(case, rows, work, out_stem):
    """Lay out and save the figure; returns (fig size, layout dict, clipped)."""
    from rck.ecdna import plot as rp

    plt, _, cnplot = rp.import_plotting(mtp.CNPLOT)
    plt.rcParams["pdf.fonttype"] = 42  # editable text in the PDF
    every = [st for _, c, _ in rows for st in steps_of(c)]
    short = [st for st in every if st.end - st.start + 1 <= mcp.LONG_STEP]
    covered = {st.chrom for st in short}
    regions = rp.axis_regions(short + [st for st in every if st.chrom not in covered])
    axis = rp.build_axis(cnplot, regions, work)
    span = {}
    for chrom, a, b in regions:
        lo, hi = span.get(chrom, (a, b))
        span[chrom] = (min(lo, a), max(hi, b))
    clipped = {}

    def drawn(kind, c):
        out = []
        for st in steps_of(c):
            lo, hi = span[st.chrom]
            if st.start < lo or st.end > hi - 1:
                clipped.setdefault(kind, []).append(c.name)
                st = rp.PStep(st.chrom, max(st.start, lo), min(st.end, hi - 1), st.hap, st.strand)
            out.append(st)
        return out

    drawn_steps = [drawn(k, c) for k, c, _ in rows]
    col_w = (FIG_W - LEFT - RIGHT - GUTTER) / 2.0
    # Walk labels: one line if it fits the column, else two lines; a clipped walk says so.
    labels = []
    for (kind, c, parts), _ in zip(rows, drawn_steps):
        parts = list(parts)
        if c.name in clipped.get(kind, []):
            parts[-1] += "; steps > {m} Mb clipped".format(m=mcp.LONG_STEP // 1_000_000)
        one = parts[0] + "; " + parts[1]
        labels.append([one] if text_width_in(one, FS_LABEL) <= col_w else parts)
    # Vertical budget per column: fixed (titles, labels, gaps) + rows * row height.
    fixed, nrows = [], []
    for kinds in COLUMNS:
        idx = [i for i, r in enumerate(rows) if r[0] in kinds]
        fixed.append(TOP + len(kinds) * TITLE_H + (len(kinds) - 1) * PANEL_GAP
                     + sum(len(labels[i]) * LABEL_LINE_H + WALK_GAP for i in idx) + BOTTOM)
        nrows.append(sum(len(drawn_steps[i]) + BASE + 0.05 for i in idx))
    row_h = min([MAX_ROW] + [(MAX_H - f) / n for f, n in zip(fixed, nrows)])
    if row_h < MIN_ROW:
        need = max(f + n * MIN_ROW for f, n in zip(fixed, nrows))
        raise SystemExit("{n}: does not fit {h} in at {s} pt; smallest height that works: {x:.2f} in".format(
            n=case["name"], h=MAX_H, s=FS_STEP, x=need))
    height = max(f + n * row_h for f, n in zip(fixed, nrows))

    fig = plt.figure(figsize=(FIG_W, height))
    first, last_ax, layout = None, [], {"row_h_in": row_h, "col_w_in": col_w, "rows_per_column": [sum(len(drawn_steps[i]) for i, r in enumerate(rows) if r[0] in k) for k in COLUMNS]}
    for col, kinds in enumerate(COLUMNS):
        x0 = LEFT + col * (col_w + GUTTER)
        y = height - TOP  # cursor, inches from the bottom
        ax = None
        for p, kind in enumerate(kinds):
            if p:
                y -= PANEL_GAP
            fig.text(x0 - LEFT + 0.02, y, TITLES[kind], ha="left", va="top", fontsize=FS_TITLE, fontweight="bold",
                     color=rp.INK, transform=fig.dpi_scale_trans)
            y -= TITLE_H
            for i, (k, c, _) in enumerate(rows):
                if k != kind:
                    continue
                for line in labels[i]:
                    fig.text(x0, y, line, ha="left", va="top", fontsize=FS_LABEL, color=rp.INK, transform=fig.dpi_scale_trans)
                    y -= LABEL_LINE_H
                h = (len(drawn_steps[i]) + BASE + 0.05) * row_h
                y -= h
                ax = fig.add_axes([x0 / FIG_W, y / height, col_w / FIG_W, h / height], sharex=first)
                first = first or ax
                draw_topology(ax, axis, drawn_steps[i], c.circular)
                ax.set_xlim(axis.ch_coords[0], axis.chr_end)
                for boundary in axis.chr_boundaries:
                    ax.axvline(boundary, color=rp.INK_MUTED, linewidth=0.6, zorder=0)
                ax.tick_params(axis="x", labelbottom=False, bottom=False)
                y -= WALK_GAP
        last_ax.append(ax)
    for ax in last_ax:  # x-axis under each column (same scale)
        ax.tick_params(axis="x", labelbottom=True, bottom=True)
        position_ticks(ax, axis)
        mcp.thin_ticks(ax, axis, min_gap=0.08)
        for t in ax.get_xticklabels():
            t.set_fontsize(FS_TICK)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor=rp.HAP_COLORS["A"], label="haplotype A"),
                        Patch(facecolor=rp.HAP_COLORS["B"], label="haplotype B"),
                        Patch(facecolor="#9a9a9a", label="no allele (CoRAL, CE)")],
               loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=FS_LEGEND, frameon=False,
               handlelength=1.6, handleheight=0.8, borderaxespad=0.2)
    fig.savefig(out_stem + ".pdf")
    fig.savefig(out_stem + ".png", dpi=300)
    report = measure(fig)
    plt.close(fig)
    layout.update(report)
    return (FIG_W, height), layout, {k: sorted(set(v)) for k, v in clipped.items()}


def measure(fig):
    """Actual figure size, smallest font, and text boxes that overlap each other or leave the canvas."""
    from matplotlib.text import Text

    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    sizes = sorted({t.get_fontsize() for t in texts})
    W, H = fig.bbox.width, fig.bbox.height
    boxes = [(t.get_text(), t.get_window_extent(renderer)) for t in texts]
    outside = [s for s, b in boxes if b.x0 < -0.5 or b.y0 < -0.5 or b.x1 > W + 0.5 or b.y1 > H + 0.5]
    overlaps = []
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            a, b = boxes[i][1], boxes[j][1]
            if a.x0 < b.x1 - 0.5 and b.x0 < a.x1 - 0.5 and a.y0 < b.y1 - 0.5 and b.y0 < a.y1 - 0.5:
                overlaps.append((boxes[i][0], boxes[j][0]))
    w, h = fig.get_size_inches()
    return {"size_in": (round(float(w), 3), round(float(h), 3)), "font_sizes_pt": sizes, "n_texts": len(texts),
            "outside": outside, "overlaps": overlaps}


def readme(set_name, case, folder, rows, pairs, clipped, layout):
    count = {k: sum(1 for r in rows if r[0] == k) for k in ["truth", "RCK MC", "CoRAL", "CE"]}
    w, h = layout["size_in"]
    lines = ["{n} ({s}): two-column print-size topology comparison (make_compare_plots_ce_2col.py, 2026-11-26)".format(
                 n=case["name"], s=set_name), "",
             "{st}.pdf (vector) / {st}.png (300 dpi): the walks of topology_compare_ce.png (same data, pairing and "
             "drawing; see README.txt) re-laid out for a printed page.".format(st=STEM), "",
             "layout: two columns sharing one genome x-axis (same regions and scale) and one row height "
             "({r:.3f} in per step); left = Truth (top), RCK MC; right = CoRAL (top), CE; columns start at the same "
             "height. One row per step in walk order (reference-contiguous steps merged), numbered on the left; "
             "arrows = strand; connectors in traversal order, dotted = closing join. Label above each walk: name -> "
             "the truth species the table scorer pairs it with, bulk copies; psi per clone (Truth, RCK MC) or "
             "cycle / linear walk (CoRAL, CE). RCK MC haplotypes after the scorer's one A/B swap per phase block. "
             "One shared legend at the bottom; no figure title (goes in the caption).".format(r=layout["row_h_in"]),
             "size: {w:.2f} x {h:.2f} in (drawn at the final printed size, full 6.5 in text width; column width "
             "{c:.2f} in).".format(w=w, h=h, c=layout["col_w_in"]),
             "fonts (pt): panel titles {t:g} bold, walk labels {l:g}, step numbers {s:g}, axis ticks {k:g}, axis label "
             "and chromosome names {a:g}, legend {g:g}; smallest {m:g}. Fonts embedded as TrueType (pdf.fonttype 42).".format(
                 t=FS_TITLE, l=FS_LABEL, s=FS_STEP, k=FS_TICK, a=FS_AXIS, g=FS_LEGEND, m=min(layout["font_sizes_pt"])),
             "walks drawn: Truth {truth}, RCK MC {RCK MC}, CoRAL {CoRAL}, CE {CE}.".format(**count), "", "pairings drawn:"]
    for label in ["RCK MC", "CoRAL", "CE"]:
        lines.append("  {l}: {p}".format(l=label, p=", ".join("{a} -> {b}".format(a=a, b=b) for a, b in sorted(pairs[label].items()))
                                          or "none"))
    if clipped:
        lines.append("steps over {m} Mb clipped to the axis: {c}".format(m=mcp.LONG_STEP // 1_000_000, c="; ".join(
            "{k} {v}".format(k=k, v=", ".join(v)) for k, v in clipped.items())))
    lines += ["", "inputs (read only):"] + ["  {w}: {p}".format(w=w_, p=p) for w_, p in mce.inputs(case)]
    with open(os.path.join(folder, STEM + "_README.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")


def run(set_name, name):
    alb.use_set(set_name)
    case = next(c for c in alb.bi.cases() if c["name"] == name)
    folder = os.path.join(OUT, name)
    os.makedirs(folder, exist_ok=True)
    rows, pairs = load(case)
    work = tempfile.mkdtemp(prefix="compare_ce_2col_")
    try:
        size, layout, clipped = build(case, rows, work, os.path.join(folder, STEM))
    finally:
        shutil.rmtree(work, ignore_errors=True)
    readme(set_name, case, folder, rows, pairs, clipped, layout)
    print("{n}: size {s} in; row {r:.3f} in; rows/column {c}; fonts {f}; {t} texts; outside {o}; overlaps {v}".format(
        n=name, s=layout["size_in"], r=layout["row_h_in"], c=layout["rows_per_column"], f=layout["font_sizes_pt"],
        t=layout["n_texts"], o=layout["outside"], v=layout["overlaps"]))
    for kind in ["truth", "RCK MC", "CoRAL", "CE"]:
        print("  {k}: {w}".format(k=kind, w=[(c.name, round(c.weight, 1), pairs.get(kind, {}).get(c.name)) for k, c, _ in rows if k == kind]))
    print("  clipped {c}".format(c=clipped))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--case", nargs=2, metavar=("SET", "CASE"), help="one case (default: make_compare_plots_ce.CASES)")
    a = p.parse_args()
    import logging
    logging.disable(logging.INFO)
    for set_name, name in ([tuple(a.case)] if a.case else mce.CASES):
        run(set_name, name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
