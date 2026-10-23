"""Figures for a finished ``rck-ecdna`` run: one page per inferred species, and the recovery matrix.

Everything is read back from the run's output folder (``decompositions.tsv``,
``flow_segments.tsv``, ``summary.json`` and, after an evaluation, ``evaluation.tsv``,
``truth_cycles.tsv`` and ``recovery_matrix.tsv``), so ``rck-ecdna-plot`` can redraw any earlier run.

Each species page has, on one cnplot ``GenomeAxis`` restricted to the segments the species (and
its matched truth species) uses:

* the topology: one row per step of the top ordering, a block at the segment's reference
  position coloured by haplotype, an arrow for the strand, and connectors in traversal order
  (CoRAL's cycles figure; the same track amplicon-allele-analysis draws for simulated truth);
* the matched truth species, drawn the same way, when the run was evaluated;
* the junctions the walk takes, as arcs coloured by haplotype (grey for a haplotype crossing,
  dotted for an imputed junction);
* allele-specific copy number, A above and B mirrored below: the observations as points
  (cnplot ``plot_scatter_1d``), the inferred ecDNA + background as a step line and this species'
  own copies (weight x multiplicity) as shaded blocks;

and, beside them, the walk as a ring, step lengths to scale except that every step gets at least
``MIN_ARC_FRACTION`` of the circle. A linear species is drawn as an open ring.

cnplot (amplicon-allele-analysis ``third_party/cnplot``) and matplotlib are optional: without
them the figures are skipped with a warning. cnplot is looked up on the normal import path,
then under ``--cnplot-path`` / ``$CNPLOT_PATH`` (the cnplot checkout or its ``src`` folder).
"""

from __future__ import annotations

import argparse
import importlib
import json
import logging
import math
import os
import sys
import tempfile
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

HAP_COLORS = {"A": "#2a78d6", "B": "#eb6834"}  # categorical slots 1 and 2: blue, orange
CROSSING_COLOR = "#52514e"
INK = "#0b0b0b"
INK_MUTED = "#52514e"
MIN_ARC_FRACTION = 0.015
PAD_FRACTION = 0.05
MIN_PAD_BP = 5_000
POINTS_PER_PAGE = 400
MATRIX_COLUMNS = [("segment_jaccard", "segment\nJaccard"), ("haplotype_jaccard", "haplotype\nJaccard"),
                  ("junction_f1", "junction\nF1"), ("phased_junction_f1", "phased\njunction F1"),
                  ("phasing_accuracy", "phasing\naccuracy"), ("order_match", "order\nmatch"), ("cn_accuracy", "1 - CN\nrel. err"),
                  ("composite", "composite"), ("best_composite", "best\ncomposite")]


@dataclass(frozen=True)
class PStep:
    chrom: str
    start: int  # inclusive
    end: int  # inclusive
    hap: str
    strand: str


@dataclass
class PSpecies:
    rank: int
    name: str
    topology: str
    status: str
    weight: float
    weight_int: Optional[float]
    length_bp: int
    call: str
    steps: List[PStep]
    joins: List[str]  # kind:aid:phasing after each step (end:<extremity>:<hap> closes a linear walk)


# ---------------------------------------------------------------------------------------- #
# loading
# ---------------------------------------------------------------------------------------- #

def import_plotting(cnplot_path: Optional[str] = None):
    """(pyplot, PdfPages, cnplot), or None when matplotlib or cnplot cannot be imported."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError as error:
        log.warning("plots skipped: %s", error)
        return None
    candidates = [p for p in (cnplot_path, os.environ.get("CNPLOT_PATH")) if p]
    try:
        return plt, PdfPages, importlib.import_module("cnplot")
    except ImportError as error:
        first_error = error
    for path in candidates:
        for entry in (os.path.join(path, "src"), path):
            if os.path.isdir(os.path.join(entry, "cnplot")) and entry not in sys.path:
                sys.path.insert(0, entry)
                try:
                    return plt, PdfPages, importlib.import_module("cnplot")
                except ImportError as error:
                    first_error = error
    log.warning("plots skipped: cnplot unavailable (%s); pip install -e amplicon-allele-analysis/third_party/cnplot "
                "or pass --cnplot-path", first_error)
    return None


def _parse_step(token: str) -> PStep:
    chrom, span, hap, strand = token.rsplit(":", 3)
    start, end = span.split("-")
    return PStep(chrom, int(start), int(end), hap, strand)


def read_species(decompositions_tsv: str, ranks: int = 1, include_low_weight: bool = False) -> List[PSpecies]:
    """Top ordering of every species in the first ``ranks`` decompositions."""
    frame = pd.read_csv(decompositions_tsv, sep="\t")
    if not len(frame):
        return []
    frame = frame[frame["decomposition_rank"] <= frame["decomposition_rank"].min() + ranks - 1]
    if "status" in frame and not include_low_weight:
        frame = frame[frame["status"] != "low_weight"]
    result = []
    for (rank, name), group in frame.groupby(["decomposition_rank", "species"], sort=True):
        row = group.sort_values("ordering_rank").iloc[0]
        weight_int = row.get("weight_int")
        result.append(PSpecies(rank=int(rank), name=str(name), topology=str(row.get("topology", "circular")),
                               status=str(row.get("status", "species")), weight=float(row["weight"]),
                               weight_int=None if weight_int is None or pd.isna(weight_int) else float(weight_int),
                               length_bp=int(row["length_bp"]), call=str(row.get("amplicon_call", "")),
                               steps=[_parse_step(t) for t in str(row["walk"]).split(",")],
                               joins=str(row["joins"]).split(",") if "joins" in row and pd.notna(row["joins"]) else []))
    return result


def read_truth_steps(path: str) -> Dict[str, Tuple[float, List[PStep]]]:
    """{species: (weight, steps)} from ``truth_cycles.tsv``."""
    frame = pd.read_csv(path, sep="\t")
    result = {}
    for name, group in frame.groupby("species", sort=False):
        group = group.sort_values("order")
        result[str(name)] = (float(group["weight"].iloc[0]),
                             [PStep(str(r["chrom"]), int(r["start"]), int(r["end"]), str(r["hap"]), str(r["strand"]))
                              for r in group.to_dict("records")])
    return result


def _read_optional(path: str) -> Optional[pd.DataFrame]:
    return pd.read_csv(path, sep="\t") if os.path.exists(path) and os.path.getsize(path) > 0 else None


# ---------------------------------------------------------------------------------------- #
# the genome axis
# ---------------------------------------------------------------------------------------- #

def axis_regions(steps: Sequence[PStep]) -> List[Tuple[str, int, int]]:
    """Half-open intervals covering ``steps``, padded and merged per chromosome, in first-use order."""
    by_chrom: Dict[str, List[Tuple[int, int]]] = {}
    for s in steps:
        by_chrom.setdefault(s.chrom, []).append((s.start, s.end + 1))
    regions = []
    for chrom, spans in by_chrom.items():
        spans.sort()
        merged: List[List[int]] = []
        for a, b in spans:
            pad = max(MIN_PAD_BP, int(PAD_FRACTION * (b - a)))
            a, b = max(a - pad, 0), b + pad
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        regions.extend((chrom, a, b) for a, b in merged)
    return regions


def build_axis(cnplot, regions: Sequence[Tuple[str, int, int]], scratch: str):
    """A cnplot GenomeAxis over ``regions`` with the gaps between them collapsed.

    GenomeAxis needs a chromosome-sizes file; the run does not know the reference, so each
    chromosome is given the length of its last region, which is all the axis uses.
    """
    bed = os.path.join(scratch, "regions.bed")
    sizes = os.path.join(scratch, "chrom.sizes")
    with open(bed, "w") as handle:
        handle.writelines("{c}\t{a}\t{b}\tr\n".format(c=c, a=a, b=b) for c, a, b in regions)
    ends: Dict[str, int] = {}
    for c, _, b in regions:
        ends[c] = max(ends.get(c, 0), b)
    with open(sizes, "w") as handle:
        handle.writelines("{c}\t{n}\n".format(c=c, n=n) for c, n in ends.items())
    return cnplot.GenomeAxis(bed, sizes, excluded_chroms=[], collapse_gaps=True)


def _intervals(axis, steps: Sequence[PStep]) -> Tuple[np.ndarray, np.ndarray]:
    frame = pd.DataFrame({"#CHR": [s.chrom for s in steps], "START": [s.start for s in steps], "END": [s.end + 1 for s in steps]})
    coords = axis.build_coordinates(frame)
    return np.asarray(coords.starts, dtype=float), np.asarray(coords.ends, dtype=float)


def _nice_step(span: float, target: int = 6) -> int:
    raw = max(span / target, 1.0)
    magnitude = 10 ** math.floor(math.log10(raw))
    for factor in (1, 2, 5, 10):
        if factor * magnitude >= raw:
            return int(factor * magnitude)
    return int(10 * magnitude)


def position_ticks(ax, axis) -> None:
    """Mb position ticks inside each axis segment (spaced for that segment, so a short one still gets
    some), and the chromosome names under them."""
    positions, labels = [], []
    for seg in axis.segments:
        shift = seg.axis_start - seg.raw_start
        step = _nice_step(seg.raw_end - seg.raw_start, target=3)
        first = int(math.ceil(seg.raw_start / step)) * step
        for p in range(first, int(seg.raw_end), step):
            positions.append(p + shift)
            labels.append("{v:g}".format(v=round(p / 1e6, 3)))
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right", rotation_mode="anchor")
    ax.set_xlabel("position (Mb)", fontsize=8, color=INK_MUTED)
    for x, name in zip(axis.xtick_chrs, axis.xlab_chrs):
        ax.text(x, -0.45, str(name), ha="center", va="top", fontsize=8, color=INK, transform=ax.get_xaxis_transform())


# ---------------------------------------------------------------------------------------- #
# tracks
# ---------------------------------------------------------------------------------------- #

def _join_kind(join: str) -> Tuple[str, str]:
    parts = join.split(":")
    return parts[0], parts[-1] if len(parts) > 2 else ""


def draw_topology(ax, axis, steps: Sequence[PStep], circular: bool, label: str, joins: Sequence[str] = ()) -> None:
    """One row per step: a block at its reference position coloured by haplotype, an arrow for the
    strand, and connectors in traversal order (dotted for the join that closes a circle)."""
    from matplotlib.patches import Rectangle

    lows, highs = _intervals(axis, steps)
    count = len(steps)
    for i, s in enumerate(steps):
        low, high = lows[i], highs[i]
        if not (np.isfinite(low) and np.isfinite(high)):
            continue
        y = count - i - 0.5
        ax.add_patch(Rectangle((low, y - 0.28), max(high - low, 1.0), 0.56, facecolor=HAP_COLORS.get(s.hap, "#9a9a9a"),
                               edgecolor="white", linewidth=0.6, zorder=3))
        forward = s.strand == "+"
        ax.annotate("", xy=(high if forward else low, y), xytext=(low if forward else high, y),
                    arrowprops={"arrowstyle": "-|>", "color": "white", "lw": 0.8, "shrinkA": 0, "shrinkB": 0}, zorder=4)
        if count <= 40:
            ax.text(axis.segments[0].axis_start, y, str(i + 1), ha="right", va="center", fontsize=6, color=INK_MUTED,
                    transform=ax.transData, clip_on=False)
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
        kind = _join_kind(joins[i])[0] if i < len(joins) else ""
        ax.plot([exit_x, exit_x, entry_x, entry_x], [y0, mid, mid, y1], color=INK_MUTED, linewidth=0.8,
                linestyle=":" if (nxt == 0 or kind == "imputed") else "-", zorder=2)
    ax.set_ylim(-0.4, count + 0.1)
    ax.set_yticks([])
    ax.set_ylabel(label, fontsize=8, color=INK, labelpad=14)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def draw_junctions(ax, axis, steps: Sequence[PStep], circular: bool, joins: Sequence[str]) -> None:
    """Arcs for the non-reference joins of the walk, coloured by the haplotype they connect."""
    from matplotlib.patches import Arc

    lows, highs = _intervals(axis, steps)
    count = len(steps)
    span = axis.chr_end - axis.ch_coords[0]
    top = 0.0
    for i in range(count if circular else count - 1):
        kind, phasing = _join_kind(joins[i]) if i < len(joins) else ("novel", "")
        if kind in ("reference", "end"):
            continue
        nxt = (i + 1) % count
        a = highs[i] if steps[i].strand == "+" else lows[i]
        b = lows[nxt] if steps[nxt].strand == "+" else highs[nxt]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        width = max(abs(b - a), span * 0.004)
        height = 0.25 + 0.75 * min(width / span, 1.0)
        top = max(top, height)
        crossing = len(phasing) == 2 and phasing[0] != phasing[1]
        color = CROSSING_COLOR if crossing or not phasing else HAP_COLORS.get(phasing[0], CROSSING_COLOR)
        ax.add_patch(Arc(((a + b) / 2, 0.0), width, 2 * height, theta1=0, theta2=180, color=color, linewidth=1.4,
                         linestyle=":" if kind == "imputed" else "-"))
    ax.set_ylim(0, max(top, 0.3) * 1.1)
    ax.set_yticks([])
    ax.set_ylabel("junctions", fontsize=8, color=INK)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)


def _species_copies(segments: pd.DataFrame, species: PSpecies) -> Dict[Tuple[str, str], float]:
    """{(sid, hap): copies this species puts on the segment} from the walk's step multiplicities."""
    index = {(r["chrom"], int(r["start"]), int(r["end"])): r["sid"] for r in segments.to_dict("records")}
    copies: Dict[Tuple[str, str], float] = {}
    for s in species.steps:
        sid = index.get((s.chrom, s.start, s.end))
        if sid is not None:
            copies[(sid, s.hap)] = copies.get((sid, s.hap), 0.0) + species.weight
    return copies


def _pseudo_bins(segments: pd.DataFrame, regions: Sequence[Tuple[str, int, int]]) -> pd.DataFrame:
    """Split each segment, clipped to the axis regions, into bins of about 1 / POINTS_PER_PAGE of the axis
    so a segment-level observation draws as a run of points along it, like the binned tracks of the
    simulator figures."""
    size = max(int(sum(b - a for _, a, b in regions) // POINTS_PER_PAGE), 1)
    rows = []
    for r in segments.to_dict("records"):
        for chrom, lo, hi in regions:
            if chrom != r["chrom"]:
                continue
            start, end = max(int(r["start"]), lo), min(int(r["end"]) + 1, hi)
            for a in range(start, end, size):
                rows.append({"#CHR": chrom, "START": a, "END": min(a + size, end), "CN-A": r["obs_A"], "CN-B": r["obs_B"],
                             "CN": r["obs_total"]})
    return pd.DataFrame(rows, columns=["#CHR", "START", "END", "CN-A", "CN-B", "CN"])


def draw_copy_number(ax, axis, cnplot, segments: pd.DataFrame, species: PSpecies, hap: str, ylim: Tuple[float, float],
                     regions: Sequence[Tuple[str, int, int]]) -> None:
    """Observed allele copy number (points), inferred ecDNA + background (step line) and this species'
    copies (shaded) for haplotype ``hap``; B is drawn mirrored under A."""
    from matplotlib.patches import Rectangle

    column = "CN-" + hap
    bins = _pseudo_bins(segments, regions)
    expected = pd.DataFrame({"#CHR": segments["chrom"], "START": segments["start"].astype(int), "END": segments["end"].astype(int) + 1,
                             "exp_" + column: segments["ecdna_" + hap] + segments["bg_" + hap]})
    if bins[column].notna().any():
        cnplot.plot_scatter_1d(ax, bins, axis, ycol=column, expected_df=expected, colors=[HAP_COLORS[hap]] * len(bins),
                               markersize=6, ylim=ylim, ylabel=column, reverse_y=hap == "B", plot_chrname=False,
                               show_legend=False, exp_color=INK, exp_linewidth=1.2, href=None)
    else:  # no allele-specific observation anywhere on the axis: only the inferred line
        cnplot.plot_scatter_1d(ax, expected.assign(**{column: np.nan}), axis, ycol=column, expected_df=expected, ylim=ylim,
                               ylabel=column, reverse_y=hap == "B", plot_chrname=False, show_legend=False, exp_color=INK)
    lows, highs = _intervals(axis, [PStep(r["chrom"], int(r["start"]), int(r["end"]), hap, "+") for r in segments.to_dict("records")])
    copies = _species_copies(segments, species)
    for i, sid in enumerate(segments["sid"]):
        value = copies.get((sid, hap), 0.0)
        if value > 0 and np.isfinite(lows[i]) and np.isfinite(highs[i]):
            ax.add_patch(Rectangle((lows[i], 0.0), max(highs[i] - lows[i], 1.0), value, facecolor=HAP_COLORS[hap], alpha=0.25,
                                   edgecolor="none", zorder=1))
    ax.set_ylabel(column, fontsize=8, color=INK)
    ax.tick_params(axis="y", labelsize=7)


def draw_ring(ax, species: PSpecies) -> None:
    """The walk as a ring: one arc per step, coloured by haplotype, clockwise in traversal order, with a
    tick at every non-reference join. A linear walk leaves a gap between its two ends."""
    from matplotlib.patches import Wedge

    steps = species.steps
    lengths = np.array([s.end - s.start + 1 for s in steps], dtype=float)
    gap = 0.06 if species.topology == "linear" else 0.0
    shares = lengths / lengths.sum()
    shares = np.maximum(shares, MIN_ARC_FRACTION)
    shares = shares / shares.sum() * (1.0 - gap)
    floored = bool((lengths / lengths.sum() < MIN_ARC_FRACTION).any())
    angle = 90.0 - 360.0 * gap / 2  # start at the top and run clockwise
    label_every = max(1, int(math.ceil(len(steps) / 24)))
    for i, (s, share) in enumerate(zip(steps, shares)):
        sweep = 360.0 * share
        ax.add_patch(Wedge((0, 0), 1.0, angle - sweep, angle, width=0.22, facecolor=HAP_COLORS.get(s.hap, "#9a9a9a"),
                           edgecolor="white", linewidth=1.0))
        middle = math.radians(angle - sweep / 2)
        if i % label_every == 0:
            ax.text(1.14 * math.cos(middle), 1.14 * math.sin(middle), str(i + 1), ha="center", va="center", fontsize=6, color=INK_MUTED)
        # strand: a chevron inside the arc pointing along (+) or against (-) the reading direction
        sign = -1 if s.strand == "+" else 1
        tip = math.radians(angle - sweep / 2 + sign * min(sweep / 4, 4))
        ax.annotate("", xy=(0.89 * math.cos(tip), 0.89 * math.sin(tip)), xytext=(0.89 * math.cos(middle), 0.89 * math.sin(middle)),
                    arrowprops={"arrowstyle": "-|>", "color": "white", "lw": 0.6, "mutation_scale": 6})
        end_angle = math.radians(angle - sweep)
        kind = _join_kind(species.joins[i])[0] if i < len(species.joins) else ""
        if kind in ("novel", "imputed") and (species.topology != "linear" or i < len(steps) - 1):
            ax.plot([0.74 * math.cos(end_angle), 1.04 * math.cos(end_angle)], [0.74 * math.sin(end_angle), 1.04 * math.sin(end_angle)],
                    color=INK, linewidth=1.0, linestyle=":" if kind == "imputed" else "-")
        angle -= sweep
    ax.text(0, 0.08, "{t}".format(t=species.topology), ha="center", va="center", fontsize=9, color=INK, fontweight="bold")
    ax.text(0, -0.10, "{w:.2f} copies\n{l:.2f} Mb".format(w=species.weight, l=species.length_bp / 1e6), ha="center", va="center",
            fontsize=7, color=INK_MUTED)
    if floored:
        ax.text(0, -1.32, "short steps widened to {p:.1%} of the ring".format(p=MIN_ARC_FRACTION), ha="center", fontsize=6,
                color=INK_MUTED)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.4, 1.3)
    ax.set_aspect("equal")
    ax.axis("off")


def _legend(fig) -> None:
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    handles = [Patch(facecolor=HAP_COLORS["A"], label="haplotype A"), Patch(facecolor=HAP_COLORS["B"], label="haplotype B"),
               Line2D([], [], color=CROSSING_COLOR, label="haplotype crossing"),
               Line2D([], [], color=INK_MUTED, linestyle=":", label="imputed / closing join"),
               Line2D([], [], color=INK, label="inferred CN (ecDNA + background)"),
               Line2D([], [], marker="o", linestyle="", color=INK_MUTED, markersize=3, label="observed CN")]
    fig.legend(handles=handles, loc="lower center", ncol=6, fontsize=7, frameon=False, bbox_to_anchor=(0.5, -0.01))


# ---------------------------------------------------------------------------------------- #
# pages
# ---------------------------------------------------------------------------------------- #

def species_page(plt, cnplot, species: PSpecies, segments: pd.DataFrame, name: str, scratch: str,
                 truth: Optional[Tuple[str, float, List[PStep]]] = None, score: Optional[Dict[str, object]] = None):
    """One figure for one inferred species (see the module docstring for the tracks)."""
    all_steps = list(species.steps) + (list(truth[2]) if truth else [])
    regions = axis_regions(all_steps)
    axis = build_axis(cnplot, regions, scratch)
    keep = np.zeros(len(segments), dtype=bool)
    for c, a, b in regions:
        keep |= (segments["chrom"] == c).to_numpy() & (segments["end"].to_numpy() >= a) & (segments["start"].to_numpy() < b)
    local = segments[keep].reset_index(drop=True)

    n = len(species.steps)
    heights = [max(1.6, 0.2 * n)]
    if truth:
        heights.append(max(1.2, 0.2 * len(truth[2])))
    heights += [1.2, 1.5, 1.5]
    height = min(sum(heights) + 1.8, 30)
    fig = plt.figure(figsize=(15, height))
    grid = fig.add_gridspec(len(heights), 2, width_ratios=[3.2, 1.0], height_ratios=heights, hspace=0.18, wspace=0.06,
                            top=1.0 - 0.75 / height, bottom=0.9 / height)
    axes = [fig.add_subplot(grid[0, 0])]
    for row in range(1, len(heights)):
        axes.append(fig.add_subplot(grid[row, 0], sharex=axes[0]))
    row = 0
    draw_topology(axes[row], axis, species.steps, species.topology != "linear", "inferred\n{s}".format(s=species.name), species.joins)
    if truth:
        row += 1
        draw_topology(axes[row], axis, truth[2], True, "truth\n{t}".format(t=truth[0]))
    row += 1
    draw_junctions(axes[row], axis, species.steps, species.topology != "linear", species.joins)
    observed = pd.concat([local["obs_A"], local["obs_B"], local["ecdna_A"] + local["bg_A"], local["ecdna_B"] + local["bg_B"]])
    peak = float(np.nanmax(observed.to_numpy(dtype=float))) if len(local) and np.isfinite(observed.to_numpy(dtype=float)).any() else 1.0
    ylim = (0.0, max(peak * 1.1, 2.0))
    draw_copy_number(axes[row + 1], axis, cnplot, local, species, "A", ylim, regions)
    draw_copy_number(axes[row + 2], axis, cnplot, local, species, "B", ylim, regions)
    for ax in axes:
        ax.set_xlim(axis.ch_coords[0], axis.chr_end)
        for boundary in axis.chr_boundaries:
            ax.axvline(boundary, color=INK_MUTED, linewidth=0.8, zorder=0)
        ax.tick_params(axis="x", labelbottom=False, bottom=False)
    axes[-1].tick_params(axis="x", labelbottom=True, bottom=True)
    position_ticks(axes[-1], axis)
    draw_ring(fig.add_subplot(grid[:, 1]), species)

    weight = "{w:.2f} copies".format(w=species.weight) + ("" if species.weight_int is None else " (integer {i:g})".format(i=species.weight_int))
    title = "{n}  |  decomposition {r} ({c})  |  {s}: {t}, {w}, {l:,} bp, {k} steps{f}".format(
        n=name, r=species.rank, c=species.call, s=species.name, t=species.topology, w=weight, l=species.length_bp, k=n,
        f=", low_weight" if species.status == "low_weight" else "")
    if truth:
        title += "\nmatched truth {t} ({w:.2f} copies)".format(t=truth[0], w=truth[1])
        if score:
            title += ": composite {c:.3f}, segment Jaccard {j:.3f}, junction F1 {f:.3f}, phasing accuracy {p:.3f}, order {o}".format(
                c=score["composite"], j=score["jaccard_total"], f=score["junction_f1"], p=score["phasing_accuracy"],
                o="matched" if bool(score["order_match"]) else "not matched")
    fig.suptitle(title, fontsize=10, color=INK, x=0.02, y=1.0 - 0.1 / height, ha="left", va="top")
    _legend(fig)
    return fig


def matrix_page(plt, matrix: pd.DataFrame, name: str):
    """Heatmap of the recovery matrix: one row per in-scope truth species, every score in [0, 1]."""
    rows = matrix[matrix["in_scope"].astype(bool)].copy()
    rows["cn_accuracy"] = (1.0 - pd.to_numeric(rows["cn_rel_err"], errors="coerce")).clip(lower=0.0)
    values = np.array([[_as_float(r.get(c)) for c, _ in MATRIX_COLUMNS] for r in rows.to_dict("records")], dtype=float)
    fig, ax = plt.subplots(figsize=(1.0 * len(MATRIX_COLUMNS) + 3.0, 0.45 * len(rows) + 1.6))
    image = ax.imshow(np.nan_to_num(values, nan=0.0), cmap="Blues", vmin=0.0, vmax=1.0, aspect="auto")
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            v = values[i, j]
            ax.text(j, i, "-" if np.isnan(v) else "{v:.2f}".format(v=v), ha="center", va="center", fontsize=8,
                    color="white" if (not np.isnan(v) and v > 0.6) else INK)
    ax.set_xticks(range(len(MATRIX_COLUMNS)))
    ax.set_xticklabels([label for _, label in MATRIX_COLUMNS], fontsize=8)
    labels = ["{t} -> {i}".format(t=r["truth"], i=r["inferred"] if isinstance(r["inferred"], str) else "missed") for r in rows.to_dict("records")]
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    mean = pd.to_numeric(rows["composite"], errors="coerce").mean()
    ax.set_title("{n}: recovery of the truth species (rank-1 decomposition; best = best over all alternatives)\n"
                 "mean composite {m}".format(n=name, m="n/a" if pd.isna(mean) else "{v:.3f}".format(v=mean)),
                 fontsize=10, color=INK, loc="left")
    fig.colorbar(image, ax=ax, fraction=0.03, pad=0.02).ax.tick_params(labelsize=7)
    fig.tight_layout()
    return fig


def _as_float(value) -> float:
    if value is None or (isinstance(value, str) and not value):
        return float("nan")
    if isinstance(value, str):
        return {"true": 1.0, "false": 0.0}.get(value.lower(), float("nan"))
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def plot_run(outdir: str, ranks: int = 1, include_low_weight: bool = False, cnplot_path: Optional[str] = None,
             png: bool = True) -> Dict[str, object]:
    """Write ``plots/<name>_ecdna.pdf`` (one page per species, then the recovery matrix when there is
    one) and, with ``png``, one PNG per page. Returns {"pdf": path, "png": [paths]}, or {} when
    plotting is unavailable."""
    imported = import_plotting(cnplot_path)
    if imported is None:
        return {}
    plt, PdfPages, cnplot = imported
    name = os.path.basename(os.path.normpath(outdir))
    summary_path = os.path.join(outdir, "summary.json")
    if os.path.exists(summary_path):
        with open(summary_path) as handle:
            name = json.load(handle).get("graph", {}).get("name", name)
    species = read_species(os.path.join(outdir, "decompositions.tsv"), ranks=ranks, include_low_weight=include_low_weight)
    segments = pd.read_csv(os.path.join(outdir, "flow_segments.tsv"), sep="\t")
    evaluation = _read_optional(os.path.join(outdir, "evaluation.tsv"))
    truth_path = os.path.join(outdir, "truth_cycles.tsv")
    truth = read_truth_steps(truth_path) if os.path.exists(truth_path) else {}
    matrix = _read_optional(os.path.join(outdir, "recovery_matrix.tsv"))

    plot_dir = os.path.join(outdir, "plots")
    os.makedirs(plot_dir, exist_ok=True)
    pdf_path = os.path.join(plot_dir, "{n}_ecdna.pdf".format(n=name))
    written: Dict[str, object] = {"pdf": pdf_path, "png": []}
    with PdfPages(pdf_path) as pdf, tempfile.TemporaryDirectory() as scratch:
        for sp in species:
            match, score = None, None
            if evaluation is not None and len(evaluation) and "inferred" in evaluation:
                hit = evaluation[(evaluation["decomposition_rank"] == sp.rank) & (evaluation["inferred"] == sp.name)]
                if len(hit) and str(hit["truth"].iloc[0]) in truth:
                    tname = str(hit["truth"].iloc[0])
                    match = (tname, truth[tname][0], truth[tname][1])
                    score = hit.iloc[0].to_dict()
                    if "composite" not in score:
                        score = None  # evaluation written before the composite score existed
            fig = species_page(plt, cnplot, sp, segments, name, scratch, truth=match, score=score)
            pdf.savefig(fig, bbox_inches="tight")
            if png:
                path = os.path.join(plot_dir, "{n}_rank{r}_{s}.png".format(n=name, r=sp.rank, s=sp.name))
                fig.savefig(path, dpi=130, bbox_inches="tight")
                written["png"].append(path)
            plt.close(fig)
        if matrix is not None and matrix["in_scope"].astype(bool).any():
            fig = matrix_page(plt, matrix, name)
            pdf.savefig(fig, bbox_inches="tight")
            if png:
                path = os.path.join(plot_dir, "{n}_recovery_matrix.png".format(n=name))
                fig.savefig(path, dpi=130, bbox_inches="tight")
                written["png"].append(path)
            plt.close(fig)
    log.info("plots: %d species page(s)%s in %s", len(species), " and the recovery matrix" if matrix is not None else "", pdf_path)
    return written


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="rck-ecdna-plot", description="Draw the species and recovery matrix of an rck-ecdna run.")
    parser.add_argument("--outdir", required=True, help="rck-ecdna output folder (with decompositions.tsv and flow_segments.tsv)")
    parser.add_argument("--plot-ranks", type=int, default=1, help="Draw the species of this many top decompositions")
    parser.add_argument("--plot-low-weight", action="store_true", help="Also draw species flagged low_weight")
    parser.add_argument("--cnplot-path", default=None, help="cnplot checkout (or its src/) when cnplot is not installed")
    parser.add_argument("--no-png", action="store_true", help="Write only the PDF")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    written = plot_run(args.outdir, ranks=args.plot_ranks, include_low_weight=args.plot_low_weight, cnplot_path=args.cnplot_path,
                       png=not args.no_png)
    return 0 if written else 1


if __name__ == "__main__":
    raise SystemExit(main())
