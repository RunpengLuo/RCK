"""Convert read depth and junction read counts into copy numbers.

Two scales are needed:

* ``depth_per_copy`` (gamma): tumor read depth contributed by one copy of a locus, averaged
  over the bulk sample. Segment copy number is ``mean tumor depth / gamma``.
* ``junction_reads_per_copy`` (theta): reads spanning one copy of a junction. Junction copy
  number is ``#NOV / theta`` (novel) or ``#REF / theta`` (reference join).

Gamma is taken from the first available source:

1. ``--depth-per-copy`` (explicit override);
2. ``--hatchet-seg`` (HATCHet ``best.seg.ucn``): median of tumor depth / fractional total copy
   number over long background segments. This is the route for cell lines (no matched
   normal) and for tumors whose dominant copy number is not 2;
3. the tumor windows: aggregate 1 kb windows into ``block_bp`` blocks, keep blocks away from
   the amplicon whose mirrored B-allele frequency says they are allelically balanced, take
   the mode of their depth and call it ``baseline_cn`` (default 2).

``depth_mode='ratio'`` divides tumor by matched-normal depth window by window before the same
mode calibration. That removes coverage bias that is shared between the two libraries (GC,
mappability), which is the usual situation for real data. In the amplicon-allele-analysis
simulator the bias field is drawn independently per sample, so the ratio only adds noise
there; ``tumor`` is the default for that reason.

Allele-specific copy numbers use the phase-oriented segment BAF already computed by the
subgraph pipeline: ``c_A = c (1 - BAF)``, ``c_B = c BAF``. Segments with too few SNPs get a
total-only observation.
"""

from __future__ import annotations

import json
import logging
import math
import os
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from rck.ecdna.subgraph_io import ECGraph, REFERENCE

log = logging.getLogger(__name__)


@dataclass
class ScaleParams:
    depth_per_copy: Optional[float] = None
    junction_reads_per_copy: Optional[float] = None
    baseline_cn: float = 2.0
    depth_mode: str = "tumor"  # tumor | ratio
    sample: Optional[str] = None
    hatchet_seg: Optional[str] = None
    block_bp: int = 1_000_000
    balanced_min_mbaf: float = 0.40
    min_block_snps: int = 20
    exclude_pad_bp: int = 1_000_000
    flank_bp: int = 1_000_000
    flank_gap_bp: int = 10_000
    min_snps_allele: int = 5
    theta_window_bp: int = 25_000
    hatchet_min_segment_bp: int = 1_000_000


@dataclass
class Scale:
    depth_per_copy: float
    junction_reads_per_copy: float
    baseline_cn: float
    depth_method: str
    junction_method: str
    depth_mode: str = "tumor"
    ratio_baseline: Optional[float] = None
    details: Dict[str, object] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        return asdict(self)


@dataclass
class SegmentObservation:
    sid: str
    depth: float
    total: float
    a: Optional[float]
    b: Optional[float]
    n_bins: float
    source: str


@dataclass
class Observations:
    scale: Scale
    segments: Dict[str, SegmentObservation]
    adjacencies: Dict[str, Optional[float]]
    background: Dict[str, float]

    def write_json(self, path: str) -> None:
        payload = {
            "scale": self.scale.to_dict(),
            "segments": {k: asdict(v) for k, v in self.segments.items()},
            "adjacencies": self.adjacencies,
            "background_total_cn": self.background,
        }
        with open(path, "w") as handle:
            json.dump(payload, handle, indent=2, default=float)


# ---------------------------------------------------------------------------------------- #
# windows
# ---------------------------------------------------------------------------------------- #

class Windows:
    """1 kb observation windows in Universal-Genotyping-Pipeline layout.

    Expects ``window.tsv.gz`` (``#CHR START END``), ``window.depth.npz`` (``mat``: windows x
    samples) and ``sample_ids.tsv`` (one row per ``mat`` column, with ``sample_type`` normal or
    tumor). SNP files (``snp.tsv.gz``, ``snp.Aallele.npz``, ``snp.Ballele.npz``) are optional
    and only used to find allelically balanced blocks.
    """

    def __init__(self, windows_dir: str, sample: Optional[str] = None):
        self.dir = windows_dir
        table = pd.read_csv(os.path.join(windows_dir, "window.tsv.gz"), sep="\t", usecols=["#CHR", "START", "END"])
        self.depth = np.load(os.path.join(windows_dir, "window.depth.npz"))["mat"].astype(float)
        if self.depth.ndim == 1:
            self.depth = self.depth[:, None]
        ids = pd.read_csv(os.path.join(windows_dir, "sample_ids.tsv"), sep="\t")
        if len(ids) != self.depth.shape[1]:
            raise ValueError("sample_ids.tsv has {n} rows but window.depth.npz has {m} columns".format(n=len(ids), m=self.depth.shape[1]))
        types = ids.get("sample_type", pd.Series(["tumor"] * len(ids))).astype(str).str.lower().tolist()
        names = [set(str(v) for v in row.values()) for row in ids.to_dict("records")]
        self.normal_col = types.index("normal") if "normal" in types else None
        if sample is not None:
            matches = [i for i, n in enumerate(names) if sample in n]
            if not matches:
                raise ValueError("sample {s} not found in sample_ids.tsv".format(s=sample))
            self.tumor_col = matches[0]
        else:
            tumors = [i for i, t in enumerate(types) if t == "tumor"]
            if not tumors:
                raise ValueError("no tumor sample in sample_ids.tsv")
            self.tumor_col = tumors[0]
        self.sample_names = ids.iloc[:, 0].astype(str).tolist()
        self.table = table
        self.bin_size = float(np.median((table["END"] - table["START"]).values)) if len(table) else 1000.0
        self._by_chrom = {}
        for chrom, sub in table.groupby("#CHR", sort=False):
            order = np.argsort(sub["START"].values)
            idx = sub.index.values[order]
            self._by_chrom[str(chrom)] = (sub["START"].values[order], sub["END"].values[order], idx)
        self._snps = None

    def _overlap(self, chrom: str, start: int, end_inclusive: int) -> Tuple[np.ndarray, np.ndarray]:
        """Row indices and overlap bp of windows intersecting [start, end_inclusive]."""
        if chrom not in self._by_chrom:
            return np.array([], dtype=int), np.array([])
        starts, ends, idx = self._by_chrom[chrom]
        lo = int(np.searchsorted(ends, start, side="right"))  # first window ending after start
        hi = np.searchsorted(starts, end_inclusive + 1, side="left")
        if hi <= lo:
            return np.array([], dtype=int), np.array([])
        s, e = starts[lo:hi], ends[lo:hi]
        ov = np.minimum(e, end_inclusive + 1) - np.maximum(s, start)
        keep = ov > 0
        return idx[lo:hi][keep], ov[keep].astype(float)

    def mean_depth(self, chrom: str, start: int, end_inclusive: int, col: Optional[int] = None) -> Tuple[float, float]:
        """Overlap-weighted mean depth and the effective number of windows."""
        col = self.tumor_col if col is None else col
        rows, ov = self._overlap(chrom, start, end_inclusive)
        if len(rows) == 0:
            return float("nan"), 0.0
        d = self.depth[rows, col]
        ok = np.isfinite(d)
        if not ok.any():
            return float("nan"), 0.0
        return float(np.sum(d[ok] * ov[ok]) / np.sum(ov[ok])), float(np.sum(ov[ok]) / self.bin_size)

    def ratio(self, chrom: str, start: int, end_inclusive: int) -> Tuple[float, float]:
        """Tumor / normal depth ratio over a span (sums, not mean of ratios)."""
        if self.normal_col is None:
            raise ValueError("depth_mode='ratio' needs a normal sample in sample_ids.tsv")
        rows, ov = self._overlap(chrom, start, end_inclusive)
        if len(rows) == 0:
            return float("nan"), 0.0
        t = self.depth[rows, self.tumor_col] * ov
        n = self.depth[rows, self.normal_col] * ov
        ok = np.isfinite(t) & np.isfinite(n)
        if not ok.any() or n[ok].sum() <= 0:
            return float("nan"), 0.0
        return float(t[ok].sum() / n[ok].sum()), float(ov[ok].sum() / self.bin_size)

    # ---- blocks ---- #
    def block_table(self, block_bp: int) -> pd.DataFrame:
        """Per-block mean tumor (and normal) depth, window count, and mirrored tumor BAF."""
        frame = pd.DataFrame({
            "chrom": self.table["#CHR"].astype(str).values,
            "block": (self.table["START"].values // block_bp).astype(int),
            "tumor": self.depth[:, self.tumor_col],
        })
        if self.normal_col is not None:
            frame["normal"] = self.depth[:, self.normal_col]
        frame = frame[np.isfinite(frame["tumor"])]
        agg = {"tumor": ["sum", "count"]}
        if "normal" in frame:
            agg["normal"] = "sum"
        blocks = frame.groupby(["chrom", "block"]).agg(agg)
        blocks.columns = ["tumor_sum", "n_bins"] + (["normal_sum"] if "normal" in frame else [])
        blocks = blocks.reset_index()
        blocks["tumor_mean"] = blocks["tumor_sum"] / blocks["n_bins"]
        blocks["start"] = blocks["block"] * block_bp
        blocks["end"] = blocks["start"] + block_bp - 1
        snp = self._snp_blocks(block_bp)
        if snp is not None:
            blocks = blocks.merge(snp, on=["chrom", "block"], how="left")
        else:
            blocks["mbaf"] = np.nan
            blocks["n_snps"] = 0
        blocks["n_snps"] = blocks["n_snps"].fillna(0).astype(int)
        return blocks

    def _snp_blocks(self, block_bp: int) -> Optional[pd.DataFrame]:
        paths = [os.path.join(self.dir, f) for f in ("snp.tsv.gz", "snp.Aallele.npz", "snp.Ballele.npz")]
        if not all(os.path.exists(p) for p in paths):
            return None
        snps = pd.read_csv(paths[0], sep="\t", usecols=["#CHR", "POS"])
        a = np.load(paths[1])["mat"]
        b = np.load(paths[2])["mat"]
        a = a[:, self.tumor_col] if a.ndim > 1 else a
        b = b[:, self.tumor_col] if b.ndim > 1 else b
        frame = pd.DataFrame({
            "chrom": snps["#CHR"].astype(str).values,
            "block": (snps["POS"].values // block_bp).astype(int),
            "minor": np.minimum(a, b).astype(float),
            "total": (a + b).astype(float),
        })
        frame = frame[frame["total"] > 0]
        grouped = frame.groupby(["chrom", "block"]).agg(minor=("minor", "sum"), total=("total", "sum"), n_snps=("minor", "size")).reset_index()
        grouped["mbaf"] = grouped["minor"] / grouped["total"]
        return grouped[["chrom", "block", "mbaf", "n_snps"]]


def _log_mode(values: np.ndarray, bandwidth: float = 0.02) -> float:
    """Mode of positive values on a log scale (smoothed histogram)."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if len(values) == 0:
        raise ValueError("no positive values to take a mode of")
    if len(values) < 3:
        return float(np.median(values))
    x = np.log(values)
    if x.max() - x.min() <= bandwidth:
        # (nearly) identical values: the histogram would be narrower than the kernel
        return float(math.exp(np.median(x)))
    lo, hi = x.min() - 3 * bandwidth, x.max() + 3 * bandwidth
    nbins = max(int(math.ceil((hi - lo) / bandwidth)), 1)
    hist, edges = np.histogram(x, bins=nbins, range=(lo, hi))
    width = 2.0
    radius = int(3 * width)
    offsets = np.arange(-radius, radius + 1)
    kernel = np.exp(-0.5 * (offsets / width) ** 2)
    # 'full' convolution sliced back to the histogram, so index i is always histogram bin i
    # (mode='same' returns len(kernel) values when the kernel is longer than the histogram)
    smooth = np.convolve(hist.astype(float), kernel, mode="full")[radius:radius + len(hist)]
    i = int(np.argmax(smooth))
    return float(math.exp(0.5 * (edges[i] + edges[i + 1])))


def _outside(chrom: str, start: int, end: int, regions: List[Tuple[str, int, int]], pad: int) -> bool:
    for c, s, e in regions:
        if c == chrom and start <= e + pad and end >= s - pad:
            return False
    return True


# ---------------------------------------------------------------------------------------- #
# gamma / theta
# ---------------------------------------------------------------------------------------- #

def estimate_depth_scale(windows: Optional[Windows], regions: List[Tuple[str, int, int]], params: ScaleParams) -> Tuple[float, str, Optional[float], Dict[str, object]]:
    """Return (gamma, method, ratio_baseline, details)."""
    details: Dict[str, object] = {}
    if params.depth_per_copy is not None:
        return float(params.depth_per_copy), "override", None, details
    if windows is None:
        raise ValueError("no windows directory and no --depth-per-copy: cannot convert depth to copy number")
    if params.hatchet_seg:
        gamma, info = _gamma_from_hatchet(windows, regions, params)
        details.update(info)
        if params.depth_mode == "ratio":
            log.warning("depth_mode='ratio' is ignored with --hatchet-seg; segment copy numbers use tumor depth / gamma")
            params.depth_mode = "tumor"
        return gamma, "hatchet", None, details

    blocks = windows.block_table(params.block_bp)
    min_bins = 0.5 * params.block_bp / windows.bin_size
    keep = blocks["n_bins"] >= min_bins
    keep &= np.array([_outside(c, s, e, regions, params.exclude_pad_bp) for c, s, e in zip(blocks["chrom"], blocks["start"], blocks["end"])])
    balanced = keep & (blocks["n_snps"] >= params.min_block_snps) & (blocks["mbaf"] >= params.balanced_min_mbaf)
    used = balanced
    if balanced.sum() < 10:
        log.warning("only %d allelically balanced blocks; using all %d background blocks for the depth mode", int(balanced.sum()), int(keep.sum()))
        used = keep
    details.update({"blocks_total": int(len(blocks)), "blocks_background": int(keep.sum()), "blocks_balanced": int(balanced.sum()),
                    "blocks_used": int(used.sum()), "block_bp": params.block_bp})
    tumor_mode = _log_mode(blocks.loc[used, "tumor_mean"].values)
    gamma = tumor_mode / params.baseline_cn
    details["tumor_depth_mode"] = tumor_mode
    ratio_baseline = None
    if params.depth_mode == "ratio":
        if "normal_sum" not in blocks:
            raise ValueError("depth_mode='ratio' needs a normal sample")
        ratios = (blocks.loc[used, "tumor_sum"] / blocks.loc[used, "normal_sum"]).values
        ratio_baseline = _log_mode(ratios)
        details["ratio_mode"] = ratio_baseline
    return gamma, "tumor_mode" if params.depth_mode == "tumor" else "ratio_mode", ratio_baseline, details


def _parse_hatchet(path: str, sample: Optional[str]) -> pd.DataFrame:
    seg = pd.read_csv(path, sep="\t")
    if "SAMPLE" in seg and sample is not None and (seg["SAMPLE"].astype(str) == sample).any():
        seg = seg[seg["SAMPLE"].astype(str) == sample]
    elif "SAMPLE" in seg:
        seg = seg[seg["SAMPLE"] == seg["SAMPLE"].iloc[0]]
    clones = [c[3:] for c in seg.columns if c.startswith("cn_")]
    frac = np.zeros(len(seg))
    for clone in clones:
        cn = seg["cn_" + clone].astype(str).str.split("|")
        total = cn.apply(lambda parts: sum(int(p) for p in parts))
        frac += total.values * seg["u_" + clone].astype(float).values
    seg = seg.assign(frac_cn=frac)
    return seg


def _gamma_from_hatchet(windows: Windows, regions, params: ScaleParams) -> Tuple[float, Dict[str, object]]:
    seg = _parse_hatchet(params.hatchet_seg, params.sample)
    ratios, lengths = [], []
    for row in seg.to_dict("records"):
        chrom, start, end = str(row["#CHR"]), int(row["START"]), int(row["END"])
        if end - start < params.hatchet_min_segment_bp or row["frac_cn"] <= 0:
            continue
        if not _outside(chrom, start, end, regions, params.exclude_pad_bp):
            continue
        depth, nb = windows.mean_depth(chrom, start, end - 1)
        if nb <= 0 or not math.isfinite(depth):
            continue
        ratios.append(depth / row["frac_cn"])
        lengths.append(end - start)
    if not ratios:
        raise ValueError("no usable HATCHet segments to fit depth per copy")
    order = np.argsort(ratios)
    cum = np.cumsum(np.array(lengths)[order])
    gamma = float(np.array(ratios)[order][np.searchsorted(cum, cum[-1] / 2.0)])  # length-weighted median
    return gamma, {"hatchet_segments_used": len(ratios), "hatchet_file": params.hatchet_seg}


def estimate_junction_scale(graph: ECGraph, windows: Optional[Windows], gamma: float, regions, params: ScaleParams,
                            extremities_path: Optional[str] = None) -> Tuple[float, str, Dict[str, object]]:
    """Return (theta, method, details): reads spanning one copy of a junction."""
    if params.junction_reads_per_copy is not None:
        return float(params.junction_reads_per_copy), "override", {}
    if extremities_path and os.path.exists(extremities_path) and windows is not None:
        ext = pd.read_csv(extremities_path, sep="\t")
        pos_col = "POS0" if "POS0" in ext else "START"
        ref_sum, cn_sum, n = 0.0, 0.0, 0
        for row in ext.to_dict("records"):
            chrom, pos = str(row["#CHR"]), int(row[pos_col])
            if not _outside(chrom, pos, pos, regions, params.exclude_pad_bp):
                continue
            depth, nb = windows.mean_depth(chrom, pos - params.theta_window_bp, pos + params.theta_window_bp)
            if nb <= 0 or not math.isfinite(depth):
                continue
            cn = depth / gamma
            if not 0.5 <= cn <= 8.0:
                continue
            ref_sum += float(row["#REF"])
            cn_sum += cn
            n += 1
        if n >= 5 and cn_sum > 0:
            return ref_sum / cn_sum, "extremities", {"extremities_used": n}
        log.warning("too few background extremities (%d) to estimate junction reads per copy", n)
    reads = [a.reads for a in graph.adjacencies.values() if a.kind == REFERENCE and a.reads is not None and a.reads > 0]
    if not reads:
        raise ValueError("cannot estimate junction reads per copy; pass --junction-reads-per-copy")
    theta = float(np.median(reads)) / params.baseline_cn
    log.warning("junction reads per copy estimated from the subgraph's own reference adjacencies (%.2f); "
                "pass --extremities or --junction-reads-per-copy for a better estimate", theta)
    return theta, "subgraph_reference_median", {"reference_adjacencies_used": len(reads)}


# ---------------------------------------------------------------------------------------- #
# observations
# ---------------------------------------------------------------------------------------- #

def build_observations(graph: ECGraph, params: ScaleParams, windows_dir: Optional[str] = None,
                       extremities_path: Optional[str] = None) -> Observations:
    """Observed copy numbers for every segment, adjacency and background flank of ``graph``."""
    windows = Windows(windows_dir, sample=params.sample) if windows_dir else None
    regions = graph.regions()
    gamma, depth_method, ratio_baseline, details = estimate_depth_scale(windows, regions, params)
    theta, junction_method, jdetails = estimate_junction_scale(graph, windows, gamma, regions, params, extremities_path)
    details.update(jdetails)
    scale = Scale(depth_per_copy=gamma, junction_reads_per_copy=theta, baseline_cn=params.baseline_cn,
                  depth_method=depth_method, junction_method=junction_method, depth_mode=params.depth_mode,
                  ratio_baseline=ratio_baseline, details=details)
    log.info("depth per copy %.3f (%s), junction reads per copy %.3f (%s)", gamma, depth_method, theta, junction_method)

    def cn_over(chrom, start, end):
        if windows is None:
            return float("nan"), float("nan"), 0.0
        if params.depth_mode == "ratio" and ratio_baseline:
            r, nb = windows.ratio(chrom, start, end)
            depth, _ = windows.mean_depth(chrom, start, end)
            return params.baseline_cn * r / ratio_baseline, depth, nb
        depth, nb = windows.mean_depth(chrom, start, end)
        return depth / gamma, depth, nb

    segments: Dict[str, SegmentObservation] = {}
    for s in graph.segments.values():
        total, depth, nb = cn_over(s.chrom, s.start, s.end)
        source = "windows"
        if not (nb > 0 and math.isfinite(total)):
            depth = s.depth
            total = s.depth / gamma
            nb = s.length / 1000.0
            source = "subgraph_weight"
        a = b = None
        if s.baf is not None and s.n_snps >= params.min_snps_allele:
            baf = min(max(float(s.baf), 0.0), 1.0)
            a, b = total * (1.0 - baf), total * baf
        segments[s.sid] = SegmentObservation(sid=s.sid, depth=float(depth), total=float(total), a=a, b=b, n_bins=float(nb), source=source)

    adjacencies: Dict[str, Optional[float]] = {}
    for a in graph.adjacencies.values():
        adjacencies[a.aid] = None if a.reads is None else float(a.reads) / theta

    background: Dict[str, float] = {}
    flank_cn: Dict[str, float] = {}
    for chrom, start, end in regions:
        values, weights = [], []
        for lo, hi in ((start - params.flank_bp, start - params.flank_gap_bp), (end + params.flank_gap_bp, end + params.flank_bp)):
            if hi <= lo or hi < 0:
                continue
            cn, _, nb = cn_over(chrom, max(lo, 0), hi)
            if nb > 0 and math.isfinite(cn):
                values.append(cn)
                weights.append(nb)
        flank_cn["{c}:{s}-{e}".format(c=chrom, s=start, e=end)] = float(np.average(values, weights=weights)) if values else params.baseline_cn
    for s in graph.segments.values():
        value = params.baseline_cn
        for chrom, start, end in regions:
            if chrom == s.chrom and start <= s.start and s.end <= end:
                value = flank_cn["{c}:{s}-{e}".format(c=chrom, s=start, e=end)]
                break
        background[s.sid] = value
    scale.details["background_flank_cn"] = flank_cn
    return Observations(scale=scale, segments=segments, adjacencies=adjacencies, background=background)
