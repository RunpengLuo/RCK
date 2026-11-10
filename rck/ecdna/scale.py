"""Convert read depth and junction read counts into copy numbers.

Two scales are needed:

* ``depth_per_copy`` (gamma): tumor read depth contributed by one copy of a locus, averaged
  over the bulk sample. Segment copy number is ``mean tumor depth / gamma``.
* ``junction_reads_per_copy`` (theta): reads spanning one copy of a junction. Junction copy
  number is ``#NOV / theta`` (novel) or ``#REF / theta`` (reference join).

Gamma is taken from the first available source:

0. the subgraph folder itself: ``scale.json`` written by amplicon-allele-analysis
   ``subgraph_construct.py`` (``rck_inputs.py``), which fits gamma, purity and tumor ploidy from
   genome-wide depth and SNPs and also stores theta, per-segment copy numbers (``CN_TOTAL``,
   ``CN_A``, ``CN_B`` edge columns) and a copy-number profile for the region flanks
   (``cn_profile.tsv``). With these, RCK needs no windows, extremities or HATCHet input
   (see :func:`load_subgraph_scale`; ``--depth-per-copy`` still overrides);
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

from rck.ecdna.subgraph_io import ECGraph, REFERENCE, SegmentBlocks

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
    block_flip: str = "off"  # off | unoriented: allele observations per phase-block piece (segment_blocks.tsv)
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
class PieceObservation:
    """Observation of the part of a segment inside one phase block (``--block-flip unoriented``).

    ``a`` / ``b`` (and ``count_a`` / ``count_b``) are in the piece's base frame: oriented by the
    subgraph when the block is ``fixed`` (orientation 1 or -1), observed when it is not
    (orientation 0 or missing). The flow may swap them in a block that is not fixed. ``a`` and
    ``b`` are None for a total-only piece (fewer than ``--min-snps-allele`` SNPs). ``share`` is
    the piece's part of its segment's data weight (by length)."""

    sid: str
    block: str
    start: int
    end: int
    n_snps: int
    total: float
    a: Optional[float]
    b: Optional[float]
    count_a: float
    count_b: float
    orientation: int
    share: float

    @property
    def fixed(self) -> bool:
        return self.orientation != 0

    @property
    def has_alleles(self) -> bool:
        return self.a is not None and self.b is not None

    def alleles(self, flip: int = 0) -> Tuple[Optional[float], Optional[float]]:
        return (self.b, self.a) if flip else (self.a, self.b)

    def counts(self, flip: int = 0) -> Tuple[float, float]:
        return (self.count_b, self.count_a) if flip else (self.count_a, self.count_b)


@dataclass
class Observations:
    scale: Scale
    segments: Dict[str, SegmentObservation]
    adjacencies: Dict[str, Optional[float]]
    background: Dict[str, float]
    pieces: Dict[str, List[PieceObservation]] = field(default_factory=dict)  # sid -> pieces, with --block-flip only
    block_flip: str = "off"
    min_snps_allele: int = 5

    def has_alleles(self, sid: str) -> bool:
        """Whether the segment has any allele-specific observation (in a piece, or the segment's own)."""
        if self.pieces.get(sid):
            return any(p.has_alleles for p in self.pieces[sid])
        so = self.segments[sid]
        return so.a is not None and so.b is not None

    def block_orientation(self) -> Dict[str, int]:
        """Subgraph orientation of every phase block that has a piece."""
        return {p.block: p.orientation for ps in self.pieces.values() for p in ps}

    def segment_alleles(self, sid: str, flips: Optional[Dict[str, int]] = None) -> Tuple[Optional[float], Optional[float]]:
        """Segment-level A / B copy number. With pieces: the segment total split by the BAF of the
        piece counts summed after ``flips`` (as subgraph_construct.py pools a segment's blocks)."""
        so = self.segments[sid]
        pieces = self.pieces.get(sid)
        if not pieces:
            return so.a, so.b
        flips = flips or {}
        a = b = 0.0
        for p in pieces:
            ca, cb = p.counts(flips.get(p.block, 0))
            a, b = a + ca, b + cb
        if sum(p.n_snps for p in pieces) < self.min_snps_allele or a + b <= 0:
            return None, None
        return so.total * a / (a + b), so.total * b / (a + b)

    def write_json(self, path: str, flips: Optional[Dict[str, int]] = None) -> None:
        payload = {
            "scale": self.scale.to_dict(),
            "segments": {k: asdict(v) for k, v in self.segments.items()},
            "adjacencies": self.adjacencies,
            "background_total_cn": self.background,
        }
        if self.block_flip != "off":
            flips = flips or {}
            payload["block_flip"] = self.block_flip
            payload["pieces"] = []
            for sid in sorted(self.pieces):
                for p in self.pieces[sid]:
                    flip = int(flips.get(p.block, 0))
                    a, b = p.alleles(flip)
                    payload["pieces"].append(dict(asdict(p), fixed=p.fixed, flipped=bool(flip), a_final=a, b_final=b))
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
# scale inputs stored in the subgraph folder
# ---------------------------------------------------------------------------------------- #

SUBGRAPH_SCALE_FORMAT = "rck-ecdna-subgraph-scale/1"


@dataclass
class SubgraphScale:
    """``scale.json`` (+ ``cn_profile.tsv``) of a subgraph folder written by subgraph_construct.py."""

    path: str
    data: Dict[str, object]
    profile: Optional["CNProfile"] = None

    def get(self, key, default=None):
        value = self.data.get(key, default)
        return default if value is None else value


def load_subgraph_scale(subgraph_dir: str) -> Optional[SubgraphScale]:
    path = os.path.join(subgraph_dir, "scale.json")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        data = json.load(handle)
    if data.get("format") != SUBGRAPH_SCALE_FORMAT:
        log.warning("%s has format %r, expected %r; reading it anyway", path, data.get("format"), SUBGRAPH_SCALE_FORMAT)
    profile = None
    name = (data.get("files") or {}).get("cn_profile") or "cn_profile.tsv"
    if os.path.exists(os.path.join(subgraph_dir, name)):
        profile = CNProfile(os.path.join(subgraph_dir, name))
    return SubgraphScale(path=path, data=data, profile=profile)


class CNProfile:
    """Binned cell-averaged copy number (``#CHR START END CN_TOTAL N_BINS``), for region flanks."""

    def __init__(self, path: str):
        frame = pd.read_csv(path, sep="\t")
        self.by_chrom = {str(c): (g["START"].values.astype(np.int64), g["END"].values.astype(np.int64), g["CN_TOTAL"].values.astype(float),
                                  g["N_BINS"].values.astype(float) if "N_BINS" in g else np.ones(len(g)))
                         for c, g in frame.groupby("#CHR", sort=False)}

    def mean_cn(self, chrom: str, start: int, end_inclusive: int) -> Tuple[float, float]:
        """(overlap-weighted mean copy number, overlapped bins in kb) over [start, end_inclusive]."""
        if chrom not in self.by_chrom:
            return float("nan"), 0.0
        s, e, cn, nb = self.by_chrom[chrom]
        ov = np.clip(np.minimum(e, end_inclusive + 1) - np.maximum(s, start), 0, None).astype(float)
        keep = (ov > 0) & np.isfinite(cn) & (nb > 0)
        if not keep.any():
            return float("nan"), 0.0
        return float(np.average(cn[keep], weights=ov[keep])), float(ov[keep].sum() / 1000.0)


# ---------------------------------------------------------------------------------------- #
# observations
# ---------------------------------------------------------------------------------------- #

def build_observations(graph: ECGraph, params: ScaleParams, windows_dir: Optional[str] = None,
                       extremities_path: Optional[str] = None, subgraph_scale: Optional[SubgraphScale] = None,
                       segment_blocks: Optional[SegmentBlocks] = None) -> Observations:
    """Observed copy numbers for every segment, adjacency and background flank of ``graph``.

    With ``subgraph_scale`` (subgraph-only mode) windows are not read: gamma and theta come from
    ``scale.json`` (overrides still win), segment copy numbers from the ``CN_*`` edge columns
    and the flanks from ``cn_profile.tsv``.

    With ``params.block_flip`` other than ``off`` and ``segment_blocks``, every segment that has
    pieces in ``segment_blocks.tsv`` also gets one observation per piece (see
    :class:`PieceObservation` and :func:`piece_observations`).
    """
    windows = Windows(windows_dir, sample=params.sample) if (windows_dir and subgraph_scale is None) else None
    regions = graph.regions()
    details: Dict[str, object] = {}
    ratio_baseline = None
    if subgraph_scale is not None:
        for key in ("sample", "purity", "tumor_ploidy", "cell_ploidy", "ploidy_class", "method", "junction_sample"):
            details["subgraph_" + key] = subgraph_scale.data.get(key)
        details["subgraph_alternatives"] = subgraph_scale.data.get("alternatives")
        if params.depth_per_copy is not None:
            gamma, depth_method = float(params.depth_per_copy), "override"
        else:
            gamma, depth_method = float(subgraph_scale.get("depth_per_copy")), "subgraph_scale"
        if params.junction_reads_per_copy is not None:
            theta, junction_method = float(params.junction_reads_per_copy), "override"
        elif subgraph_scale.get("junction_reads_per_copy"):
            theta, junction_method = float(subgraph_scale.get("junction_reads_per_copy")), "subgraph_scale"
        else:
            theta, junction_method, jdetails = estimate_junction_scale(graph, None, gamma, regions, params, None)
            details.update(jdetails)
    else:
        gamma, depth_method, ratio_baseline, details = estimate_depth_scale(windows, regions, params)
        theta, junction_method, jdetails = estimate_junction_scale(graph, windows, gamma, regions, params, extremities_path)
        details.update(jdetails)
    profile = subgraph_scale.profile if subgraph_scale is not None else None
    scale = Scale(depth_per_copy=gamma, junction_reads_per_copy=theta, baseline_cn=params.baseline_cn,
                  depth_method=depth_method, junction_method=junction_method, depth_mode=params.depth_mode,
                  ratio_baseline=ratio_baseline, details=details)
    log.info("depth per copy %.3f (%s), junction reads per copy %.3f (%s)", gamma, depth_method, theta, junction_method)

    def cn_over(chrom, start, end):
        if windows is None:
            if profile is not None:
                cn, kb = profile.mean_cn(chrom, start, end)
                return cn, cn * gamma, kb
            return float("nan"), float("nan"), 0.0
        if params.depth_mode == "ratio" and ratio_baseline:
            r, nb = windows.ratio(chrom, start, end)
            depth, _ = windows.mean_depth(chrom, start, end)
            return params.baseline_cn * r / ratio_baseline, depth, nb
        depth, nb = windows.mean_depth(chrom, start, end)
        return depth / gamma, depth, nb

    segments: Dict[str, SegmentObservation] = {}
    for s in graph.segments.values():
        if subgraph_scale is not None and s.cn_total is not None and params.depth_per_copy is None:
            total, nb, source = float(s.cn_total), s.length / 1000.0, "subgraph_cn"
            depth = total * gamma
            a = b = None
            if s.cn_a is not None and s.cn_b is not None and s.n_snps >= params.min_snps_allele:
                a, b = float(s.cn_a), float(s.cn_b)
            elif s.baf is not None and s.n_snps >= params.min_snps_allele:
                baf = min(max(float(s.baf), 0.0), 1.0)
                a, b = total * (1.0 - baf), total * baf
            segments[s.sid] = SegmentObservation(sid=s.sid, depth=float(depth), total=total, a=a, b=b, n_bins=float(nb), source=source)
            continue
        total, depth, nb = cn_over(s.chrom, s.start, s.end) if windows is not None else (float("nan"), float("nan"), 0.0)
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
    obs = Observations(scale=scale, segments=segments, adjacencies=adjacencies, background=background,
                       block_flip=params.block_flip, min_snps_allele=params.min_snps_allele)
    if params.block_flip != "off" and segment_blocks is not None:
        sub_gamma = float(subgraph_scale.get("depth_per_copy")) if subgraph_scale is not None and subgraph_scale.get("depth_per_copy") else None
        obs.pieces = piece_observations(graph, segment_blocks, segments, params, gamma, sub_gamma,
                                        cn_over if windows is not None else None)
        log.info("block flips: %d pieces on %d segments, %d phase blocks (%d not oriented by the subgraph)",
                 sum(len(v) for v in obs.pieces.values()), len(obs.pieces), len(obs.block_orientation()),
                 sum(1 for o in obs.block_orientation().values() if o == 0))
    return obs


def piece_observations(graph: ECGraph, blocks: SegmentBlocks, segments: Dict[str, SegmentObservation], params: ScaleParams,
                       gamma: float, sub_gamma: Optional[float], cn_over=None) -> Dict[str, List[PieceObservation]]:
    """Per-piece observations: total, A / B in the base frame, and the share of the segment weight.

    The total is the piece's ``CN_TOTAL`` (rescaled when ``--depth-per-copy`` overrides the
    subgraph's depth per copy), or the windows' copy number over the piece in windows mode. A / B
    come from the ``CN_A`` / ``CN_B`` columns in subgraph mode, else from the piece total split by
    its BAF; pieces with fewer than ``--min-snps-allele`` SNPs are total-only. Pieces of a block
    with orientation -1 are swapped to the oriented frame and marked fixed."""
    by_sid: Dict[str, List] = {}
    for bp in blocks.pieces:
        by_sid.setdefault(bp.sid, []).append(bp)
    result: Dict[str, List[PieceObservation]] = {}
    for sid, rows in by_sid.items():
        s, so = graph.segments[sid], segments[sid]
        rows = sorted(rows, key=lambda r: r.start)
        span = float(sum(r.length for r in rows))
        out = []
        for r in rows:
            total, a, b = None, None, None
            if cn_over is not None:
                cn, _, nb = cn_over(s.chrom, r.start, r.end - 1)
                total = cn if nb > 0 and math.isfinite(cn) else None
            elif r.cn_total is not None:
                total = r.cn_total if params.depth_per_copy is None or not sub_gamma else r.cn_total * sub_gamma / gamma
                if params.depth_per_copy is None and r.cn_a is not None and r.cn_b is not None:
                    a, b = r.cn_a, r.cn_b
            if total is None:
                total = so.total
            if r.n_snps < params.min_snps_allele or r.count_a + r.count_b <= 0:
                a = b = None
            elif a is None or b is None:
                baf = r.count_b / (r.count_a + r.count_b)
                a, b = total * (1.0 - baf), total * baf
            orientation = int(blocks.orientation.get(r.block, 0))
            ca, cb = r.count_a, r.count_b
            if orientation == -1:
                a, b, ca, cb = b, a, cb, ca
            out.append(PieceObservation(sid=sid, block=r.block, start=r.start, end=r.end, n_snps=r.n_snps, total=float(total),
                                        a=None if a is None else float(a), b=None if b is None else float(b),
                                        count_a=float(ca), count_b=float(cb), orientation=orientation,
                                        share=r.length / span if span > 0 else 1.0 / len(rows)))
        result[sid] = out
    return result
