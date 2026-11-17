"""RCK's simulated-input degradations (Aganezov & Raphael 2020) applied to a multi-clone DIAG folder.

The multi-clone route reads ``<subgraph>/rck_input/`` (RCK's native files; ``rck.ecdna.multiclone``). Built
with ``subgraph_construct.py --clone-profiles truth`` it holds the true clone- and allele-specific copy
numbers, which is what RCK's simulation study gave RCK before degrading it ("Evaluation of RCK on simulated
data"; Methods "Simulating rearranged cancer samples"). This module writes the degraded copies
(results/bench_IW_allele/PLAN.md, Part B):

    chil      clone and haplotype information loss: fragments = segments; each fragment's A / B swapped with
              probability 0.5 (the same swap in every clone); junction haplotype tags and labelling groups
              removed (molecule groups kept: RCK's Methods "Third-generation sequencing technologies and novel
              adjacency groups")
    chil_avg  copy numbers averaged over 50 kb fragments (length-weighted, per clone and allele; see
              ``window_fragments`` for where windows are cut), then CHIL
    chil_err  ±50 bp on both breakends of a random half of the novel adjacencies, 10% spurious novel
              adjacencies, 50 kb averaging, ±1 on a random 5% of the fragment copy-number values (per clone
              and allele), then CHIL

A moved or new breakend is aligned to the segments by splitting the segment it falls in, copying its copy
numbers and bounds (RCK aligns segment extremities to novel adjacency positions the same way: Methods
"Deriving extremities and novel adjacencies from data"). The multi-clone ILP takes junction haplotype
information only through labelling groups (``agt=L`` in rck.adjgroups.tsv), so removing them is the label
loss; the ``junction_hap`` / ``hap_obs_end*`` tags are dropped too so no file still carries them.

Coordinates: segments are 1-based inclusive; a "-" extremity sits at a segment start, a "+" at a segment end.
"""
from __future__ import annotations

import ast
import json
import os
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

AVG_BP = 50_000
SHIFT_BP = 50
SHIFT_SHARE = 0.5
SPURIOUS_SHARE = 0.10
PERTURB_SHARE = 0.05
FLIP_P = 0.5
LABEL_TAGS = ("junction_hap", "hap_obs_end1", "hap_obs_end2")
VARIANTS = ("chil", "chil_avg", "chil_err")


# ---------------------------------------------------------------------------------------------- #
# files
# ---------------------------------------------------------------------------------------------- #

@dataclass
class Segment:
    chrom: str
    start: int
    end: int
    extra: Dict[str, object]  # parsed key=value pairs; cn / cnb as dicts

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_extra(text: str) -> Dict[str, object]:
    out: Dict[str, object] = {}
    for part in (text or "").split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        out[key] = ast.literal_eval(value) if key in ("cn", "cnb") else value
    return out


def format_extra(extra: Dict[str, object]) -> str:
    return ";".join("{k}={v}".format(k=k, v=v) for k, v in extra.items())


def read_segments(path: str) -> List[Segment]:
    out = []
    with open(path) as handle:
        handle.readline()
        for line in handle:
            f = line.rstrip("\n").split("\t")
            out.append(Segment(f[0], int(f[1]), int(f[2]), parse_extra(f[3] if len(f) > 3 else "")))
    return out


def write_segments(path: str, segments: Sequence[Segment]) -> None:
    with open(path, "w") as handle:
        handle.write("chr\tstart\tend\textra\n")
        for s in segments:
            handle.write("{c}\t{a}\t{b}\t{e}\n".format(c=s.chrom, a=s.start, b=s.end, e=format_extra(s.extra)))


@dataclass
class Adjacency:
    aid: str
    chrom1: str
    coord1: int
    strand1: str
    chrom2: str
    coord2: int
    strand2: str
    extra: Dict[str, str] = field(default_factory=dict)


def read_adjacencies(path: str) -> List[Adjacency]:
    out = []
    with open(path) as handle:
        handle.readline()
        for line in handle:
            f = line.rstrip("\n").split("\t")
            out.append(Adjacency(f[0], f[1], int(f[2]), f[3], f[4], int(f[5]), f[6], parse_extra(f[7] if len(f) > 7 else "")))
    return out


def write_adjacencies(path: str, adjacencies: Sequence[Adjacency]) -> None:
    with open(path, "w") as handle:
        handle.write("aid\tchr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\textra\n")
        for a in adjacencies:
            handle.write("\t".join(str(x) for x in (a.aid, a.chrom1, a.coord1, a.strand1, a.chrom2, a.coord2, a.strand2,
                                                    format_extra(a.extra))) + "\n")


# ---------------------------------------------------------------------------------------------- #
# segment surgery
# ---------------------------------------------------------------------------------------------- #

def _copy(s: Segment, start: int, end: int) -> Segment:
    return Segment(s.chrom, start, end, {k: (ast.literal_eval(repr(v)) if isinstance(v, dict) else v) for k, v in s.extra.items()})


def split_at(segments: List[Segment], chrom: str, coord: int, strand: str) -> List[Segment]:
    """Segments with an extremity at (chrom, coord, strand): a "-" end needs a segment starting at coord, a "+"
    end one ending there. The segment holding coord is split in two, both halves keeping its copy numbers."""
    out = []
    for s in segments:
        if s.chrom == chrom and s.start <= coord <= s.end:
            if strand == "-" and coord > s.start:
                out += [_copy(s, s.start, coord - 1), _copy(s, coord, s.end)]
                continue
            if strand == "+" and coord < s.end:
                out += [_copy(s, s.start, coord), _copy(s, coord + 1, s.end)]
                continue
        out.append(s)
    return out


def covered(segments: Sequence[Segment], chrom: str, coord: int) -> bool:
    return any(s.chrom == chrom and s.start <= coord <= s.end for s in segments)


def window_fragments(segments: List[Segment], width: int = AVG_BP) -> Tuple[List[Segment], List[List[int]]]:
    """RCK's 50 kb fragments, where they matter. Windows are [k*width + 1, (k+1)*width]. A window that
    overlaps two or more segments is one fragment: the segments are split at its edges and averaged over it.
    A stretch of one segment covering whole windows alone stays one segment and one fragment (averaging over
    it changes nothing; cutting it into 50 kb fragments would only add fragment flips the model can undo for
    free, and on these DIAGs 99.7% of the bases are in segments > 50 kb, mostly 1 Mb).

    Returns (segments, fragments as lists of segment indices)."""
    def win(x):
        return (x - 1) // width

    shared = set()  # (chrom, window) overlapped by >= 2 segments
    count: Dict[Tuple[str, int], int] = {}
    for s in segments:
        for k in range(win(s.start), win(s.end) + 1):
            count[(s.chrom, k)] = count.get((s.chrom, k), 0) + 1
    shared = {key for key, n in count.items() if n > 1}
    for chrom, k in sorted(shared):
        for edge in (k * width, (k + 1) * width):  # "+" ends at the window's edges
            if edge > 0:
                segments = split_at(segments, chrom, edge, "+")
    segments = sorted(segments, key=lambda x: (x.chrom, x.start))
    fragments: List[List[int]] = []
    for i, s in enumerate(segments):
        key = (s.chrom, win(s.start))
        if fragments and key in shared and win(s.end) == key[1]:
            prev = segments[fragments[-1][-1]]
            if prev.chrom == s.chrom and prev.end + 1 == s.start and (prev.chrom, win(prev.start)) == key and win(prev.end) == key[1]:
                fragments[-1].append(i)
                continue
        fragments.append([i])
    return segments, fragments


def segment_fragments(segments: Sequence[Segment]) -> List[List[int]]:
    return [[i] for i in range(len(segments))]


# ---------------------------------------------------------------------------------------------- #
# copy-number transforms
# ---------------------------------------------------------------------------------------------- #

def _num(x: float):
    r = round(float(x), 3)
    return int(r) if r == int(r) else r


def average(segments: List[Segment], fragments: Sequence[Sequence[int]]) -> None:
    """Each fragment's copy numbers (per clone and allele) replaced by their length-weighted mean."""
    for frag in fragments:
        total = sum(segments[i].length for i in frag)
        clones = segments[frag[0]].extra["cn"]
        mean = {c: {h: _num(sum(segments[i].extra["cn"][c][h] * segments[i].length for i in frag) / total) for h in hv}
                for c, hv in clones.items()}
        for i in frag:
            segments[i].extra["cn"] = {c: dict(hv) for c, hv in mean.items()}


def perturb(segments: List[Segment], fragments: Sequence[Sequence[int]], rng: np.random.Generator,
            share: float = PERTURB_SHARE) -> int:
    """±1 on a random ``share`` of the fragment copy-number values (fragment x clone x allele); a value
    below 1 only goes up. Returns the number changed."""
    entries = [(k, c, h) for k, frag in enumerate(fragments) for c, hv in segments[frag[0]].extra["cn"].items() for h in hv]
    n = int(round(share * len(entries)))
    for j in rng.choice(len(entries), size=n, replace=False):
        k, c, h = entries[int(j)]
        value = segments[fragments[k][0]].extra["cn"][c][h]
        step = 1 if value < 1 or rng.random() < 0.5 else -1
        for i in fragments[k]:
            segments[i].extra["cn"][c][h] = _num(value + step)
    return n


def chil_flip(segments: List[Segment], fragments: Sequence[Sequence[int]], rng: np.random.Generator,
              p: float = FLIP_P) -> List[int]:
    """Swap A and B (copy numbers and bounds, every clone) of each fragment with probability p; returns the
    flipped fragment indices."""
    flipped = []
    for k, frag in enumerate(fragments):
        if rng.random() < p:
            flipped.append(k)
            for i in frag:
                for key in ("cn", "cnb"):
                    if key in segments[i].extra:
                        segments[i].extra[key] = {c: {"A": hv.get("B"), "B": hv.get("A")} for c, hv in segments[i].extra[key].items()}
    return flipped


# ---------------------------------------------------------------------------------------------- #
# adjacency transforms
# ---------------------------------------------------------------------------------------------- #

def _ends(adjacencies: Sequence[Adjacency]) -> set:
    return {(a.chrom1, a.coord1, a.strand1) for a in adjacencies} | {(a.chrom2, a.coord2, a.strand2) for a in adjacencies}


def _clear(end, taken, gap: int) -> bool:
    return not any(c == end[0] and abs(x - end[1]) <= gap for c, x, _ in taken)


def shift_breakends(adjacencies: List[Adjacency], segments: Sequence[Segment], rng: np.random.Generator,
                    share: float = SHIFT_SHARE, bp: int = SHIFT_BP) -> List[str]:
    """Both breakends of a random ``share`` of the adjacencies moved by an integer in [-bp, bp] \\ {0}, kept
    inside the segments and off every other breakend; returns the moved aids."""
    n = int(round(share * len(adjacencies)))
    chosen = sorted(int(j) for j in rng.choice(len(adjacencies), size=n, replace=False))
    moved = []
    for j in chosen:
        a = adjacencies[j]
        others = _ends(adjacencies) - {(a.chrom1, a.coord1, a.strand1), (a.chrom2, a.coord2, a.strand2)}
        new = []
        for chrom, coord in ((a.chrom1, a.coord1), (a.chrom2, a.coord2)):
            for _ in range(100):
                d = int(rng.integers(1, bp + 1)) * (1 if rng.random() < 0.5 else -1)
                if covered(segments, chrom, coord + d) and _clear((chrom, coord + d), others | set(new), 2 * bp):
                    break
            else:
                d = 0
            new.append((chrom, coord + d, None))
        a.coord1, a.coord2 = new[0][1], new[1][1]
        moved.append(a.aid)
    return moved


def spurious_adjacencies(adjacencies: Sequence[Adjacency], segments: Sequence[Segment], rng: np.random.Generator,
                         share: float = SPURIOUS_SHARE, min_gap: int = 1000) -> List[Adjacency]:
    """``max(1, round(share * n))`` new novel adjacencies with both breakends at random positions of the
    amplicon's (non-flank) segments, random strands, at least ``min_gap`` from every other breakend; read
    support (#NOV, #REF) copied from a random true adjacency."""
    n = max(1, int(round(share * len(adjacencies))))
    pool = [s for s in segments if s.extra.get("kind") != "flank" and s.length > 2 * min_gap] or list(segments)
    weights = np.array([s.length for s in pool], dtype=float)
    weights /= weights.sum()
    taken = set(_ends(adjacencies))
    out = []
    for k in range(n):
        ends = []
        while len(ends) < 2:
            s = pool[int(rng.choice(len(pool), p=weights))]
            coord = int(rng.integers(s.start + min_gap, s.end - min_gap + 1))
            end = (s.chrom, coord, "+" if rng.random() < 0.5 else "-")
            if _clear(end, taken, min_gap):
                ends.append(end)
                taken.add(end)
        donor = adjacencies[int(rng.integers(0, len(adjacencies)))]
        extra = {"at": "N", "#NOV": donor.extra.get("#NOV", "20"), "#REF": donor.extra.get("#REF", "0,0"),
                 "imputed": "0", "edge_type": "variant", "spurious": "1"}
        out.append(Adjacency("spurious{k:05d}".format(k=k + 1), ends[0][0], ends[0][1], ends[0][2],
                             ends[1][0], ends[1][1], ends[1][2], extra))
    return out


# ---------------------------------------------------------------------------------------------- #
# one variant
# ---------------------------------------------------------------------------------------------- #

def degrade(src: str, dst: str, variant: str, seed: Sequence[int]) -> Dict[str, object]:
    """Write ``dst`` (an rck_input folder) from ``src`` for ``variant``; returns the record written to
    ``dst/transforms.json``."""
    if variant not in VARIANTS:
        raise ValueError("variant must be one of {v}".format(v=", ".join(VARIANTS)))
    rng = np.random.default_rng(list(seed))
    os.makedirs(dst, exist_ok=True)
    for name in os.listdir(src):
        if name not in ("diag_segments.tsv",) and os.path.isfile(os.path.join(src, name)):
            shutil.copy2(os.path.join(src, name), os.path.join(dst, name))
    segments = read_segments(os.path.join(src, "rck.scnt.tsv"))
    adjacencies = read_adjacencies(os.path.join(src, "rck.adj.tsv"))
    record: Dict[str, object] = {"variant": variant, "seed": list(seed), "source": os.path.abspath(src),
                                 "segments_in": len(segments), "adjacencies_in": len(adjacencies)}

    if variant == "chil_err":
        record["shifted"] = shift_breakends(adjacencies, segments, rng)
        spurious = spurious_adjacencies(adjacencies, segments, rng)
        record["spurious"] = [a.aid for a in spurious]
        adjacencies = adjacencies + spurious
    for a in adjacencies:  # align every breakend (moved and new) to a segment extremity
        segments = split_at(segments, a.chrom1, a.coord1, a.strand1)
        segments = split_at(segments, a.chrom2, a.coord2, a.strand2)
    segments = sorted(segments, key=lambda s: (s.chrom, s.start))

    if variant == "chil":
        fragments = segment_fragments(segments)
    else:
        segments, fragments = window_fragments(segments)
        average(segments, fragments)
        if variant == "chil_err":
            record["perturbed"] = perturb(segments, fragments, rng)
    record["flipped_fragments"] = chil_flip(segments, fragments, rng)
    record.update(segments_out=len(segments), fragments=len(fragments), adjacencies_out=len(adjacencies))

    for a in adjacencies:
        for tag in LABEL_TAGS:
            a.extra.pop(tag, None)
    write_segments(os.path.join(dst, "rck.scnt.tsv"), segments)
    write_segments(os.path.join(dst, "rck.scnb.tsv"),
                   [Segment(s.chrom, s.start, s.end, {"cnb": s.extra["cnb"]}) for s in segments if "cnb" in s.extra])
    with open(os.path.join(dst, "rck.fragments.tsv"), "w") as handle:
        handle.write("chr\tstart\tend\textra\n")
        for k, frag in enumerate(fragments):
            first, last = segments[frag[0]], segments[frag[-1]]
            handle.write("{c}\t{a}\t{b}\tfragment=f{k:05d};flipped={f}\n".format(
                c=first.chrom, a=first.start, b=last.end, k=k, f=int(k in set(record["flipped_fragments"]))))
    write_adjacencies(os.path.join(dst, "rck.adj.tsv"), adjacencies)
    groups = os.path.join(src, "rck.adjgroups.tsv")
    if os.path.exists(groups):
        kept = dropped = 0
        with open(groups) as handle, open(os.path.join(dst, "rck.adjgroups.tsv"), "w") as out:
            out.write(handle.readline())
            for line in handle:
                if "agt=L" in line.split("\t")[-1]:
                    dropped += 1
                    continue
                out.write(line)
                kept += 1
        record.update(groups_kept=kept, labelling_groups_dropped=dropped)
    record["flipped_fragments"] = len(record["flipped_fragments"])
    with open(os.path.join(dst, "transforms.json"), "w") as handle:
        json.dump(record, handle, indent=2)
    return record
