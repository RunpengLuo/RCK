"""Score inferred ecDNA structures against ground truth.

Truth comes either from an explicit table (``--truth-cycles``; columns ``species order chrom
start end strand hap weight``, 0-based inclusive coordinates) or from an
amplicon-allele-analysis simulator run directory (``--truth-dir``). The simulator writes its
walks over species-local segment ids without coordinates or haplotypes
(``truth/walks.walk``), so they are rebuilt here:

1. the species' segments are the rows of ``truth/segments.tsv`` whose copy number exceeds the
   clone-mixture background (``truth/clones.tsv`` x ``truth/proportions.tsv``); contiguous rows
   form the species' intervals;
2. the interval order (which defines the segment numbering) is the one under which every join
   of the walk is either a reference join or a junction of that species
   (``truth/adjacency_clone.tsv`` + ``obs/adjacency.tsv``);
3. haplotypes come from the junction phasings, propagated across reference joins.

Both truth and inferred cycles are compared after merging reference-contiguous steps into
runs, so a different segmentation of the same molecule scores as identical:

* ``jaccard_hap`` / ``jaccard_total`` - base-pair Jaccard of the multiplicity-weighted
  coverage, with and without haplotype;
* junction recall / precision (unphased and phased) within ``tol_bp``;
* ``order_match`` - same cyclic run sequence up to rotation and reversal (coordinates within
  ``tol_bp``), checked against every reported ordering;
* ``cn_rel_err`` - |weight - truth weight| / truth weight;
* ``phasing_accuracy`` - share of the bases both cycles cover that also carry the same
  haplotype (``phasing_accuracy_swap_invariant`` takes the better of A/B as given or swapped);
* ``composite`` - one recovery score in [0, 1], the mean of ``jaccard_total``, ``junction_f1``,
  ``phasing_accuracy``, ``order_match`` and max(0, 1 - ``cn_rel_err``); a missed species scores 0.

``recovery_matrix.tsv`` has one row per truth species: these scores for the rank-1
decomposition, plus the best composite over every reported decomposition.
"""

from __future__ import annotations

import argparse
import itertools
import json
import logging
import os
import re
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

log = logging.getLogger(__name__)


@dataclass(frozen=True)
class CStep:
    chrom: str
    start: int  # inclusive
    end: int  # inclusive
    hap: str
    strand: str


@dataclass
class Cycle:
    name: str
    weight: float
    orderings: List[List[CStep]]
    circular: bool = True  # False for a linear walk (its last step does not join its first)
    weight_int: Optional[float] = None  # inferred species: integer weight before the copy-number refit

    @property
    def steps(self) -> List[CStep]:
        return self.orderings[0]


# ---------------------------------------------------------------------------------------- #
# geometry
# ---------------------------------------------------------------------------------------- #

def _contiguous(a: CStep, b: CStep, tol: int = 2) -> bool:
    if a.chrom != b.chrom or a.hap != b.hap or a.strand != b.strand:
        return False
    if a.strand == "+":
        return abs(b.start - (a.end + 1)) <= tol
    return abs(a.start - (b.end + 1)) <= tol


def runs(steps: Sequence[CStep], tol: int = 2, circular: bool = True) -> List[CStep]:
    """Merge reference-contiguous consecutive steps (cyclically for a circle) into maximal runs."""
    steps = list(steps)
    if not steps:
        return []
    if not circular:
        result = [steps[0]]
        for s in steps[1:]:
            last = result[-1]
            if _contiguous(last, s, tol):
                result[-1] = CStep(last.chrom, min(last.start, s.start), max(last.end, s.end), last.hap, last.strand)
            else:
                result.append(s)
        return result
    n = len(steps)
    start = 0
    for i in range(n):
        if not _contiguous(steps[i - 1], steps[i], tol):
            start = i
            break
    else:
        # every join is contiguous: a circle made of one reference interval
        lo = min(s.start for s in steps)
        hi = max(s.end for s in steps)
        return [CStep(steps[0].chrom, lo, hi, steps[0].hap, steps[0].strand)]
    ordered = steps[start:] + steps[:start]
    result = [ordered[0]]
    for s in ordered[1:]:
        last = result[-1]
        if _contiguous(last, s, tol):
            result[-1] = CStep(last.chrom, min(last.start, s.start), max(last.end, s.end), last.hap, last.strand)
        else:
            result.append(s)
    return result


def _coverage(steps: Sequence[CStep], with_hap: bool) -> Dict[Tuple[str, ...], List[Tuple[int, int]]]:
    events: Dict[Tuple[str, ...], List[Tuple[int, int]]] = {}
    for s in steps:
        key = (s.chrom, s.hap) if with_hap else (s.chrom,)
        events.setdefault(key, []).extend([(s.start, 1), (s.end + 1, -1)])
    return events


def _overlap_bp(a: Sequence[CStep], b: Sequence[CStep], with_hap: bool) -> Tuple[int, int, int]:
    """(overlap, total a, total b) of multiplicity-weighted coverage."""
    ca, cb = _coverage(a, with_hap), _coverage(b, with_hap)
    overlap = 0
    for key in set(ca) & set(cb):
        points = sorted(set(p for p, _ in ca[key]) | set(p for p, _ in cb[key]))
        da, db = {}, {}
        for p, d in ca[key]:
            da[p] = da.get(p, 0) + d
        for p, d in cb[key]:
            db[p] = db.get(p, 0) + d
        xa = xb = 0
        for p, q in zip(points, points[1:]):
            xa += da.get(p, 0)
            xb += db.get(p, 0)
            overlap += min(xa, xb) * (q - p)
    total_a = sum(s.end - s.start + 1 for s in a)
    total_b = sum(s.end - s.start + 1 for s in b)
    return overlap, total_a, total_b


def jaccard(a: Sequence[CStep], b: Sequence[CStep], with_hap: bool = True) -> float:
    ov, ta, tb = _overlap_bp(a, b, with_hap)
    union = ta + tb - ov
    return ov / union if union > 0 else 0.0


Breakend = Tuple[str, int, str]


def junctions(steps: Sequence[CStep], tol: int = 2, circular: bool = True) -> List[Tuple[Breakend, Breakend, str, str]]:
    """Non-contiguous joins of a walk as (exit breakend, entry breakend, hap1, hap2); a linear
    walk has no join from its last step back to its first."""
    result = []
    n = len(steps)
    for i, s in enumerate(steps):
        if not circular and i == n - 1:
            break
        t = steps[(i + 1) % n]
        if _contiguous(s, t, tol):
            continue
        exit_be = (s.chrom, s.end, "+") if s.strand == "+" else (s.chrom, s.start, "-")
        entry_be = (t.chrom, t.start, "-") if t.strand == "+" else (t.chrom, t.end, "+")
        result.append((exit_be, entry_be, s.hap, t.hap))
    return result


def _be_match(x: Breakend, y: Breakend, tol: int) -> bool:
    return x[0] == y[0] and x[2] == y[2] and abs(x[1] - y[1]) <= tol


def _junction_match(j1, j2, tol: int, phased: bool) -> bool:
    (a1, b1, h1, k1), (a2, b2, h2, k2) = j1, j2
    if _be_match(a1, a2, tol) and _be_match(b1, b2, tol):
        return not phased or (h1, k1) == (h2, k2)
    if _be_match(a1, b2, tol) and _be_match(b1, a2, tol):
        return not phased or (h1, k1) == (k2, h2)
    return False


def junction_scores(truth: Sequence[CStep], inferred: Sequence[CStep], tol: int, phased: bool,
                    truth_circular: bool = True, inferred_circular: bool = True) -> Tuple[float, float]:
    """(recall, precision) of distinct junctions."""
    jt = _unique(junctions(truth, circular=truth_circular), tol)
    ji = _unique(junctions(inferred, circular=inferred_circular), tol)
    if not jt and not ji:
        return 1.0, 1.0
    recall = sum(any(_junction_match(t, i, tol, phased) for i in ji) for t in jt) / len(jt) if jt else 1.0
    precision = sum(any(_junction_match(i, t, tol, phased) for t in jt) for i in ji) / len(ji) if ji else 1.0
    return recall, precision


def _unique(js, tol):
    out = []
    for j in js:
        if not any(_junction_match(j, k, 0, True) for k in out):
            out.append(j)
    return out


def _flip(s: CStep) -> CStep:
    return CStep(s.chrom, s.start, s.end, s.hap, "-" if s.strand == "+" else "+")


def _runs_equal(a: Sequence[CStep], b: Sequence[CStep], tol: int) -> bool:
    return len(a) == len(b) and all(x.chrom == y.chrom and x.hap == y.hap and x.strand == y.strand and
                                    abs(x.start - y.start) <= tol and abs(x.end - y.end) <= tol for x, y in zip(a, b))


def order_match(truth: Sequence[CStep], inferred: Sequence[CStep], tol: int,
                truth_circular: bool = True, inferred_circular: bool = True) -> bool:
    """Same run sequence up to reversal (and rotation, for circles); topologies must agree."""
    if truth_circular != inferred_circular:
        return False
    rt, ri = runs(truth, circular=truth_circular), runs(inferred, circular=inferred_circular)
    if len(rt) != len(ri):
        return False
    candidates = [ri, [_flip(s) for s in reversed(ri)]]
    for seq in candidates:
        for r in (range(len(seq)) if truth_circular else (0,)):
            if _runs_equal(rt, seq[r:] + seq[:r], tol):
                return True
    return False


# ---------------------------------------------------------------------------------------- #
# reading
# ---------------------------------------------------------------------------------------- #

def _parse_step(token: str) -> CStep:
    chrom, span, hap, strand = token.rsplit(":", 3)
    start, end = span.split("-")
    return CStep(chrom, int(start), int(end), hap, strand)


def read_inferred(decompositions_tsv: str, include_low_weight: bool = False) -> Dict[int, Dict[str, object]]:
    """{rank: {"objective": float, "cycles": [Cycle], "n_low_weight": int}} from ``decompositions.tsv``.

    Species flagged ``low_weight`` are left out unless ``include_low_weight``.
    """
    frame = pd.read_csv(decompositions_tsv, sep="\t")
    result: Dict[int, Dict[str, object]] = {}
    for rank, group in frame.groupby("decomposition_rank", sort=True):
        low = 0
        if "status" in group:
            low = int(group.loc[group["status"] == "low_weight", "species"].nunique())
        result[int(rank)] = {"objective": float(group["objective"].iloc[0]), "cycles": [], "n_low_weight": low}
    if "status" in frame and not include_low_weight:
        frame = frame[frame["status"] != "low_weight"]
    for (rank, species), group in frame.groupby(["decomposition_rank", "species"], sort=True):
        group = group.sort_values("ordering_rank")
        orderings = [[_parse_step(t) for t in str(w).split(",")] for w in group["walk"]]
        entry = result[int(rank)]
        circular = str(group["topology"].iloc[0]) != "linear" if "topology" in group else True
        weight_int = None
        if "weight_int" in group and pd.notna(group["weight_int"].iloc[0]):
            weight_int = float(group["weight_int"].iloc[0])
        entry["cycles"].append(Cycle(name=str(species), weight=float(group["weight"].iloc[0]), orderings=orderings, circular=circular,
                                     weight_int=weight_int))
    return result


def read_truth_table(path: str) -> List[Cycle]:
    frame = pd.read_csv(path, sep="\t")
    cycles = []
    for species, group in frame.groupby("species", sort=False):
        group = group.sort_values("order")
        steps = [CStep(str(r["chrom"]), int(r["start"]), int(r["end"]), str(r["hap"]), str(r["strand"])) for r in group.to_dict("records")]
        cycles.append(Cycle(name=str(species), weight=float(group["weight"].iloc[0]), orderings=[steps]))
    return cycles


def write_truth_table(path: str, cycles: Sequence[Cycle]) -> None:
    rows = []
    for c in cycles:
        for i, s in enumerate(c.steps, start=1):
            rows.append({"species": c.name, "order": i, "chrom": s.chrom, "start": s.start, "end": s.end, "strand": s.strand,
                         "hap": s.hap, "weight": c.weight})
    pd.DataFrame(rows, columns=["species", "order", "chrom", "start", "end", "strand", "hap", "weight"]).to_csv(path, sep="\t", index=False)


def _read_walks(path: str) -> Dict[str, List[Tuple[int, str]]]:
    walks, current = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:]
                name = header.split()[0]
                current = name if "circular=true" in header else None
                continue
            if current is None:
                continue
            steps = []
            for token in line.split(","):
                m = re.fullmatch(r"seg(\d+)([+-])", token.strip())
                if m is None:
                    raise ValueError("walk {w}: step {t} is not a species segment".format(w=current, t=token))
                steps.append((int(m.group(1)), m.group(2)))
            walks[current] = steps
            current = None
    return walks


def _read_csv_opt(path: str) -> Optional[pd.DataFrame]:
    return pd.read_csv(path, sep="\t") if os.path.exists(path) else None


def resolve_sample(props: Optional[pd.DataFrame], sample: Optional[str]) -> Optional[str]:
    """Map a sample name (``tumor1``, or the windows name ``sim_tumor1``) onto proportions.tsv."""
    if props is None or "SAMPLE" not in props:
        return sample
    names = [str(s) for s in props["SAMPLE"].unique() if str(s) not in ("normal", ".", "nan")]
    if sample:
        if sample in names:
            return sample
        hits = [n for n in names if sample.endswith(n) or n in sample]
        if hits:
            return max(hits, key=len)
        log.warning("truth sample %s not found in proportions.tsv (have %s); using %s", sample, names, names[0] if names else None)
    return names[0] if names else None


def species_weights(props: Optional[pd.DataFrame], ecdna: Optional[pd.DataFrame], sample: Optional[str]) -> Dict[str, float]:
    """Copies per cell of each species in ``sample`` (bulk average), else the mean over samples."""
    weights: Dict[str, float] = {}
    if ecdna is not None and "PHI" in ecdna:
        for r in ecdna.to_dict("records"):
            weights[str(r["SPECIES"])] = float(r["PHI"])
    if props is not None and {"KIND", "SAMPLE", "ID", "PHI"} <= set(props.columns):
        rows = props[(props["KIND"] == "ecdna_burden") & (props["SAMPLE"].astype(str) == str(sample))]
        for r in rows.to_dict("records"):
            if pd.notna(r["PHI"]):
                weights[str(r["ID"])] = float(r["PHI"])
    return weights


def species_junctions(adj_clone: Optional[pd.DataFrame], adjs: Optional[pd.DataFrame]) -> Dict[str, Dict[frozenset, Tuple[Tuple[str, str], str]]]:
    """{species: {frozenset(EID pair): ((EID-1, EID-2), PHASING)}} for the ecDNA junctions."""
    result: Dict[str, Dict[frozenset, Tuple[Tuple[str, str], str]]] = {}
    if adj_clone is None or adjs is None:
        return result
    eid_pair = {str(r["AID"]): (str(r["EID-1"]), str(r["EID-2"])) for r in adjs.to_dict("records")}
    for r in adj_clone.to_dict("records"):
        species = str(r.get("SPECIES", "."))
        if species in (".", "nan", ""):
            continue
        pair = eid_pair.get(str(r["AID"]))
        if pair:
            result.setdefault(species, {})[frozenset(pair)] = (pair, str(r["PHASING"]))
    return result


def _majority_hap(junction_phasing) -> str:
    counts: Dict[str, int] = {}
    for _, phasing in junction_phasing.values():
        for h in phasing:
            if h in ("A", "B"):
                counts[h] = counts.get(h, 0) + 1
    return max(sorted(counts), key=lambda h: counts[h]) if counts else "A"


def read_aa_cycles(path: str) -> List[Tuple[List[dict], List[Tuple[int, str]]]]:
    """AmpliconArchitect / ecSimulator cycles file -> [(segment rows, walk)].

    Segment rows use the simulator's 0-based half-open coordinates (``START = start - 1``,
    ``END = end``), the same convention as truth/segments.tsv, so step extremity names match
    obs/adjacency.tsv. The walk is a list of (row index, strand) with 1-based indices.
    """
    segments: Dict[str, Tuple[str, int, int]] = {}
    cycles = []
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "Segment" and len(fields) >= 5:
                segments[fields[1]] = (fields[2], int(fields[3]) - 1, int(fields[4]))
            elif line.startswith("Cycle="):
                parts = dict(p.split("=", 1) for p in line.strip().split(";") if "=" in p)
                tokens = [t.strip() for t in parts.get("Segments", "").split(",") if t.strip()]
                cycles.append(tokens)
    result = []
    for tokens in cycles:
        rows, walk, index = [], [], {}
        for t in tokens:
            sid, strand = t[:-1], t[-1]
            if strand not in "+-" or sid not in segments or sid == "0":
                raise ValueError("{p}: cycle step {t} is not a listed segment".format(p=path, t=t))
            if sid not in index:
                chrom, start, end = segments[sid]
                rows.append({"#CHR": chrom, "START": start, "END": end})
                index[sid] = len(rows)
            walk.append((index[sid], strand))
        if walk:
            result.append((rows, walk))
    return result


def truth_from_simulation(run_dir: str, sample: Optional[str] = None, max_permutations: int = 50_000) -> List[Cycle]:
    """Rebuild coordinate- and haplotype-resolved truth cycles from a simulator run directory.

    For each species (truth/ecdna.tsv, else the files in truth/aa/), the walk with coordinates
    comes from ``truth/aa/<species>_cycles.txt`` when it exists; otherwise it is rebuilt from
    ``truth/walks.walk`` + ``truth/segments.tsv`` (older runs). Haplotypes come from the species'
    junction phasings (``truth/adjacency_clone.tsv`` + ``obs/adjacency.tsv``), propagated across
    reference joins. The weight is the species' copies per cell in ``sample``
    (``ecdna_burden`` rows of truth/proportions.tsv), else ecdna.tsv PHI.
    """
    tdir, odir = os.path.join(run_dir, "truth"), os.path.join(run_dir, "obs")
    if not os.path.isdir(tdir):
        raise ValueError("{d} has no truth/ directory".format(d=run_dir))
    props = _read_csv_opt(os.path.join(tdir, "proportions.tsv"))
    ecdna = _read_csv_opt(os.path.join(tdir, "ecdna.tsv"))
    adj_clone = _read_csv_opt(os.path.join(tdir, "adjacency_clone.tsv"))
    adjs = _read_csv_opt(os.path.join(odir, "adjacency.tsv"))
    sample = resolve_sample(props, sample)
    weights = species_weights(props, ecdna, sample)
    junctions_by_species = species_junctions(adj_clone, adjs)
    aa_dir = os.path.join(tdir, "aa")

    names: List[str] = []
    if ecdna is not None and "SPECIES" in ecdna:
        names = [str(n) for n in ecdna["SPECIES"]]
    elif os.path.isdir(aa_dir):
        names = sorted(f[:-len("_cycles.txt")] for f in os.listdir(aa_dir) if f.endswith("_cycles.txt"))
    walks = None
    cycles: List[Cycle] = []
    for name in names:
        jp = junctions_by_species.get(name, {})
        default_hap = _majority_hap(jp)
        weight = weights.get(name, float("nan"))
        aa_path = os.path.join(aa_dir, "{n}_cycles.txt".format(n=name))
        if os.path.exists(aa_path):
            parsed = read_aa_cycles(aa_path)
            for k, (seglist, walk) in enumerate(parsed, start=1):
                steps = _label_haplotypes(seglist, walk, jp, default_hap)
                cname = name if len(parsed) == 1 else "{n}.{k}".format(n=name, k=k)
                cycles.append(Cycle(name=cname, weight=weight, orderings=[steps]))
            continue
        if walks is None:
            walks = _read_walks(os.path.join(tdir, "walks.walk"))
        if name not in walks:
            log.warning("species %s has no circular walk in walks.walk and no truth/aa cycles; skipped", name)
            continue
        seglist = _segments_from_profiles(tdir, props, sample, walks[name], ecdna, name, jp, max_permutations)
        if seglist is None:
            raise ValueError("could not rebuild species {s}: no truth/aa/{s}_cycles.txt, and its segments could not be "
                             "mapped onto truth/segments.tsv; pass --truth-cycles instead".format(s=name))
        cycles.append(Cycle(name=name, weight=weight, orderings=[_label_haplotypes(seglist, walks[name], jp, default_hap)]))
    if not cycles:
        raise ValueError("no ecDNA species found under {d}".format(d=tdir))
    return cycles


def _segments_from_profiles(tdir, props, sample, walk, ecdna, name, junction_phasing, max_permutations):
    """Older runs without truth/aa: the species' segments are the segments.tsv rows above the
    clone-mixture background, in the interval order that makes every walk join consistent."""
    segs = _read_csv_opt(os.path.join(tdir, "segments.tsv"))
    clones = _read_csv_opt(os.path.join(tdir, "clones.tsv"))
    if segs is None or clones is None or props is None:
        return None
    clone_rows = props[props["KIND"] == "clone"]
    u = {str(r["ID"]): float(r["U"]) for r in clone_rows[clone_rows["SAMPLE"].astype(str) == str(sample)].to_dict("records")}
    merged = segs.merge(clones, on=["SID", "#CHR", "START", "END"], how="left")
    bg = sum(u_k * (merged["A_" + k].fillna(0) + merged["B_" + k].fillna(0)) for k, u_k in u.items() if "A_" + k in merged)
    merged["excess"] = merged["CN_A"] + merged["CN_B"] - bg
    cand = merged[merged["excess"] > 0.5].sort_values(["#CHR", "START"]).to_dict("records")
    intervals: List[List[dict]] = []
    for row in cand:
        if intervals and intervals[-1][-1]["#CHR"] == row["#CHR"] and intervals[-1][-1]["END"] == row["START"]:
            intervals[-1].append(row)
        else:
            intervals.append([row])
    n = None
    if ecdna is not None and "N_SEGMENTS" in ecdna:
        match = ecdna[ecdna["SPECIES"].astype(str) == name]
        if len(match):
            n = int(match["N_SEGMENTS"].iloc[0])
    n = n or max(idx for idx, _ in walk)
    return _order_intervals(intervals, walk, n, junction_phasing, max_permutations)


def _step_eids(row, strand):
    chrom, start, end = str(row["#CHR"]), int(row["START"]), int(row["END"])
    head, tail = "{c}:{p}+".format(c=chrom, p=end - 1), "{c}:{p}-".format(c=chrom, p=start)
    return (tail, head) if strand == "+" else (head, tail)  # (entry, exit)


def _concordant(e1: str, e2: str) -> bool:
    def parse(e):
        chrom, rest = e.rsplit(":", 1)
        return chrom, int(rest[:-1]), rest[-1]

    (c1, p1, s1), (c2, p2, s2) = parse(e1), parse(e2)
    if c1 != c2 or s1 == s2:
        return False
    plus, minus = (p1, p2) if s1 == "+" else (p2, p1)
    return minus == plus + 1


def _joins(seglist, walk):
    out = []
    for i, (idx, strand) in enumerate(walk):
        nidx, nstrand = walk[(i + 1) % len(walk)]
        exit_eid = _step_eids(seglist[idx - 1], strand)[1]
        entry_eid = _step_eids(seglist[nidx - 1], nstrand)[0]
        out.append((exit_eid, entry_eid))
    return out


def _order_intervals(intervals, walk, n, junction_phasing, max_permutations):
    tried = 0
    for r in range(1, len(intervals) + 1):
        for combo in itertools.combinations(intervals, r):
            if sum(len(iv) for iv in combo) != n:
                continue
            for perm in itertools.permutations(combo):
                tried += 1
                if tried > max_permutations:
                    return None
                seglist = [row for iv in perm for row in iv]
                if all(frozenset(j) in junction_phasing or _concordant(*j) for j in _joins(seglist, walk)):
                    return seglist
    return None


def _label_haplotypes(seglist, walk, junction_phasing, default_hap: str = "A") -> List[CStep]:
    n = len(walk)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    labels: Dict[int, str] = {}
    conflicts = unexplained = 0
    pending = []
    for i, (exit_eid, entry_eid) in enumerate(_joins(seglist, walk)):
        j = (i + 1) % n
        key = frozenset((exit_eid, entry_eid))
        if key in junction_phasing:
            (e1, e2), phasing = junction_phasing[key]
            h_exit, h_entry = (phasing[0], phasing[1]) if exit_eid == e1 else (phasing[1], phasing[0])
            pending.append((i, h_exit))
            pending.append((j, h_entry))
        elif _concordant(exit_eid, entry_eid):
            parent[find(i)] = find(j)
        else:
            unexplained += 1
    for i, h in pending:
        root = find(i)
        if root in labels and labels[root] != h:
            conflicts += 1
        labels.setdefault(root, h)
    if conflicts:
        log.warning("%d haplotype label conflict(s) while labeling truth walk", conflicts)
    if unexplained:
        log.warning("%d truth walk join(s) are neither a listed junction nor a reference join; "
                    "haplotypes there default to %s", unexplained, default_hap)
    steps = []
    for i, (idx, strand) in enumerate(walk):
        row = seglist[idx - 1]
        steps.append(CStep(str(row["#CHR"]), int(row["START"]), int(row["END"]) - 1, labels.get(find(i), default_hap), strand))
    return steps


# ---------------------------------------------------------------------------------------- #
# scoring
# ---------------------------------------------------------------------------------------- #

@dataclass
class PairScore:
    truth: str
    inferred: Optional[str]
    jaccard_hap: float = 0.0
    jaccard_total: float = 0.0
    segment_recall: float = 0.0  # covered truth bases / truth bases (with multiplicity, haplotype ignored)
    segment_precision: float = 0.0  # covered inferred bases / inferred bases
    junction_recall: float = 0.0
    junction_precision: float = 0.0
    junction_f1: float = 0.0
    phased_junction_recall: float = 0.0
    phased_junction_precision: float = 0.0
    phased_junction_f1: float = 0.0
    phasing_accuracy: float = 0.0  # of the bases both cover, the share on the same haplotype
    phasing_accuracy_swap_invariant: float = 0.0  # the same, A/B as given or swapped, whichever is better
    order_match: bool = False
    order_match_ordering: Optional[int] = None
    truth_weight: float = 0.0
    inferred_weight: Optional[float] = None
    cn_rel_err: Optional[float] = None
    inferred_weight_int: Optional[float] = None
    cn_rel_err_int: Optional[float] = None  # the same error with the integer weight (before the refit)
    inferred_topology: Optional[str] = None
    topology_match: Optional[bool] = None
    composite: float = 0.0


def _f1(recall: float, precision: float) -> float:
    return 2 * recall * precision / (recall + precision) if recall + precision > 0 else 0.0


def _swap_hap(steps: Sequence[CStep]) -> List[CStep]:
    other = {"A": "B", "B": "A"}
    return [CStep(s.chrom, s.start, s.end, other.get(s.hap, s.hap), s.strand) for s in steps]


def composite_score(s: PairScore) -> float:
    """Equal-weight mean of the COMPOSITE_TERMS; cn_accuracy is max(0, 1 - cn_rel_err)."""
    if s.inferred is None:
        return 0.0
    cn_accuracy = max(0.0, 1.0 - s.cn_rel_err) if s.cn_rel_err is not None else 0.0
    terms = [s.jaccard_total, s.junction_f1, s.phasing_accuracy, float(s.order_match), cn_accuracy]
    return sum(terms) / len(terms)


def score_pair(t: Cycle, c: Cycle, tol: int) -> PairScore:
    s = PairScore(truth=t.name, inferred=c.name, truth_weight=t.weight, inferred_weight=c.weight)
    s.jaccard_hap = jaccard(t.steps, c.steps, True)
    s.jaccard_total = jaccard(t.steps, c.steps, False)
    s.junction_recall, s.junction_precision = junction_scores(t.steps, c.steps, tol, False, t.circular, c.circular)
    s.phased_junction_recall, s.phased_junction_precision = junction_scores(t.steps, c.steps, tol, True, t.circular, c.circular)
    s.junction_f1 = _f1(s.junction_recall, s.junction_precision)
    s.phased_junction_f1 = _f1(s.phased_junction_recall, s.phased_junction_precision)
    covered, total_t, total_c = _overlap_bp(t.steps, c.steps, False)
    s.segment_recall = covered / total_t if total_t else 0.0
    s.segment_precision = covered / total_c if total_c else 0.0
    if covered:
        same = _overlap_bp(t.steps, c.steps, True)[0]
        swapped = _overlap_bp(t.steps, _swap_hap(c.steps), True)[0]
        s.phasing_accuracy = same / covered
        s.phasing_accuracy_swap_invariant = max(same, swapped) / covered
    s.inferred_topology = "circular" if c.circular else "linear"
    s.topology_match = t.circular == c.circular
    for k, steps in enumerate(c.orderings, start=1):
        if order_match(t.steps, steps, tol, t.circular, c.circular):
            s.order_match, s.order_match_ordering = True, k
            break
    s.cn_rel_err = abs(c.weight - t.weight) / t.weight if t.weight else None
    if c.weight_int is not None:
        s.inferred_weight_int = c.weight_int
        s.cn_rel_err_int = abs(c.weight_int - t.weight) / t.weight if t.weight else None
    s.composite = composite_score(s)
    return s


def score_decomposition(truth: Sequence[Cycle], inferred: Sequence[Cycle], tol: int) -> List[PairScore]:
    """One-to-one greedy matching of truth species to inferred species by haplotype Jaccard."""
    pairs = []
    for t in truth:
        for c in inferred:
            pairs.append((jaccard(t.steps, c.steps, True), t.name, c.name))
    pairs.sort(reverse=True)
    used_t, used_c, chosen = set(), set(), {}
    for score, tn, cn in pairs:
        if tn in used_t or cn in used_c or score <= 0:
            continue
        used_t.add(tn)
        used_c.add(cn)
        chosen[tn] = cn
    by_name = {c.name: c for c in inferred}
    result = []
    for t in truth:
        if t.name in chosen:
            result.append(score_pair(t, by_name[chosen[t.name]], tol))
        else:
            result.append(PairScore(truth=t.name, inferred=None, truth_weight=t.weight))
    return result


Region = Tuple[str, int, int]


def parse_region(text: str) -> Region:
    chrom, span = text.rsplit(":", 1)
    a, b = span.split("-")
    return chrom, int(a), int(b)


def fraction_in_regions(steps: Sequence[CStep], regions: Sequence[Region], pad: int = 0) -> float:
    """Share of a cycle's bases (with multiplicity) that fall inside ``regions``."""
    total = sum(s.end - s.start + 1 for s in steps)
    if total <= 0:
        return 0.0
    inside = 0
    for s in steps:
        for chrom, a, b in regions:
            if chrom == s.chrom:
                lo, hi = max(s.start, a - pad), min(s.end, b + pad)
                if hi >= lo:
                    inside += hi - lo + 1
    return min(inside / total, 1.0)


def evaluate(decompositions_tsv: str, truth: Sequence[Cycle], tol: int = 10_000, include_low_weight: bool = False,
             regions: Optional[Sequence[Region]] = None, min_truth_fraction: float = 0.5) -> Tuple[pd.DataFrame, Dict[str, object]]:
    """Score every decomposition against the truth species that belong to this subgraph.

    With ``regions`` (the subgraph's segment spans), only truth species with at least
    ``min_truth_fraction`` of their bases inside the regions are scored; the others (species
    built by a different amplicon subgraph) are listed in the summary as out of scope instead
    of counting as missed.
    """
    inferred = read_inferred(decompositions_tsv, include_low_weight=include_low_weight)
    species_info = []
    scored = []
    for t in truth:
        frac = fraction_in_regions(t.steps, regions, pad=tol) if regions else 1.0
        in_scope = frac >= min_truth_fraction
        species_info.append({"species": t.name, "weight": t.weight, "length_bp": sum(s.end - s.start + 1 for s in t.steps),
                             "n_steps": len(t.steps), "n_junctions": len(_unique(junctions(t.steps, circular=t.circular), tol)),
                             "fraction_in_subgraph": round(frac, 4), "in_scope": in_scope})
        if in_scope:
            scored.append(t)
    rows = []
    best_order_rank, best_jaccard, best_jaccard_rank = None, -1.0, None
    for rank in sorted(inferred):
        cycles = inferred[rank]["cycles"]
        scores = score_decomposition(scored, cycles, tol)
        mean_j = sum(s.jaccard_hap for s in scores) / len(scores) if scores else 0.0
        if scores and all(s.order_match for s in scores) and best_order_rank is None:
            best_order_rank = rank
        if scores and mean_j > best_jaccard:
            best_jaccard, best_jaccard_rank = mean_j, rank
        for s in scores:
            row = asdict(s)
            row.update({"decomposition_rank": rank, "objective": inferred[rank]["objective"], "n_species_inferred": len(cycles),
                        "n_low_weight": inferred[rank]["n_low_weight"], "n_species_truth": len(scored)})
            rows.append(row)
    frame = pd.DataFrame(rows)
    top = frame[frame["decomposition_rank"] == min(inferred)] if (len(frame) and inferred) else frame
    summary = {
        "n_decompositions": len(inferred),
        "n_truth_species": len(scored),
        "n_truth_species_total": len(truth),
        "truth_species": species_info,
        "first_rank_with_all_orders_matched": best_order_rank,
        "best_mean_jaccard_hap": best_jaccard if scored else None,
        "best_mean_jaccard_rank": best_jaccard_rank,
        "top_decomposition": _json_records(top.drop(columns=["objective"])) if len(top) else [],
        "tol_bp": tol,
        "min_truth_fraction": min_truth_fraction if regions else None,
    }
    matrix = recovery_matrix(frame, species_info, min(inferred) if inferred else None)
    scored_rows = matrix[matrix["in_scope"]]
    summary["mean_composite"] = float(scored_rows["composite"].mean()) if len(scored_rows) else None
    summary["recovery_matrix"] = _json_records(matrix)
    if not scored:
        summary["note"] = "no truth species has at least {f:.0%} of its bases in this subgraph".format(f=min_truth_fraction)
    return frame, summary


def _json_records(frame: pd.DataFrame) -> List[Dict[str, object]]:
    """Rows as plain Python values (no NaN, no numpy scalars), so evaluation.json stays strict JSON."""
    records = []
    for record in frame.to_dict("records"):
        records.append({k: (None if pd.api.types.is_scalar(v) and pd.isna(v) else (v.item() if hasattr(v, "item") else v))
                        for k, v in record.items()})
    return records


RECOVERY_COLUMNS = ["truth", "in_scope", "fraction_in_subgraph", "truth_length_bp", "truth_n_junctions", "truth_weight",
                    "inferred", "inferred_topology", "inferred_weight",
                    "segment_recall", "segment_precision", "segment_jaccard", "haplotype_jaccard",
                    "junction_recall", "junction_precision", "junction_f1",
                    "phased_junction_recall", "phased_junction_precision", "phased_junction_f1",
                    "phasing_accuracy", "phasing_accuracy_swap_invariant", "order_match", "topology_match", "cn_rel_err",
                    "composite", "best_rank", "best_inferred", "best_composite"]

_RENAMED = {"jaccard_total": "segment_jaccard", "jaccard_hap": "haplotype_jaccard"}


def recovery_matrix(frame: pd.DataFrame, species_info: Sequence[Dict[str, object]], top_rank: Optional[int]) -> pd.DataFrame:
    """One row per truth species: its scores against the rank-1 (``top_rank``) decomposition, and the
    decomposition rank whose matched species scores the highest composite (an upper bound over the
    alternatives). Out-of-scope species are listed with ``in_scope`` False and no scores."""
    rows = []
    for info in species_info:
        row = {"truth": info["species"], "in_scope": bool(info["in_scope"]), "fraction_in_subgraph": info["fraction_in_subgraph"],
               "truth_length_bp": info["length_bp"], "truth_n_junctions": info["n_junctions"], "truth_weight": info["weight"]}
        mine = frame[frame["truth"] == info["species"]] if len(frame) else frame
        if info["in_scope"]:
            top = mine[mine["decomposition_rank"] == top_rank] if len(mine) else mine
            if len(top):
                record = top.iloc[0].to_dict()
                for column in RECOVERY_COLUMNS:
                    source = next((k for k, v in _RENAMED.items() if v == column), column)
                    if column not in row and source in record:
                        row[column] = record[source]
            else:
                row["composite"] = 0.0
            if len(mine):
                best = mine.sort_values(["composite", "decomposition_rank"], ascending=[False, True]).iloc[0]
                row.update({"best_rank": int(best["decomposition_rank"]), "best_inferred": best["inferred"],
                            "best_composite": float(best["composite"])})
        rows.append(row)
    matrix = pd.DataFrame(rows, columns=RECOVERY_COLUMNS)
    matrix["best_rank"] = matrix["best_rank"].astype("Int64")
    return matrix


def regions_from_summary(outdir: str) -> Optional[List[Region]]:
    """Subgraph regions recorded by rck-ecdna in summary.json, if present."""
    path = os.path.join(outdir, "summary.json")
    if not os.path.exists(path):
        return None
    try:
        with open(path) as handle:
            texts = json.load(handle)["graph"]["regions"]
        return [parse_region(t) for t in texts]
    except (KeyError, ValueError, TypeError):
        return None


def write_evaluation(outdir: str, frame: pd.DataFrame, summary: Dict[str, object]) -> Dict[str, str]:
    paths = {"table": os.path.join(outdir, "evaluation.tsv"), "summary": os.path.join(outdir, "evaluation.json"),
             "recovery_matrix": os.path.join(outdir, "recovery_matrix.tsv")}
    frame.to_csv(paths["table"], sep="\t", index=False)
    matrix = pd.DataFrame(summary.get("recovery_matrix", []), columns=RECOVERY_COLUMNS)
    matrix["best_rank"] = matrix["best_rank"].astype("Int64")
    matrix.to_csv(paths["recovery_matrix"], sep="\t", index=False)
    with open(paths["summary"], "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    return paths


def load_truth(truth_dir: Optional[str] = None, truth_cycles: Optional[str] = None, sample: Optional[str] = None) -> List[Cycle]:
    if truth_cycles:
        return read_truth_table(truth_cycles)
    if truth_dir:
        return truth_from_simulation(truth_dir, sample=sample)
    raise ValueError("give --truth-dir or --truth-cycles")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="rck-ecdna-eval", description="Score rck-ecdna output against ground truth.")
    parser.add_argument("--decompositions", required=True, help="decompositions.tsv written by rck-ecdna")
    parser.add_argument("--truth-dir", default=None, help="Simulator run directory (with truth/ and obs/)")
    parser.add_argument("--truth-cycles", default=None, help="Explicit truth table (species order chrom start end strand hap weight)")
    parser.add_argument("--truth-sample", default=None)
    parser.add_argument("--tol-bp", type=int, default=10_000)
    parser.add_argument("--include-low-weight", action="store_true", help="Also score species flagged low_weight")
    parser.add_argument("--min-truth-fraction", type=float, default=0.5,
                        help="Score only truth species with at least this share of their bases inside the subgraph")
    parser.add_argument("--all-truth-species", action="store_true", help="Score every truth species, inside the subgraph or not")
    parser.add_argument("--outdir", default=None, help="Default: the directory of --decompositions")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    truth = load_truth(args.truth_dir, args.truth_cycles, args.truth_sample)
    outdir = args.outdir or os.path.dirname(os.path.abspath(args.decompositions))
    os.makedirs(outdir, exist_ok=True)
    write_truth_table(os.path.join(outdir, "truth_cycles.tsv"), truth)
    regions = None if args.all_truth_species else regions_from_summary(os.path.dirname(os.path.abspath(args.decompositions)))
    frame, summary = evaluate(args.decompositions, truth, tol=args.tol_bp, include_low_weight=args.include_low_weight,
                              regions=regions, min_truth_fraction=args.min_truth_fraction)
    write_evaluation(outdir, frame, summary)
    log_summary(summary)
    return 0


def log_summary(summary: Dict[str, object]) -> None:
    for sp in summary.get("truth_species", []):
        log.info("truth %s: %.0f bp, weight %.2f, %.0f%% inside the subgraph -> %s", sp["species"], sp["length_bp"], sp["weight"],
                 100 * sp["fraction_in_subgraph"], "scored" if sp["in_scope"] else "out of scope")
    if summary.get("note"):
        log.warning("evaluation: %s", summary["note"])
        return
    for row in summary.get("recovery_matrix", []):
        if row["in_scope"]:
            log.info("recovery %s: composite %.3f with %s (best %s at rank %s)", row["truth"], row["composite"] or 0.0,
                     row.get("inferred") or "nothing", "n/a" if row.get("best_composite") is None else "{b:.3f}".format(b=row["best_composite"]),
                     row.get("best_rank"))
    best = summary.get("best_mean_jaccard_hap")
    log.info("evaluation: first rank with every truth species in order: %s; best mean haplotype Jaccard %s (rank %s)",
             summary["first_rank_with_all_orders_matched"], "n/a" if best is None else "{b:.3f}".format(b=best),
             summary["best_mean_jaccard_rank"])


if __name__ == "__main__":
    raise SystemExit(main())
