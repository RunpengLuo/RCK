"""Read an ecDNA breakpoint subgraph and normalize it into a DIAG.

Input is the directory written by amplicon-allele-analysis ``subgraph_construct.py`` for one
amplicon (``nodes.tsv``, ``edges.tsv``, optionally ``phase_blocks.tsv``). Its edge types map
onto RCK's interval adjacency graph as follows:

=================  ======================================================================
``EDGE_TYPE``      DIAG element
=================  ======================================================================
``reference``      segment edge (an interval), weight = mean depth
``extremity``      reference adjacency between the two sides of an extremity, weight =
                   ``#REF`` (reads spanning the reference join)
``variant``        observed novel adjacency, weight = ``#NOV`` (split / spanning reads)
``imputed_variant`` candidate novel adjacency proposed by proximity imputation
=================  ======================================================================

Three normalizations are applied, so the solver sees a clean DIAG in which every
non-terminal extremity has exactly one segment edge:

1. **Phantom segments.** Simulated data names both sides of a breakpoint as separate
   extremities (``x+`` and ``x+1-``), which produces a 1 bp segment between them that only
   carries reference flow. Each such segment is contracted: its two reference adjacencies
   become one, the base is absorbed into the right-hand segment, and any other edge attached
   to the phantom's tail (head) is moved to the true tail (head) on the same side.
2. **Parallel imputed edges.** An imputed edge that joins the same two extremities as an
   existing reference adjacency (``inter`` edges at imputed extremities, anchors that loop
   back onto their own extremity) is the reference join itself. It is merged into that
   reference adjacency; its weight is kept only as weak evidence. Imputed edges parallel to
   an observed variant are merged into the variant.
3. **Terminals.** Extremities with no segment edge (the outer flanks of each region),
   chromosome boundary nodes, and segment ends without any reference adjacency are marked
   terminal. Only the chromosomal background layer may start or end there; the ecDNA layer
   is strictly circular.

Allowed haplotype labelings (RCK ``Phasing``) are attached to every adjacency here:
reference adjacencies are AA/BB only; observed junctions follow their read-based
``junction_hap`` label (``A`` -> AA, ``B`` -> BB, ``mixed`` -> AB/BA, i.e. a haplotype
crossing); imputed candidates are AA/BB. Novel self-loops cannot be AB/BA, as in RCK.
"""

from __future__ import annotations

import logging
import math
import os
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

from rck.core.structures import Haplotype, Phasing

log = logging.getLogger(__name__)

HAPLOTYPES = (Haplotype.A, Haplotype.B)
ALL_PHASINGS = (Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB)
SAME_HAP_PHASINGS = (Phasing.AA, Phasing.BB)
CROSS_HAP_PHASINGS = (Phasing.AB, Phasing.BA)

REFERENCE = "reference"
NOVEL = "novel"
IMPUTED = "imputed"

HEAD = "+"  # right end of a segment (the subgraph's "_+" node side)
TAIL = "-"  # left end of a segment (the subgraph's "_-" node side)


@dataclass
class Node:
    """One extremity side of the DIAG (an RCK position)."""

    nid: str
    chrom: str
    pos: int
    side: Optional[str]  # HEAD, TAIL, or None for chromosome boundary nodes
    source: str = ""
    segment: Optional[str] = None  # sid of the incident segment edge, if any
    terminal: bool = False

    @property
    def label(self) -> str:
        return "{c}:{p}{s}".format(c=self.chrom, p=self.pos, s=self.side or "")


@dataclass
class Segment:
    """A segment edge. Coordinates are the inclusive tail and head positions."""

    sid: str
    chrom: str
    start: int
    end: int
    tail: str
    head: str
    depth: float
    baf: Optional[float] = None
    n_snps: int = 0
    seg_hap: Optional[str] = None
    absorbed_bp: int = 0
    cn_total: Optional[float] = None  # cell-averaged copy number written by subgraph_construct.py (CN_TOTAL)
    cn_a: Optional[float] = None  # ... and its allele split (CN_A / CN_B)
    cn_b: Optional[float] = None

    @property
    def length(self) -> int:
        return max(self.end - self.start + 1, 1)

    def other_end(self, nid: str) -> str:
        if nid == self.tail:
            return self.head
        if nid == self.head:
            return self.tail
        raise KeyError("{n} is not an end of segment {s}".format(n=nid, s=self.sid))

    @property
    def label(self) -> str:
        return "{c}:{s}-{e}".format(c=self.chrom, s=self.start, e=self.end)


@dataclass
class Adjacency:
    """An adjacency edge between two extremities, ``u`` and ``v`` (``u == v`` is a self-loop).

    Phasing ``(h1, h2)`` means the copy of ``u`` on haplotype ``h1`` is joined to the copy of
    ``v`` on haplotype ``h2``, as in RCK.
    """

    aid: str
    u: str
    v: str
    kind: str
    reads: Optional[float] = None  # #REF or #NOV; None when unobserved (imputed)
    reads_a: Optional[float] = None
    reads_b: Optional[float] = None
    junction_hap: Optional[str] = None
    allowed: Tuple[Phasing, ...] = SAME_HAP_PHASINGS
    preferred: Tuple[Phasing, ...] = SAME_HAP_PHASINGS
    imputed_support: Optional[float] = None
    sources: List[str] = field(default_factory=list)

    @property
    def is_self_loop(self) -> bool:
        return self.u == self.v

    @property
    def is_novel_like(self) -> bool:
        return self.kind in (NOVEL, IMPUTED)

    def other_end(self, nid: str) -> str:
        if nid == self.u:
            return self.v
        if nid == self.v:
            return self.u
        raise KeyError("{n} is not an end of adjacency {a}".format(n=nid, a=self.aid))

    def hap_at(self, phasing: Phasing, nid: str) -> List[Haplotype]:
        """Haplotypes this phasing occupies at ``nid`` (two entries for a self-loop)."""
        h1, h2 = phasing.value
        result = []
        if nid == self.u:
            result.append(h1)
        if nid == self.v:
            result.append(h2)
        return result


class ECGraph:
    """The normalized DIAG of one amplicon."""

    def __init__(self, nodes: Dict[str, Node], segments: Dict[str, Segment], adjacencies: Dict[str, Adjacency],
                 name: str = "amplicon", log_messages: Optional[List[str]] = None, node_alias: Optional[Dict[str, str]] = None):
        self.name = name
        self.nodes = nodes
        self.segments = segments
        self.adjacencies = adjacencies
        self.log_messages = log_messages or []
        self.node_alias = node_alias or {}
        self._adj_at = defaultdict(list)
        for a in self.adjacencies.values():
            self._adj_at[a.u].append(a.aid)
            if a.v != a.u:
                self._adj_at[a.v].append(a.aid)

    # ---- queries ---- #
    def adjacencies_at(self, nid: str) -> List[Adjacency]:
        return [self.adjacencies[aid] for aid in self._adj_at.get(nid, [])]

    def segment_at(self, nid: str) -> Optional[Segment]:
        sid = self.nodes[nid].segment
        return None if sid is None else self.segments[sid]

    def terminals(self) -> List[str]:
        return sorted(n.nid for n in self.nodes.values() if n.terminal)

    def resolve(self, nid: str) -> Optional[str]:
        """Map an original subgraph node id (possibly a contracted phantom) to a graph node."""
        seen = set()
        while nid in self.node_alias and nid not in seen:
            seen.add(nid)
            nid = self.node_alias[nid]
        return nid if nid in self.nodes else None

    def regions(self, max_gap: int = 1_000_000) -> List[Tuple[str, int, int]]:
        """Spans covered by segments, as (chrom, start, end); segments further apart than
        ``max_gap`` on the same chromosome start a new span."""
        result: List[Tuple[str, int, int]] = []
        for s in sorted(self.segments.values(), key=lambda x: (x.chrom, x.start)):
            if result and result[-1][0] == s.chrom and s.start - result[-1][2] <= max_gap:
                chrom, start, end = result[-1]
                result[-1] = (chrom, start, max(end, s.end))
            else:
                result.append((s.chrom, s.start, s.end))
        return result

    def summary(self) -> Dict[str, object]:
        kinds = defaultdict(int)
        for a in self.adjacencies.values():
            kinds[a.kind] += 1
        return {
            "name": self.name,
            "nodes": len(self.nodes),
            "segments": len(self.segments),
            "adjacencies": dict(kinds),
            "terminals": len(self.terminals()),
            "regions": ["{c}:{s}-{e}".format(c=c, s=s, e=e) for c, s, e in self.regions()],
            "normalization": list(self.log_messages),
        }

    def validate(self) -> List[str]:
        """Structural problems; an empty list means every node is a proper DIAG extremity."""
        problems = []
        seg_count = defaultdict(int)
        for s in self.segments.values():
            seg_count[s.tail] += 1
            seg_count[s.head] += 1
        for nid, node in self.nodes.items():
            if seg_count[nid] > 1:
                problems.append("node {n} has {k} segment edges".format(n=nid, k=seg_count[nid]))
            if seg_count[nid] == 0 and not node.terminal:
                problems.append("node {n} has no segment edge and is not terminal".format(n=nid))
        for a in self.adjacencies.values():
            for x in (a.u, a.v):
                if x not in self.nodes:
                    problems.append("adjacency {a} references unknown node {n}".format(a=a.aid, n=x))
        return problems


# ---------------------------------------------------------------------------------------- #
# reading
# ---------------------------------------------------------------------------------------- #

def _num(value, default=None):
    try:
        if value is None:
            return default
        f = float(value)
        if math.isnan(f):
            return default
        return f
    except (TypeError, ValueError):
        return default


def _str(value) -> Optional[str]:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    s = str(value).strip()
    return s if s and s.lower() != "nan" else None


def _junction_phasings(junction_hap: Optional[str], self_loop: bool, soft: bool) -> Tuple[Tuple[Phasing, ...], Tuple[Phasing, ...]]:
    """(allowed, preferred) labelings for an observed junction."""
    label = (junction_hap or "").strip().lower()
    if label == "a":
        preferred = (Phasing.AA,)
    elif label == "b":
        preferred = (Phasing.BB,)
    elif label == "mixed":
        preferred = CROSS_HAP_PHASINGS
    else:
        preferred = ALL_PHASINGS
    allowed = ALL_PHASINGS if soft else preferred
    if self_loop:
        allowed = tuple(p for p in allowed if p in SAME_HAP_PHASINGS)
        preferred = tuple(p for p in preferred if p in SAME_HAP_PHASINGS)
        if not allowed:
            # a "mixed" self-loop cannot be realized; fall back to same-haplotype labelings
            allowed = SAME_HAP_PHASINGS
        if not preferred:
            preferred = allowed
    return allowed, preferred


def read_subgraph(subgraph_dir: str, phantom_max_bp: int = 2, soft_junction_hap: bool = False,
                  name: Optional[str] = None) -> ECGraph:
    """Read ``nodes.tsv`` and ``edges.tsv`` from ``subgraph_dir`` into a normalized :class:`ECGraph`.

    Args:
        subgraph_dir: Directory with the subgraph tables of one amplicon.
        phantom_max_bp: Segments of at most this many bp whose two ends each have exactly one
            reference adjacency are contracted; other edges on them are moved to the true
            extremities (see module docstring). 0 disables contraction.
        soft_junction_hap: When True every labeling is allowed for observed junctions and the
            read-based label only sets the preferred ones (penalized in the MILP); when False
            (default) the read-based label is a hard constraint.
        name: Name used in outputs; defaults to the directory name.
    """
    nodes_path = os.path.join(subgraph_dir, "nodes.tsv")
    edges_path = os.path.join(subgraph_dir, "edges.tsv")
    nodes_df = pd.read_csv(nodes_path, sep="\t", dtype={"NODE_ID": str})
    edges_df = pd.read_csv(edges_path, sep="\t", dtype={"NODE_1": str, "NODE_2": str})
    return build_graph(nodes_df, edges_df, phantom_max_bp=phantom_max_bp, soft_junction_hap=soft_junction_hap,
                       name=name or os.path.basename(os.path.normpath(subgraph_dir)))


def build_graph(nodes_df: pd.DataFrame, edges_df: pd.DataFrame, phantom_max_bp: int = 2, soft_junction_hap: bool = False,
                name: str = "amplicon") -> ECGraph:
    """Build an :class:`ECGraph` from subgraph tables already loaded as data frames."""
    messages: List[str] = []
    nodes: Dict[str, Node] = {}
    for row in nodes_df.to_dict("records"):
        nid = str(row["NODE_ID"])
        orient = _str(row.get("ORIENT"))
        side = orient if orient in (HEAD, TAIL) else None
        nodes[nid] = Node(nid=nid, chrom=str(row["#CHR"]), pos=int(row["POS"]), side=side, source=_str(row.get("SOURCE")) or "")

    def ensure_node(nid, chrom, pos):
        if nid not in nodes:
            side = HEAD if nid.endswith("_+") else TAIL if nid.endswith("_-") else None
            nodes[nid] = Node(nid=nid, chrom=str(chrom), pos=int(pos), side=side, source="edge_only")
            messages.append("node {n} appears only in edges.tsv; added".format(n=nid))

    segments: Dict[str, Segment] = {}
    raw_adjacencies: List[Adjacency] = []
    counter = defaultdict(int)

    for row in edges_df.to_dict("records"):
        etype = str(row["EDGE_TYPE"]).strip()
        n1, n2 = str(row["NODE_1"]), str(row["NODE_2"])
        ensure_node(n1, row["CHROM_1"], row["POS_1"])
        ensure_node(n2, row["CHROM_2"], row["POS_2"])
        weight = _num(row.get("WEIGHT"), 0.0)
        if etype == "reference":
            a, b = (n1, n2) if nodes[n1].pos <= nodes[n2].pos else (n2, n1)
            counter["s"] += 1
            sid = "s{i:04d}".format(i=counter["s"])
            baf = _num(row.get("BAF"))
            n_snps = int(_num(row.get("n_snps"), 0) or 0)
            segments[sid] = Segment(sid=sid, chrom=nodes[a].chrom, start=nodes[a].pos, end=nodes[b].pos, tail=a, head=b,
                                    depth=weight, baf=baf if n_snps > 0 else None, n_snps=n_snps, seg_hap=_str(row.get("seg_hap")),
                                    cn_total=_num(row.get("CN_TOTAL")), cn_a=_num(row.get("CN_A")), cn_b=_num(row.get("CN_B")))
        elif etype == "extremity":
            counter["r"] += 1
            imputed_ext = n1.startswith("imputed_ext") or n2.startswith("imputed_ext")
            raw_adjacencies.append(Adjacency(aid="ref{i:04d}".format(i=counter["r"]), u=n1, v=n2, kind=REFERENCE,
                                             reads=None if imputed_ext else weight, sources=["extremity:{a}-{b}".format(a=n1, b=n2)]))
        elif etype == "variant":
            aid = _str(row.get("AID")) or "nov{i:04d}".format(i=counter["n"] + 1)
            counter["n"] += 1
            jhap = _str(row.get("junction_hap"))
            allowed, preferred = _junction_phasings(jhap, n1 == n2, soft_junction_hap)
            raw_adjacencies.append(Adjacency(aid=aid, u=n1, v=n2, kind=NOVEL, reads=weight,
                                             reads_a=_num(row.get("NOV_A")), reads_b=_num(row.get("NOV_B")),
                                             junction_hap=jhap, allowed=allowed, preferred=preferred, sources=[aid]))
        elif etype == "imputed_variant":
            aid = _str(row.get("AID")) or "imp{i:04d}".format(i=counter["i"] + 1)
            counter["i"] += 1
            raw_adjacencies.append(Adjacency(aid=aid, u=n1, v=n2, kind=IMPUTED, reads=None, imputed_support=weight,
                                             junction_hap="imputed", sources=[aid]))
        else:
            messages.append("ignored edge {a}-{b} of unknown type {t}".format(a=n1, b=n2, t=etype))

    alias: Dict[str, str] = {}
    segments, raw_adjacencies = _contract_phantoms(nodes, segments, raw_adjacencies, phantom_max_bp, alias, messages)
    adjacencies = _merge_parallel(raw_adjacencies, messages)

    # attach segments to nodes
    for node in nodes.values():
        node.segment = None
    for s in segments.values():
        nodes[s.tail].segment = s.sid
        nodes[s.head].segment = s.sid

    # drop nodes no longer used by anything
    used = set()
    for s in segments.values():
        used.update((s.tail, s.head))
    for a in adjacencies.values():
        used.update((a.u, a.v))
    for nid in list(nodes):
        if nid not in used:
            del nodes[nid]

    # terminals
    ref_at = defaultdict(int)
    for a in adjacencies.values():
        if a.kind == REFERENCE:
            ref_at[a.u] += 1
            ref_at[a.v] += 1
    for node in nodes.values():
        boundary = node.source == "chromosome_boundary" or node.side is None
        node.terminal = node.segment is None or boundary or ref_at[node.nid] == 0

    graph = ECGraph(nodes=nodes, segments=segments, adjacencies=adjacencies, name=name, log_messages=messages, node_alias=alias)
    for problem in graph.validate():
        log.warning("subgraph %s: %s", name, problem)
        messages.append("WARNING " + problem)
    return graph


def _contract_phantoms(nodes, segments, adjacencies, phantom_max_bp, alias, messages):
    """Contract 1-2 bp segments that only connect two reference adjacencies."""
    if phantom_max_bp <= 0:
        return segments, adjacencies
    refs_at = defaultdict(list)
    for a in adjacencies:
        if a.kind == REFERENCE:
            refs_at[a.u].append(a)
            refs_at[a.v].append(a)
    seg_at = {}
    for s in segments.values():
        seg_at[s.tail] = s.sid
        seg_at[s.head] = s.sid

    removed_adjs = set()
    new_adjs: List[Adjacency] = []
    contracted = 0
    for sid in sorted(segments, key=lambda k: (segments[k].chrom, segments[k].start)):
        s = segments[sid]
        if s.length > phantom_max_bp:
            continue
        u_refs = [a for a in refs_at[s.tail] if id(a) not in removed_adjs]
        v_refs = [a for a in refs_at[s.head] if id(a) not in removed_adjs]
        if len(u_refs) != 1 or len(v_refs) != 1:
            continue
        ru, rv = u_refs[0], v_refs[0]
        p = ru.other_end(s.tail)  # head of the left segment
        q = rv.other_end(s.head)  # tail of the right segment
        if p == s.head or q == s.tail or p not in seg_at or q not in seg_at:
            continue
        removed_adjs.update((id(ru), id(rv)))
        reads = [r for r in (ru.reads, rv.reads) if r is not None]
        merged = Adjacency(aid="ref_c{i:03d}".format(i=contracted + 1), u=p, v=q, kind=REFERENCE,
                           reads=(sum(reads) / len(reads)) if reads else None, sources=ru.sources + rv.sources)
        new_adjs.append(merged)
        refs_at[p] = [a for a in refs_at[p] if a is not ru] + [merged]
        refs_at[q] = [a for a in refs_at[q] if a is not rv] + [merged]
        # the phantom base joins the right-hand segment; edges on the phantom tail go to the
        # true tail q, edges on the phantom head go to the true head p
        right = segments[seg_at[q]]
        if right.start > s.start:
            right.absorbed_bp += right.start - s.start
            right.start = s.start
        alias[s.tail] = q
        alias[s.head] = p
        del segments[sid]
        del seg_at[s.tail]
        del seg_at[s.head]
        contracted += 1
    if contracted:
        messages.append("contracted {k} phantom segment(s) of <= {b} bp".format(k=contracted, b=phantom_max_bp))

    result = []
    for a in adjacencies:
        if id(a) in removed_adjs:
            continue
        u, v = alias.get(a.u, a.u), alias.get(a.v, a.v)
        if (u, v) != (a.u, a.v):
            messages.append("moved {k} {aid} from phantom extremity to {u}--{v}".format(k=a.kind, aid=a.aid, u=u, v=v))
            a.u, a.v = u, v
            if a.kind == NOVEL:
                # self-loop status may have changed; recompute, keeping soft mode if it was on
                was_soft = a.allowed == ALL_PHASINGS and a.preferred != ALL_PHASINGS
                a.allowed, a.preferred = _junction_phasings(a.junction_hap, a.u == a.v, was_soft)
        result.append(a)
    result.extend(new_adjs)
    return segments, result


def _merge_parallel(adjacencies: Iterable[Adjacency], messages: List[str]) -> Dict[str, Adjacency]:
    """Merge imputed edges into a reference adjacency or variant on the same node pair."""
    by_pair = defaultdict(list)
    for a in adjacencies:
        by_pair[frozenset((a.u, a.v))].append(a)
    result: Dict[str, Adjacency] = {}
    for pair, group in by_pair.items():
        refs = [a for a in group if a.kind == REFERENCE]
        novs = [a for a in group if a.kind == NOVEL]
        imps = [a for a in group if a.kind == IMPUTED]
        if len(refs) > 1:
            keep = refs[0]
            reads = [r.reads for r in refs if r.reads is not None]
            keep.reads = (sum(reads) / len(reads)) if reads else None
            for r in refs[1:]:
                keep.sources.extend(r.sources)
            messages.append("merged {k} parallel reference adjacencies into {a}".format(k=len(refs), a=keep.aid))
            refs = [keep]
        if len(novs) > 1:
            keep = novs[0]
            keep.reads = sum((n.reads or 0.0) for n in novs)
            for n in novs[1:]:
                keep.sources.extend(n.sources)
            messages.append("merged {k} parallel variants into {a}".format(k=len(novs), a=keep.aid))
            novs = [keep]
        targets = refs or novs
        if imps and targets:
            target = targets[0]
            support = [i.imputed_support for i in imps if i.imputed_support is not None]
            target.imputed_support = max(support) if support else None
            target.sources.extend(s for i in imps for s in i.sources)
            messages.append("merged imputed {ids} into {k} {a}".format(ids=",".join(i.aid for i in imps), k=target.kind, a=target.aid))
            imps = []
        elif len(imps) > 1:
            keep = imps[0]
            support = [i.imputed_support for i in imps if i.imputed_support is not None]
            keep.imputed_support = max(support) if support else None
            for i in imps[1:]:
                keep.sources.extend(i.sources)
            imps = [keep]
        for a in refs + novs + imps:
            if a.aid in result:
                a.aid = "{a}_{i}".format(a=a.aid, i=len(result))
            result[a.aid] = a
    return result
