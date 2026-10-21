"""Turn a balanced phased edge multiset into ordered circular or linear walks.

A walk alternates segment and adjacency edges: after arriving at an extremity on an
adjacency, it must leave on that extremity's segment (and vice versa). Because every phased
node has exactly one segment edge and a balanced multiset gives it as many adjacency ends as
segment copies, a connected balanced multiset always has such an alternating closed walk.

* :func:`alternating_circuit` - Hierholzer's algorithm with alternation: build a greedy
  closed trail, then splice sub-trails in before any segment traversal whose segment still
  has copies left. Splicing before a segment step keeps alternation intact.
* :func:`enumerate_circuits` - depth-first enumeration of distinct circuits (up to rotation and
  reversal) for the same multiset. Different circuits have identical copy numbers, so read
  depth cannot tell them apart; long-read paths can (:func:`path_support`).

A linear species (two linear ends, see :mod:`rck.ecdna.phased`) is handled by joining its two
ends with a pseudo segment, extracting a closed circuit that starts on that pseudo segment, and
dropping it again: the result is an open walk whose first step carries ``entry`` (the end it
starts at) and whose last step's ``join`` is the end it stops at.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from rck.ecdna.phased import (ADJ, END, END_A, END_B, PSEUDO, PSEUDO_KEY, SEG, EdgeKey, PhasedNode, endpoints, incidence,
                              segment_key_at)
from rck.ecdna.subgraph_io import REFERENCE, ECGraph


@dataclass(frozen=True)
class Step:
    """One segment traversal and the adjacency taken right after it."""

    sid: str
    hap: str
    strand: str  # "+" tail -> head, "-" head -> tail
    join: EdgeKey  # ("J", aid, phasing) leading to the next step; ("T", nid, hap) ends a linear walk
    entry: Optional[EdgeKey] = None  # first step of a linear walk: ("T", nid, hap) it starts at


def is_linear(steps: Sequence[Step]) -> bool:
    return bool(steps) and steps[-1].join[0] == END


def _flip(strand: str) -> str:
    return "-" if strand == "+" else "+"


def _other(graph: ECGraph, key: EdgeKey, node: PhasedNode) -> PhasedNode:
    x, y = endpoints(graph, key)
    return y if node == x else x


class _State:
    def __init__(self, graph: ECGraph, counts: Dict[EdgeKey, int]):
        self.graph = graph
        self.rem = {k: int(v) for k, v in counts.items() if v > 0}
        inc = incidence(graph, self.rem)
        self.inc = {n: sorted(v) for n, v in inc.items()}

    def seg_at(self, node: PhasedNode) -> Optional[EdgeKey]:
        return segment_key_at(self.graph, node)


def _trail(state: _State, start: PhasedNode) -> List[Tuple[str, EdgeKey, PhasedNode, PhasedNode]]:
    items = []
    node = start
    while True:
        sk = state.seg_at(node)
        if sk is None or state.rem.get(sk, 0) <= 0:
            if node != start:
                raise ValueError("edge multiset is not balanced at {n}".format(n=node))
            break
        state.rem[sk] -= 1
        other = _other(state.graph, sk, node)
        items.append((SEG, sk, node, other))
        choice = next((ak for ak in state.inc.get(other, []) if state.rem.get(ak, 0) > 0), None)
        if choice is None:
            raise ValueError("edge multiset is not balanced at {n}".format(n=other))
        state.rem[choice] -= 1
        nxt = _other(state.graph, choice, other)
        items.append((ADJ, choice, other, nxt))
        node = nxt
    return items


def _as_end(key: EdgeKey) -> EdgeKey:
    return (END, key[1], key[2]) if key[0] in (END_A, END_B) else key


def _to_steps(graph: ECGraph, items) -> List[Step]:
    steps = []
    entry = None
    for i in range(0, len(items), 2):
        _, sk, frm, _ = items[i]
        _, ak, _, _ = items[i + 1]
        if sk[0] == PSEUDO:
            entry = _as_end(ak)  # the end the linear walk starts at
            continue
        tail = (graph.segments[sk[1]].tail, sk[2])
        steps.append(Step(sid=sk[1], hap=sk[2], strand="+" if frm == tail else "-", join=_as_end(ak),
                          entry=entry if not steps else None))
    return steps


def _close_linear(counts: Dict[EdgeKey, int]) -> Tuple[Dict[EdgeKey, int], bool]:
    """Join the two linear ends of a multiset with a pseudo segment; (counts, is_linear)."""
    ends = []
    for key in sorted(counts):
        if key[0] == END:
            ends.extend([key] * int(counts[key]))
    if not ends:
        return counts, False
    if len(ends) != 2:
        raise ValueError("a linear walk has exactly two ends, got {n}".format(n=len(ends)))
    closed = {k: v for k, v in counts.items() if k[0] != END}
    a, b = ends
    closed[(END_A, a[1], a[2])] = 1
    closed[(END_B, b[1], b[2])] = 1
    closed[PSEUDO_KEY] = 1
    return closed, True


def alternating_circuit(graph: ECGraph, counts: Dict[EdgeKey, int], start: Optional[EdgeKey] = None) -> List[Step]:
    """One alternating walk using every edge of ``counts`` exactly its multiplicity times:
    closed for a circular multiset, open from end to end for a linear one."""
    counts, linear = _close_linear(counts)
    state = _State(graph, counts)
    seg_keys = sorted(k for k in state.rem if k[0] == SEG)
    if not seg_keys:
        return []
    start_key = PSEUDO_KEY if linear else (start or seg_keys[0])
    circuit = _trail(state, endpoints(graph, start_key)[0])
    i = 0
    while i < len(circuit):
        kind, key, frm, _ = circuit[i]
        if kind == SEG and state.rem.get(key, 0) > 0:
            circuit[i:i] = _trail(state, frm)
            continue
        i += 1
    if any(v > 0 for v in state.rem.values()):
        raise ValueError("edge multiset is not connected; it is more than one cycle")
    return _to_steps(graph, circuit)


def enumerate_circuits(graph: ECGraph, counts: Dict[EdgeKey, int], limit: int = 20, max_expansions: int = 200_000) -> List[List[Step]]:
    """Distinct alternating walks (up to rotation / reversal for circles, reversal for linear
    walks), at most ``limit`` of them."""
    counts, linear = _close_linear(counts)
    state = _State(graph, counts)
    seg_keys = sorted(k for k in state.rem if k[0] == SEG)
    if not seg_keys:
        return []
    start = endpoints(graph, PSEUDO_KEY if linear else seg_keys[0])[0]
    remaining = [sum(state.rem.values())]
    path: List = []
    found: Dict[Tuple, List[Step]] = {}
    budget = [max_expansions]

    def dfs(node: PhasedNode):
        if budget[0] <= 0 or len(found) >= limit:
            return
        sk = state.seg_at(node)
        if sk is None or state.rem.get(sk, 0) <= 0:
            return
        state.rem[sk] -= 1
        other = _other(graph, sk, node)
        tried = set()
        for ak in state.inc.get(other, []):
            if state.rem.get(ak, 0) <= 0 or ak in tried:
                continue
            tried.add(ak)
            budget[0] -= 1
            state.rem[ak] -= 1
            nxt = _other(graph, ak, other)
            path.append((SEG, sk, node, other))
            path.append((ADJ, ak, other, nxt))
            remaining[0] -= 2
            if remaining[0] == 0:
                if nxt == start:
                    steps = _to_steps(graph, path)
                    found.setdefault(canonical(steps), steps)
            else:
                dfs(nxt)
            remaining[0] += 2
            path.pop()
            path.pop()
            state.rem[ak] += 1
            if budget[0] <= 0 or len(found) >= limit:
                break
        state.rem[sk] += 1

    dfs(start)
    return [found[k] for k in sorted(found)]


# ---- canonical form ---- #

def _tokens(steps: Sequence[Step]) -> List[Tuple[str, str, str, str, str]]:
    return [(s.sid, s.hap, s.strand, s.join[1], s.join[2]) for s in steps]


def reverse_walk(steps: Sequence[Step]) -> List[Step]:
    n = len(steps)
    if is_linear(steps):
        out = []
        for i in range(n - 1, -1, -1):
            st = steps[i]
            join = steps[i - 1].join if i > 0 else st.entry
            out.append(Step(st.sid, st.hap, _flip(st.strand), join, entry=st.join if i == n - 1 else None))
        return out
    return [Step(steps[i].sid, steps[i].hap, _flip(steps[i].strand), steps[(i - 1) % n].join) for i in range(n - 1, -1, -1)]


def canonical(steps: Sequence[Step]) -> Tuple:
    """Rotation- and reversal-invariant key of a circular walk; reversal-invariant for a linear one."""
    if not steps:
        return ()
    if is_linear(steps):
        forward = tuple(_tokens(steps)) + (("entry",) + tuple(steps[0].entry or ()),)
        rev = reverse_walk(steps)
        backward = tuple(_tokens(rev)) + (("entry",) + tuple(rev[0].entry or ()),)
        return ("linear",) + min(forward, backward)
    best = None
    for seq in (_tokens(steps), _tokens(reverse_walk(steps))):
        for r in range(len(seq)):
            cand = tuple(seq[r:] + seq[:r])
            if best is None or cand < best:
                best = cand
    return best


def walk_counts(steps: Sequence[Step]) -> Dict[EdgeKey, int]:
    counts: Dict[EdgeKey, int] = {}
    for s in steps:
        k = (SEG, s.sid, s.hap)
        counts[k] = counts.get(k, 0) + 1
        counts[s.join] = counts.get(s.join, 0) + 1
        if s.entry is not None:
            counts[s.entry] = counts.get(s.entry, 0) + 1
    return counts


# ---- description ---- #

def describe_step(graph: ECGraph, step: Step) -> str:
    seg = graph.segments[step.sid]
    return "{c}:{s}-{e}:{h}:{st}".format(c=seg.chrom, s=seg.start, e=seg.end, h=step.hap, st=step.strand)


def describe_join(graph: ECGraph, step: Step) -> str:
    if step.join[0] == END:
        return "end:{n}:{h}".format(n=step.join[1], h=step.join[2])
    a = graph.adjacencies[step.join[1]]
    return "{k}:{a}:{p}".format(k=a.kind, a=a.aid, p=step.join[2])


def walk_length(graph: ECGraph, steps: Sequence[Step]) -> int:
    return int(sum(graph.segments[s.sid].length for s in steps))


# ---- long-read path support ---- #

@dataclass
class ReadPath:
    pid: str
    support: float
    junctions: List[Tuple[str, str]]  # (exit node, entry node) for each non-reference join


def _junction_sequence(graph: ECGraph, steps: Sequence[Step]) -> List[Tuple[str, str]]:
    seq = []
    for i, s in enumerate(steps):
        if s.join[0] == END:
            continue  # a linear walk does not wrap around
        a = graph.adjacencies[s.join[1]]
        if a.kind == REFERENCE:
            continue
        seg = graph.segments[s.sid]
        exit_node = seg.head if s.strand == "+" else seg.tail
        nxt = steps[(i + 1) % len(steps)]
        nseg = graph.segments[nxt.sid]
        entry_node = nseg.tail if nxt.strand == "+" else nseg.head
        seq.append((exit_node, entry_node))
    return seq


def read_paths(path_tsv: str, graph: ECGraph) -> List[ReadPath]:
    """Long-read paths (amplicon-allele-analysis ``obs/path.tsv``) that touch only nodes of ``graph``.

    Only joins between consecutive path steps that are not reference adjacencies of the graph
    are kept, so paths are compared at junction resolution and tolerate a coarser or finer
    segmentation than the one the path was emitted on.
    """
    frame = pd.read_csv(path_tsv, sep="\t")
    ref_pairs = {frozenset((a.u, a.v)) for a in graph.adjacencies.values() if a.kind == REFERENCE}
    result = []
    for row in frame.to_dict("records"):
        tokens = [t.strip() for t in str(row["PATH"]).split(",") if t.strip()]
        nodes = []
        for t in tokens:
            eid, _, orient = t.rpartition("@")
            nodes.append(graph.resolve("{e}_{o}".format(e=eid, o=orient)))
        junctions, ok = [], True
        for i in range(1, len(nodes) - 1, 2):
            x, y = nodes[i], nodes[i + 1]
            if x is None or y is None:
                ok = False
                break
            if frozenset((x, y)) in ref_pairs:
                continue
            junctions.append((x, y))
        if ok and junctions:
            result.append(ReadPath(pid=str(row.get("PID", "")), support=float(row.get("#NOV", 0) or 0), junctions=junctions))
    return result


def _contains(seq: List[Tuple[str, str]], pattern: List[Tuple[str, str]], cyclic: bool = True) -> bool:
    n, m = len(seq), len(pattern)
    if m == 0 or n == 0 or m > n:
        return False
    starts = range(n) if cyclic else range(n - m + 1)
    for r in starts:
        if all(seq[(r + i) % n] == pattern[i] for i in range(m)):
            return True
    return False


def path_support(graph: ECGraph, steps: Sequence[Step], paths: Iterable[ReadPath]) -> Tuple[float, int, int]:
    """(support of matched multi-junction paths, matched paths, applicable paths)."""
    seq = _junction_sequence(graph, steps)
    rev = [(b, a) for a, b in reversed(seq)]
    cyclic = not is_linear(steps)
    nodes = {x for pair in seq for x in pair}
    score, matched, applicable = 0.0, 0, 0
    for p in paths:
        if not all(x in nodes and y in nodes for x, y in p.junctions):
            continue
        applicable += 1
        if _contains(seq, p.junctions, cyclic) or _contains(rev, p.junctions, cyclic):
            matched += 1
            if len(p.junctions) >= 2:
                score += p.support
    return score, matched, applicable
