"""Turn a balanced phased edge multiset into ordered circular walks.

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
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd

from rck.ecdna.phased import ADJ, SEG, EdgeKey, PhasedNode, endpoints, incidence, segment_key_at
from rck.ecdna.subgraph_io import REFERENCE, ECGraph


@dataclass(frozen=True)
class Step:
    """One segment traversal and the adjacency taken right after it."""

    sid: str
    hap: str
    strand: str  # "+" tail -> head, "-" head -> tail
    join: EdgeKey  # ("J", aid, phasing) leading to the next step


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


def _to_steps(graph: ECGraph, items) -> List[Step]:
    steps = []
    for i in range(0, len(items), 2):
        _, sk, frm, _ = items[i]
        _, ak, _, _ = items[i + 1]
        tail = (graph.segments[sk[1]].tail, sk[2])
        steps.append(Step(sid=sk[1], hap=sk[2], strand="+" if frm == tail else "-", join=ak))
    return steps


def alternating_circuit(graph: ECGraph, counts: Dict[EdgeKey, int], start: Optional[EdgeKey] = None) -> List[Step]:
    """One alternating closed walk using every edge of ``counts`` exactly its multiplicity times."""
    state = _State(graph, counts)
    seg_keys = sorted(k for k in state.rem if k[0] == SEG)
    if not seg_keys:
        return []
    start_key = start or seg_keys[0]
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
    """Distinct alternating closed walks (up to rotation / reversal), at most ``limit`` of them."""
    state = _State(graph, counts)
    seg_keys = sorted(k for k in state.rem if k[0] == SEG)
    if not seg_keys:
        return []
    start = endpoints(graph, seg_keys[0])[0]
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
    return [Step(steps[i].sid, steps[i].hap, _flip(steps[i].strand), steps[(i - 1) % n].join) for i in range(n - 1, -1, -1)]


def canonical(steps: Sequence[Step]) -> Tuple:
    """Rotation- and reversal-invariant key of a circular walk."""
    if not steps:
        return ()
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
    return counts


# ---- description ---- #

def describe_step(graph: ECGraph, step: Step) -> str:
    seg = graph.segments[step.sid]
    return "{c}:{s}-{e}:{h}:{st}".format(c=seg.chrom, s=seg.start, e=seg.end, h=step.hap, st=step.strand)


def describe_join(graph: ECGraph, step: Step) -> str:
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


def _contains(cyclic: List[Tuple[str, str]], pattern: List[Tuple[str, str]]) -> bool:
    n, m = len(cyclic), len(pattern)
    if m == 0 or n == 0 or m > n:
        return False
    for r in range(n):
        if all(cyclic[(r + i) % n] == pattern[i] for i in range(m)):
            return True
    return False


def path_support(graph: ECGraph, steps: Sequence[Step], paths: Iterable[ReadPath]) -> Tuple[float, int, int]:
    """(support of matched multi-junction paths, matched paths, applicable paths)."""
    seq = _junction_sequence(graph, steps)
    rev = [(b, a) for a, b in reversed(seq)]
    nodes = {x for pair in seq for x in pair}
    score, matched, applicable = 0.0, 0, 0
    for p in paths:
        if not all(x in nodes and y in nodes for x, y in p.junctions):
            continue
        applicable += 1
        if _contains(seq, p.junctions) or _contains(rev, p.junctions):
            matched += 1
            if len(p.junctions) >= 2:
                score += p.support
    return score, matched, applicable
