"""Haplotype-specific ("phased") edges and nodes of an :class:`~rck.ecdna.subgraph_io.ECGraph`.

A phased node is ``(node_id, "A" | "B")``. A phased edge key is either

* ``("S", sid, hap)`` - one haplotype copy of a segment, joining ``(tail, hap)`` to ``(head, hap)``;
* ``("J", aid, phasing)`` - one labeling of an adjacency, joining ``(u, phasing[0])`` to
  ``(v, phasing[1])``, with ``phasing`` an RCK ``Phasing`` name such as ``"AB"``.

The flow MILP, the cycle decomposition and the walk extraction all speak this vocabulary.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple

from rck.ecdna.subgraph_io import ECGraph

SEG = "S"
ADJ = "J"

PhasedNode = Tuple[str, str]
EdgeKey = Tuple[str, str, str]


def seg_key(sid: str, hap) -> EdgeKey:
    return (SEG, sid, str(hap))


def adj_key(aid: str, phasing) -> EdgeKey:
    return (ADJ, aid, phasing if isinstance(phasing, str) else str(phasing))


def endpoints(graph: ECGraph, key: EdgeKey) -> Tuple[PhasedNode, PhasedNode]:
    kind, ident, label = key
    if kind == SEG:
        s = graph.segments[ident]
        return (s.tail, label), (s.head, label)
    a = graph.adjacencies[ident]
    return (a.u, label[0]), (a.v, label[1])


def segment_key_at(graph: ECGraph, node: PhasedNode) -> Optional[EdgeKey]:
    sid = graph.nodes[node[0]].segment
    return None if sid is None else seg_key(sid, node[1])


def incidence(graph: ECGraph, keys: Iterable[EdgeKey]) -> Dict[PhasedNode, List[EdgeKey]]:
    """Adjacency-edge keys incident to each phased node (a self-loop on one phased node is listed twice)."""
    result: Dict[PhasedNode, List[EdgeKey]] = {}
    for key in keys:
        if key[0] != ADJ:
            continue
        x, y = endpoints(graph, key)
        result.setdefault(x, []).append(key)
        result.setdefault(y, []).append(key)
    return result


def multiplicity_at(graph: ECGraph, key: EdgeKey, node: PhasedNode) -> int:
    """How many edge ends of ``key`` sit on ``node`` (2 for a self-loop on that phased node)."""
    x, y = endpoints(graph, key)
    return int(x == node) + int(y == node)


def is_balanced(graph: ECGraph, counts: Dict[EdgeKey, int]) -> List[str]:
    """Copy-number balance violations of an edge multiset (empty list = balanced)."""
    problems = []
    adj_load: Dict[PhasedNode, int] = {}
    for key, value in counts.items():
        if value <= 0 or key[0] != ADJ:
            continue
        x, y = endpoints(graph, key)
        adj_load[x] = adj_load.get(x, 0) + value
        adj_load[y] = adj_load.get(y, 0) + value
    nodes = set(adj_load)
    for key, value in counts.items():
        if value > 0 and key[0] == SEG:
            nodes.update(endpoints(graph, key))
    for node in nodes:
        sk = segment_key_at(graph, node)
        seg_value = counts.get(sk, 0) if sk is not None else 0
        if seg_value != adj_load.get(node, 0):
            problems.append("{n}: segment {s} vs adjacencies {a}".format(n=node, s=seg_value, a=adj_load.get(node, 0)))
    return problems
