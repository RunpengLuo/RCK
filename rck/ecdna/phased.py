"""Haplotype-specific ("phased") edges and nodes of an :class:`~rck.ecdna.subgraph_io.ECGraph`.

A phased node is ``(node_id, "A" | "B")``. A phased edge key is either

* ``("S", sid, hap)`` - one haplotype copy of a segment, joining ``(tail, hap)`` to ``(head, hap)``;
* ``("J", aid, phasing)`` - one labeling of an adjacency, joining ``(u, phasing[0])`` to
  ``(v, phasing[1])``, with ``phasing`` an RCK ``Phasing`` name such as ``"AB"``;
* ``("T", nid, hap)`` - a linear end: one copy of the molecule stops at ``(nid, hap)``. It joins
  that extremity to a virtual end node, and sits on the adjacency side of the extremity, so a
  segment copy that ends there needs no adjacency.

A circular species uses no ``T`` edges; a linear species uses exactly two (its two ends). For
walk extraction the two ends of one linear species are re-labelled ``T1`` / ``T2`` and joined by
a pseudo segment ``P``, which turns the open walk into a closed alternating circuit.

The flow MILP, the cycle decomposition and the walk extraction all speak this vocabulary.
"""

from __future__ import annotations

from typing import Dict, Iterable, List, Optional, Tuple

from rck.ecdna.subgraph_io import ECGraph

SEG = "S"
ADJ = "J"
END = "T"  # linear end of a molecule
END_A, END_B = "T1", "T2"  # the two ends of one linear walk, during walk extraction
PSEUDO = "P"  # pseudo segment joining END_A and END_B

END_NODE = ("__end__", "*")
END_NODE_A = ("__end_a__", "*")
END_NODE_B = ("__end_b__", "*")
PSEUDO_KEY = (PSEUDO, "*", "*")
ADJ_LIKE = (ADJ, END, END_A, END_B)
VIRTUAL_NODES = (END_NODE, END_NODE_A, END_NODE_B)

PhasedNode = Tuple[str, str]
EdgeKey = Tuple[str, str, str]


def seg_key(sid: str, hap) -> EdgeKey:
    return (SEG, sid, str(hap))


def adj_key(aid: str, phasing) -> EdgeKey:
    return (ADJ, aid, phasing if isinstance(phasing, str) else str(phasing))


def end_key(nid: str, hap) -> EdgeKey:
    return (END, nid, str(hap))


def is_adj_like(key: EdgeKey) -> bool:
    """Edges on the adjacency side of an extremity: adjacencies and linear ends."""
    return key[0] in ADJ_LIKE


def endpoints(graph: ECGraph, key: EdgeKey) -> Tuple[PhasedNode, PhasedNode]:
    kind, ident, label = key
    if kind == END:
        return (ident, label), END_NODE
    if kind == END_A:
        return (ident, label), END_NODE_A
    if kind == END_B:
        return (ident, label), END_NODE_B
    if kind == PSEUDO:
        return END_NODE_A, END_NODE_B
    if kind == SEG:
        s = graph.segments[ident]
        return (s.tail, label), (s.head, label)
    a = graph.adjacencies[ident]
    return (a.u, label[0]), (a.v, label[1])


def segment_key_at(graph: ECGraph, node: PhasedNode) -> Optional[EdgeKey]:
    if node in (END_NODE_A, END_NODE_B):
        return PSEUDO_KEY
    if node == END_NODE:
        return None
    sid = graph.nodes[node[0]].segment
    return None if sid is None else seg_key(sid, node[1])


def incidence(graph: ECGraph, keys: Iterable[EdgeKey]) -> Dict[PhasedNode, List[EdgeKey]]:
    """Adjacency-side edge keys (adjacencies and linear ends) incident to each phased node
    (a self-loop on one phased node is listed twice)."""
    result: Dict[PhasedNode, List[EdgeKey]] = {}
    for key in keys:
        if not is_adj_like(key):
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
    """Copy-number balance violations of an edge multiset (empty list = balanced).

    Linear ends count on the adjacency side of their extremity; virtual end nodes are not checked.
    """
    problems = []
    adj_load: Dict[PhasedNode, int] = {}
    for key, value in counts.items():
        if value <= 0 or not is_adj_like(key):
            continue
        x, y = endpoints(graph, key)
        adj_load[x] = adj_load.get(x, 0) + value
        adj_load[y] = adj_load.get(y, 0) + value
    nodes = set(adj_load)
    for key, value in counts.items():
        if value > 0 and key[0] == SEG:
            nodes.update(endpoints(graph, key))
    for node in nodes:
        if node in VIRTUAL_NODES:
            continue
        sk = segment_key_at(graph, node)
        seg_value = counts.get(sk, 0) if sk is not None else 0
        if seg_value != adj_load.get(node, 0):
            problems.append("{n}: segment {s} vs adjacencies {a}".format(n=node, s=seg_value, a=adj_load.get(node, 0)))
    return problems


def n_ends(counts: Dict[EdgeKey, int]) -> int:
    """Number of linear ends used by an edge multiset (0 for a circle, 2 for a linear walk)."""
    return sum(v for k, v in counts.items() if k[0] == END and v > 0)
