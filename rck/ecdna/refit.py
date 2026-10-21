"""Refit species copy numbers against the observations: fractional weights, independent of the flow.

The decomposition (:mod:`rck.ecdna.decompose`) finds species *structures* by reproducing the
integer ecDNA-layer flow. That flow is fitted edge by edge to noisy observations, so it can be
uneven: 15 copies on most edges of a molecule and 28 on one segment. Reproducing it exactly
then needs leftover species, and an integer leftover can take copies that belong to the real
molecule (demo4 tumor2: 6 + 1 copies against a true 7.44).

Here every decomposition keeps its structures (the multiplicities ``m_ke``) and gets new
copy numbers ``w_k >= 0``, as real numbers (CoRAL also uses fractional cycle weights), that fit
the observations directly with the flow MILP's own data terms:

* each segment haplotype (or the segment total without allele data):
  ``seg_weight x |bg + sum_k w_k m_k,seg - observed| / sqrt(observed)``;
* each adjacency with read support:
  ``alpha x |bg + sum_k w_k m_k,adj - observed| / sqrt(observed)``;

with the background layer held at the flow solution, plus ``weight_penalty`` per copy as a
tie-break. It is a small linear program (solved with scipy's HiGHS, no Gurobi licence needed).
A species refitted below ``min_copies`` is not supported by the data: the lightest such species
is dropped and the fit repeated, so leftovers that only absorbed flow unevenness disappear.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from rck.ecdna.phased import EdgeKey, adj_key, seg_key
from rck.ecdna.subgraph_io import HAPLOTYPES, REFERENCE, ECGraph

log = logging.getLogger(__name__)


@dataclass
class RefitParams:
    enabled: bool = True
    min_copies: float = 0.5  # species refitted below this many copies are dropped
    weight_penalty: float = 1e-3  # per copy: among equal fits, prefer fewer copies
    digits: int = 2  # reported precision


@dataclass
class Row:
    cost: float
    target: float  # observed copy number minus the background
    keys: Tuple[EdgeKey, ...]  # phased edges whose ecDNA copies add up to this observation


@dataclass
class RefitResult:
    weights: List[float]  # aligned with ``kept``
    kept: List[int]  # indices into the species passed in
    dropped: List[Tuple[int, float]]  # (index, weight when dropped)
    objective: float  # data misfit of the final fit
    objective_int: float  # data misfit with the integer weights
    status: str


def observation_rows(graph: ECGraph, obs, flow, cparams) -> List[Row]:
    """The flow MILP's data terms with the background fixed at ``flow`` (a FlowSolution)."""
    rows: List[Row] = []
    for sid in sorted(graph.segments):
        s, so = graph.segments[sid], obs.segments[sid]
        w = max(min(s.length / float(cparams.seg_length_cap), 1.0), cparams.min_seg_weight)
        norm = 1.0 / math.sqrt(max(so.total, 1.0))
        if so.a is not None and so.b is not None:
            for h, value in zip(HAPLOTYPES, (so.a, so.b)):
                rows.append(Row(w * norm, value - flow.seg_bg.get((sid, h), 0), (seg_key(sid, h),)))
        else:
            bg = sum(flow.seg_bg.get((sid, h), 0) for h in HAPLOTYPES)
            rows.append(Row(w * norm, so.total - bg, tuple(seg_key(sid, h) for h in HAPLOTYPES)))
    for aid in sorted(graph.adjacencies):
        a, value = graph.adjacencies[aid], obs.adjacencies.get(aid)
        if value is None:
            continue
        alpha = cparams.alpha_ref if a.kind == REFERENCE else cparams.alpha_novel
        bg = sum(flow.adj_bg.get((aid, ph), 0) for ph in a.allowed)
        rows.append(Row(alpha / math.sqrt(max(value, 1.0)), value - bg, tuple(adj_key(aid, ph) for ph in a.allowed)))
    return rows


def _coefficients(rows: Sequence[Row], counts: Sequence[Dict[EdgeKey, int]]) -> List[List[float]]:
    return [[float(sum(c.get(key, 0) for key in r.keys)) for c in counts] for r in rows]


def misfit(rows: Sequence[Row], counts: Sequence[Dict[EdgeKey, int]], weights: Sequence[float]) -> float:
    coef = _coefficients(rows, counts)
    return float(sum(r.cost * abs(sum(a * w for a, w in zip(coef[i], weights)) - r.target) for i, r in enumerate(rows)))


def _solve_lp(rows: Sequence[Row], counts: Sequence[Dict[EdgeKey, int]], penalty: float) -> Tuple[List[float], str]:
    import numpy as np
    from scipy.optimize import linprog

    K, N = len(counts), len(rows)
    coef = np.array(_coefficients(rows, counts), dtype=float).reshape(N, K)
    target = np.array([r.target for r in rows], dtype=float)
    # variables: w (K), d (N); minimise penalty * sum(w) + sum(cost * d), d >= |coef w - target|
    c = np.concatenate([np.full(K, penalty), np.array([r.cost for r in rows], dtype=float)])
    eye = np.eye(N)
    a_ub = np.vstack([np.hstack([coef, -eye]), np.hstack([-coef, -eye])])
    b_ub = np.concatenate([target, -target])
    res = linprog(c, A_ub=a_ub, b_ub=b_ub, bounds=[(0, None)] * (K + N), method="highs")
    if not res.success:
        return [], str(res.message)
    return [float(v) for v in res.x[:K]], "optimal"


def refit_weights(rows: Sequence[Row], counts: Sequence[Dict[EdgeKey, int]], weights_int: Sequence[float],
                  params: Optional[RefitParams] = None) -> RefitResult:
    params = params or RefitParams()
    objective_int = misfit(rows, counts, weights_int)
    kept = list(range(len(counts)))
    dropped: List[Tuple[int, float]] = []
    while True:
        weights, status = _solve_lp(rows, [counts[i] for i in kept], params.weight_penalty)
        if not weights:  # keep the integer answer if the LP fails
            log.warning("weight refit failed (%s); keeping integer weights", status)
            return RefitResult([float(weights_int[i]) for i in kept], kept, dropped, objective_int, objective_int, status)
        j = min(range(len(kept)), key=lambda i: weights[i])
        if len(kept) > 1 and weights[j] < params.min_copies:
            dropped.append((kept[j], weights[j]))
            kept.pop(j)
            continue
        objective = misfit(rows, [counts[i] for i in kept], weights)
        return RefitResult(weights, kept, dropped, objective, objective_int, status)


def refit_decompositions(decompositions, graph: ECGraph, obs, flow, cparams, dparams, params: Optional[RefitParams] = None):
    """Refit every decomposition in place: fractional weights, unsupported species dropped, calls updated."""
    from rck.ecdna.decompose import _flag_low_weight, classify

    params = params or RefitParams()
    if not params.enabled:
        return decompositions
    rows = observation_rows(graph, obs, flow, cparams)
    for d in decompositions:
        if not d.species:
            continue
        counts = [s.counts for s in d.species]
        result = refit_weights(rows, counts, [s.weight for s in d.species], params)
        d.data_objective, d.data_objective_int = round(result.objective, 6), round(result.objective_int, 6)
        d.pruned = [{"species_index": i + 1, "weight_int": d.species[i].weight_int, "refit_weight": round(w, params.digits),
                     "topology": d.species[i].topology, "length_bp": d.species[i].length_bp} for i, w in result.dropped]
        kept = []
        for i, w in zip(result.kept, result.weights):
            s = d.species[i]
            s.weight = round(w, params.digits)
            kept.append(s)
        kept.sort(key=lambda s: (-s.weight, -sum(s.counts.values())))
        _flag_low_weight(kept, dparams)
        d.species = kept
        d.call, d.circular_fraction = classify(kept, dparams)
    return decompositions
