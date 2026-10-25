"""Decompose the ecDNA flow into the fewest species (circular or linear), with alternatives.

Given the ecDNA-layer flow ``F_e`` on phased edges (from :mod:`rck.ecdna.ilp_circular`), find
up to ``K`` species slots. Slot ``k`` has an integer copy weight ``w_k`` and a small integer
multiplicity ``m_ke`` (0..``max_multiplicity``) on every edge, such that

* ``sum_k w_k * m_ke`` reproduces ``F_e`` (deviations are penalized, not forbidden);
* each slot is balanced with alternation at every phased node
  (segment multiplicity == sum of incident adjacency multiplicities, self-loops twice);
* each slot's support is connected (lazy cuts in a Gurobi callback), so with balance it is
  one alternating walk;
* each slot uses either no linear ends (a circle) or exactly two (a linear walk); linear ends
  share one virtual node, so connectivity also keeps a linear slot to a single walk.

The product ``w_k * m_ke`` is linearized exactly through the binary expansion of ``m_ke``
(McCormick on binary x integer). Objective: ``sum_e omega_e |residual_e| + species_penalty *
#species``. Slots are ordered by weight to break symmetry.

**Walk complexity (``step_penalty``).** Nested species can be explained equally well by
separate circles or by walks that traverse the shared sequence again (an outer circle plus a
walk through outer and nested together, or one walk round the outer cycle twice and the nested
one once). These reproduce the flow, and the copy numbers, exactly, so only a prior can choose.
A step is one segment traversal of one species (the walk length in segments, not weighted by
copies). Splitting a species into pieces keeps the total step count about the same, while
merging nested species adds the re-traversed steps, so ``step_penalty`` favours separate simple
walks without pushing towards splits. It must stay below ``species_penalty`` per step that a
true internal duplication needs.

The penalty is applied in the ranking, to significant species only, not in the MILP: in the
MILP it would also charge the light leftover species that absorb flow unevenness, whose walks
are often long, and the solver then reshaped real species to shorten them (demo3, demo7); a
weight threshold inside the MILP is gamed by splitting a merged walk into light copies. With
``step_penalty`` > 0 the solution pool is widened to every solution within ``ambiguity_window``
of the best objective (and never less than ``pool_gap``), so that the alternatives the prior
chooses between, which can need more species, reach the ranking. ``multiplicity_penalty``
remains in the MILP as a much smaller tie-break on all edge traversals.

**Ambiguity.** After ranking, the top decomposition is checked for alternatives that the data
cannot tell apart (:meth:`CycleDecomposer.same_fit_alternatives`). The ecDNA flow covered by
each edge is fixed to the top decomposition's, two species may not share the same support (which
excludes splitting one walk into copies), and alternatives within ``ambiguity_window`` of the
objective are enumerated. Any alternative whose significant species structures differ means the
choice was made by the priors (species count, steps, tie-break) or by the reads, not by the copy
numbers; the gap says by how much.

Alternatives come from Gurobi's solution pool (``PoolSearchMode=2``) or, with
``mode='nogood'``, from re-solving with a no-good cut on the slot supports. Pool solutions
are re-checked for balance and connectivity and deduplicated by their species multisets.

Each decomposition gets an amplicon call from its significant (not ``low_weight``) species:
``ecDNA`` when circular species carry at least 80% of the copy-number x length, ``linear`` at
most 20%, ``mixed`` in between, ``none`` without significant species.

**Long reads (CoRAL-style).** Every long-read path that spans at least ``min_path_junctions``
junctions (and has ``min_path_support`` reads) whose junctions all carry ecDNA flow becomes a
read constraint. For read ``p`` and slot ``k`` a binary ``z_kp`` equals 1 exactly when slot
``k`` uses every junction of ``p`` (any labeling), and ``x_p = OR_k z_kp`` marks the read as
explained. Unexplained reads cost ``path_weight`` x their share of the read support (so at most
``path_weight`` in total), and ``path_fraction`` optionally requires that share of the support
to be explained (CoRAL uses 0.9). ``z`` and ``x`` are fully determined by the slot supports, so
they do not multiply the solution pool. Containing every junction is necessary but not
sufficient for a read to fit a walk; after solving, each read is checked against the actual
orderings (:func:`rck.ecdna.walk.path_matches`) and the explained share is reported.

Species weights here are integers (the structure search). :mod:`rck.ecdna.refit` then refits
them as fractional copy numbers against the observations.
"""

from __future__ import annotations

import logging
import math
from functools import reduce
from dataclasses import dataclass, field, replace
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

from rck.ecdna.gurobi_env import close_model, default_threads
from rck.ecdna.phased import (ADJ, END, END_NODE, SEG, EdgeKey, PhasedNode, endpoints, is_adj_like, is_balanced, n_ends,
                              segment_key_at)
from rck.ecdna.subgraph_io import IMPUTED, NOVEL, REFERENCE, ECGraph
from rck.ecdna.walk import ReadPath, Step, alternating_circuit, enumerate_circuits, path_matches, path_support, walk_length

log = logging.getLogger(__name__)

GRB_INFINITY = 1e100


@dataclass
class DecompositionParams:
    max_species: int = 4
    max_multiplicity: int = 3
    species_penalty: float = 1.0
    multiplicity_penalty: float = 1e-4  # tie-break: prefer the walk with fewer edge traversals
    step_penalty: float = 0.0  # cost per segment traversal of a species (walk complexity); 0 = tie-break only
    ambiguity_window: float = 2.0  # objective window for same-fit alternatives to the top decomposition; 0 = off
    ambiguity_pool: int = 50  # same-fit alternatives enumerated (Gurobi pool size)
    min_species_weight: Optional[int] = None  # species below this weight are flagged low_weight
    min_species_weight_frac: float = 0.10  # default threshold: this fraction of the heaviest species
    min_ecdna_copies: int = 2  # circular species lighter than this do not make an amplicon ecDNA
    ecdna_call_fraction: float = 0.8  # circular share of copy-number x length for an "ecDNA" call
    adj_residual_weight: float = 1.0
    seg_length_cap: int = 50_000
    min_seg_weight: float = 0.05
    num_alternatives: int = 10
    pool_gap: float = 0.10
    mode: str = "pool"  # pool | nogood
    max_orderings: int = 20
    path_weight: float = 0.5  # cost of leaving all multi-junction read support unexplained (scaled by share)
    path_fraction: float = 0.0  # required explained share of that support (CoRAL: 0.9); 0 = off
    min_path_support: float = 1.0  # reads behind a path before it constrains the decomposition
    min_path_junctions: int = 2  # single-junction paths are already counted by the junction reads
    max_expansions: int = 200_000
    time_limit: Optional[float] = None
    threads: Optional[int] = None
    mip_gap: float = 1e-4
    log_file: Optional[str] = None
    log_to_console: bool = False


@dataclass
class Ordering:
    steps: List[Step]
    path_score: float = 0.0
    paths_matched: int = 0
    paths_applicable: int = 0


@dataclass
class Species:
    weight: float  # copies per cell: integer from the structure search, fractional after refit
    counts: Dict[EdgeKey, int]
    orderings: List[Ordering] = field(default_factory=list)
    length_bp: int = 0
    imputed: List[str] = field(default_factory=list)
    junctions: List[str] = field(default_factory=list)
    low_weight: bool = False  # below the minimum species weight: likely flow residue, not a molecule
    weight_int: Optional[int] = None  # integer weight from the structure search

    @property
    def topology(self) -> str:
        return "linear" if n_ends(self.counts) else "circular"

    @property
    def status(self) -> str:
        return "low_weight" if self.low_weight else "species"

    @property
    def signature(self) -> Tuple[float, FrozenSet]:
        return (self.weight if self.weight_int is None else self.weight_int), frozenset(self.counts.items())


@dataclass
class Decomposition:
    rank: int
    objective: float
    residual: float
    species: List[Species]
    call: str = "none"
    circular_fraction: Optional[float] = None
    reads_total: float = 0.0  # support of the read constraints (multi-junction reads on the flow)
    reads_explained: float = 0.0  # share of it found, in order, in some species' top ordering
    data_objective: Optional[float] = None  # observation misfit after the weight refit
    data_objective_int: Optional[float] = None  # the same with the integer weights
    pruned: List[Dict[str, object]] = field(default_factory=list)  # species the refit dropped
    steps: int = 0  # segment traversals of the significant species
    ambiguity: Optional[Dict[str, object]] = None  # same-fit alternatives (top decomposition only)

    @property
    def signature(self) -> FrozenSet:
        return frozenset(s.signature for s in self.species)


@dataclass
class ReadConstraint:
    path: ReadPath  # representative read (orientation as first seen)
    support: float  # reads behind this junction sequence (merged over identical paths)
    groups: List[Tuple[EdgeKey, ...]]  # per distinct junction: the phased adjacency edges that can carry it
    pids: List[str] = field(default_factory=list)


def read_constraints(graph: ECGraph, edges: Sequence[EdgeKey], paths: Optional[Sequence[ReadPath]],
                     min_support: float = 1.0, min_junctions: int = 2) -> Tuple[List[ReadConstraint], int]:
    """Read constraints for the decomposition, and the number of reads left out because a junction has no flow."""
    by_pair: Dict[FrozenSet[str], List[EdgeKey]] = {}
    for e in edges:
        if e[0] == ADJ and graph.adjacencies[e[1]].kind != REFERENCE:
            a = graph.adjacencies[e[1]]
            by_pair.setdefault(frozenset((a.u, a.v)), []).append(e)
    merged: Dict[Tuple, ReadConstraint] = {}
    skipped = 0
    for p in paths or []:
        if len(p.junctions) < min_junctions or p.support < min_support:
            continue
        groups: List[Tuple[EdgeKey, ...]] = []
        for x, y in p.junctions:
            group = tuple(sorted(by_pair.get(frozenset((x, y)), [])))
            if not group:
                groups = []
                break
            if group not in groups:
                groups.append(group)
        if not groups:
            skipped += 1
            continue
        forward = tuple(p.junctions)
        key = min(forward, tuple((b, a) for a, b in reversed(p.junctions)))
        if key in merged:
            merged[key].support += p.support
            merged[key].pids.append(p.pid)
        else:
            merged[key] = ReadConstraint(path=p, support=p.support, groups=groups, pids=[p.pid])
    return sorted(merged.values(), key=lambda r: (-r.support, r.pids[0])), skipped


def _components(graph: ECGraph, keys: Sequence[EdgeKey]) -> List[List[EdgeKey]]:
    parent: Dict[PhasedNode, PhasedNode] = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for k in keys:
        x, y = endpoints(graph, k)
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
    groups: Dict[PhasedNode, List[EdgeKey]] = {}
    for k in keys:
        groups.setdefault(find(endpoints(graph, k)[0]), []).append(k)
    return list(groups.values())


class CycleDecomposer:
    def __init__(self, graph: ECGraph, flow: Dict[EdgeKey, int], params: Optional[DecompositionParams] = None,
                 paths: Optional[Sequence[ReadPath]] = None):
        self.graph = graph
        self.flow = {k: int(v) for k, v in flow.items() if v > 0}
        self.params = params or DecompositionParams()
        self.edges = sorted(self.flow)
        self.index = {e: i for i, e in enumerate(self.edges)}
        self.model = None
        self.discarded = 0
        p = self.params
        self.reads, self.reads_skipped = read_constraints(graph, self.edges, paths, p.min_path_support, p.min_path_junctions)
        if not (p.path_weight > 0 or p.path_fraction > 0):
            self.reads = []

    def _omega(self, key: EdgeKey) -> float:
        p = self.params
        if key[0] == SEG:
            s = self.graph.segments[key[1]]
            return max(min(s.length / float(p.seg_length_cap), 1.0), p.min_seg_weight)
        return p.adj_residual_weight

    def build(self, distinct_supports: bool = False):
        """Build the MILP. With ``distinct_supports`` no two species may use the same set of edges
        (used by :meth:`same_fit_alternatives` to exclude one walk split into copies)."""
        import gurobipy as gp
        from gurobipy import GRB

        p, g, E = self.params, self.graph, self.edges
        m = gp.Model("rck-ecdna-decompose")
        self.model = m
        m.Params.OutputFlag = 1 if (p.log_to_console or p.log_file) else 0
        m.Params.LogToConsole = 1 if p.log_to_console else 0
        if p.log_file:
            m.Params.LogFile = p.log_file
        if p.time_limit:
            m.Params.TimeLimit = p.time_limit
        m.Params.Threads = default_threads(p.threads)
        m.Params.MIPGap = p.mip_gap
        m.Params.LazyConstraints = 1

        W = max(self.flow.values())
        K = max(1, p.max_species)
        nbits = max(1, int(math.ceil(math.log2(p.max_multiplicity + 1))))
        self.K, self.W, self.nbits = K, W, nbits
        self.u = {k: m.addVar(vtype=GRB.BINARY, name="u[{k}]".format(k=k)) for k in range(K)}
        self.w = {k: m.addVar(lb=0, ub=W, vtype=GRB.INTEGER, name="w[{k}]".format(k=k)) for k in range(K)}
        self.beta, self.zeta, self.y, self.mult, self.linear = {}, {}, {}, {}, {}
        for k in range(K):
            m.addConstr(self.w[k] <= W * self.u[k])
            m.addConstr(self.w[k] >= self.u[k])
            for e in E:
                bits = []
                for b in range(nbits):
                    beta = m.addVar(vtype=GRB.BINARY, name="beta[{k},{i},{b}]".format(k=k, i=self.index[e], b=b))
                    zeta = m.addVar(lb=0, ub=W, name="zeta[{k},{i},{b}]".format(k=k, i=self.index[e], b=b))
                    m.addConstr(zeta <= W * beta)
                    m.addConstr(zeta <= self.w[k])
                    m.addConstr(zeta >= self.w[k] - W * (1 - beta))
                    self.beta[k, e, b], self.zeta[k, e, b] = beta, zeta
                    bits.append(beta)
                mult = gp.quicksum((2 ** b) * bits[b] for b in range(nbits))
                self.mult[k, e] = mult
                y = m.addVar(vtype=GRB.BINARY, name="y[{k},{i}]".format(k=k, i=self.index[e]))
                self.y[k, e] = y
                m.addConstr(mult <= p.max_multiplicity)
                m.addConstr(mult <= p.max_multiplicity * y)
                m.addConstr(y <= mult)
                m.addConstr(y <= self.u[k])
            m.addConstr(gp.quicksum(self.y[k, e] for e in E) >= 2 * self.u[k])
            # a slot is a circle (no ends) or one linear walk (exactly two ends)
            ends = [e for e in E if e[0] == END]
            if ends:
                lin = m.addVar(vtype=GRB.BINARY, name="linear[{k}]".format(k=k))
                self.linear[k] = lin
                m.addConstr(gp.quicksum(self.mult[k, e] for e in ends) == 2 * lin)
                m.addConstr(lin <= self.u[k])
            # balance with alternation, per phased node (linear ends sit on the adjacency side)
            nodes: Dict[PhasedNode, List[EdgeKey]] = {}
            for e in E:
                if is_adj_like(e):
                    x, z = endpoints(g, e)
                    nodes.setdefault(x, []).append(e)
                    nodes.setdefault(z, []).append(e)
                else:
                    for x in endpoints(g, e):
                        nodes.setdefault(x, [])
            nodes.pop(END_NODE, None)  # the virtual end node has no segment to balance against
            for j, (node, adjs) in enumerate(sorted(nodes.items())):
                sk = segment_key_at(g, node)
                seg_expr = self.mult[k, sk] if (sk is not None and sk in self.flow) else 0
                m.addConstr(seg_expr == gp.quicksum(self.mult[k, a] for a in adjs), name="bal[{k},{j}]".format(k=k, j=j))
        for k in range(K - 1):
            m.addConstr(self.w[k] >= self.w[k + 1])
            m.addConstr(self.u[k] >= self.u[k + 1])

        if distinct_supports:
            for k in range(K):
                for l in range(k + 1, K):
                    diffs = []
                    for e in E:
                        d = m.addVar(lb=0, ub=1, name="xor[{k},{l},{i}]".format(k=k, l=l, i=self.index[e]))
                        a, b = self.y[k, e], self.y[l, e]
                        m.addConstr(d >= a - b)
                        m.addConstr(d >= b - a)
                        m.addConstr(d <= a + b)
                        m.addConstr(d <= 2 - a - b)
                        diffs.append(d)
                    m.addConstr(gp.quicksum(diffs) >= self.u[l], name="distinct[{k},{l}]".format(k=k, l=l))

        self.res, self.cover = {}, {}
        terms = []
        for e in E:
            cover = gp.quicksum((2 ** b) * self.zeta[k, e, b] for k in range(K) for b in range(nbits))
            self.cover[e] = cover
            r = m.addVar(lb=0, name="r[{i}]".format(i=self.index[e]))
            m.addConstr(r >= cover - self.flow[e])
            m.addConstr(r >= self.flow[e] - cover)
            self.res[e] = r
            terms.append(self._omega(e) * r)
        terms.append(p.species_penalty * gp.quicksum(self.u.values()))
        terms.extend(self._read_terms(gp, GRB))
        if p.multiplicity_penalty:
            terms.append(p.multiplicity_penalty * gp.quicksum(self.mult[k, e] for k in range(K) for e in E))
        m.setObjective(gp.quicksum(terms), GRB.MINIMIZE)
        m.update()
        return m

    def _read_terms(self, gp, GRB) -> list:
        """Read constraints (see module docstring); returns the objective terms."""
        p, m = self.params, self.model
        self.o, self.z, self.x = {}, {}, {}
        if not self.reads:
            return []
        total = sum(r.support for r in self.reads)
        for k in range(self.K):
            for ri, r in enumerate(self.reads):
                ors = []
                for group in r.groups:
                    if (k, group) not in self.o:
                        if len(group) == 1:
                            self.o[k, group] = self.y[k, group[0]]
                        else:  # several labelings of the junction: OR of their supports
                            o = m.addVar(vtype=GRB.BINARY, name="o[{k},{r},{n}]".format(k=k, r=ri, n=len(self.o)))
                            for e in group:
                                m.addConstr(o >= self.y[k, e])
                            m.addConstr(o <= gp.quicksum(self.y[k, e] for e in group))
                            self.o[k, group] = o
                    ors.append(self.o[k, group])
                z = m.addVar(vtype=GRB.BINARY, name="z[{k},{r}]".format(k=k, r=ri))
                for o in ors:
                    m.addConstr(z <= o)
                m.addConstr(z >= gp.quicksum(ors) - (len(ors) - 1))
                self.z[k, ri] = z
        for ri in range(len(self.reads)):
            x = m.addVar(vtype=GRB.BINARY, name="x[{r}]".format(r=ri))
            zs = [self.z[k, ri] for k in range(self.K)]
            for z in zs:
                m.addConstr(x >= z)
            m.addConstr(x <= gp.quicksum(zs))
            self.x[ri] = x
        explained = gp.quicksum(r.support * self.x[ri] for ri, r in enumerate(self.reads))
        if p.path_fraction > 0:
            m.addConstr(explained >= p.path_fraction * total, name="path_fraction")
        if p.path_weight > 0:
            # cost of the unexplained share; kept non-negative so the pool's relative gap is unaffected
            return [p.path_weight * (1 - explained / total)]
        return []

    def _callback(self):
        import gurobipy as gp
        from gurobipy import GRB

        E, g = self.edges, self.graph

        def cb(model, where):
            if where != GRB.Callback.MIPSOL:
                return
            for k in range(self.K):
                values = model.cbGetSolution([self.y[k, e] for e in E])
                used = [e for e, v in zip(E, values) if v > 0.5]
                if not used:
                    continue
                comps = _components(g, used)
                if len(comps) <= 1:
                    continue
                for i in range(len(comps)):
                    side = {x for e in comps[i] for x in endpoints(g, e)}
                    cut = [f for f in E if sum(x in side for x in endpoints(g, f)) == 1]
                    for j in range(len(comps)):
                        if i == j:
                            continue
                        model.cbLazy(gp.quicksum(self.y[k, f] for f in cut) >= self.y[k, comps[i][0]] + self.y[k, comps[j][0]] - 1)

        return cb

    def _read(self, getter) -> Tuple[List[Species], bool]:
        species, valid = [], True
        for k in range(self.K):
            if getter(self.u[k]) < 0.5:
                continue
            weight = int(round(getter(self.w[k])))
            counts = {}
            for e in self.edges:
                value = sum((2 ** b) * getter(self.beta[k, e, b]) for b in range(self.nbits))
                value = int(round(value))
                if value > 0:
                    counts[e] = value
            if not counts or weight <= 0:
                continue
            # a walk that goes round the same circle g times is that circle at g times the weight
            g = reduce(math.gcd, counts.values())
            if g > 1 and not n_ends(counts):
                counts = {e: v // g for e, v in counts.items()}
                weight *= g
            if is_balanced(self.graph, counts) or len(_components(self.graph, list(counts))) > 1 or n_ends(counts) not in (0, 2):
                valid = False
            species.append(Species(weight=weight, counts=counts, weight_int=weight))
        return species, valid

    def _residual(self, species: List[Species]) -> float:
        total = 0.0
        for e in self.edges:
            cover = sum(s.weight * s.counts.get(e, 0) for s in species)
            total += self._omega(e) * abs(cover - self.flow[e])
        return total

    def solve(self) -> List[Tuple[float, List[Species]]]:
        """Raw (objective, species) solutions, best first, deduplicated and validated."""
        import gurobipy as gp

        if not self.edges:
            return [(0.0, [])]
        if self.model is None:
            self.build()
        m, p = self.model, self.params
        cb = self._callback()
        solutions: List[Tuple[float, List[Species]]] = []
        seen = set()

        def collect(objective, species):
            sig = frozenset(s.signature for s in species)
            if sig in seen:
                return
            seen.add(sig)
            solutions.append((objective, species))

        if p.mode == "pool":
            m.Params.PoolSearchMode = 2
            m.Params.PoolSolutions = max(20, 4 * p.num_alternatives)
            m.Params.PoolGap = p.pool_gap
            if p.step_penalty > 0:
                # the ranking prior may prefer solutions with more species: keep an absolute window
                m.Params.PoolGap = GRB_INFINITY
                m.Params.PoolGapAbs = p.ambiguity_window if p.ambiguity_window > 0 else 2.0 * p.species_penalty
                m.Params.PoolSolutions = max(100, 4 * p.num_alternatives)
            m.optimize(cb)
            if m.SolCount == 0:
                status = m.Status
                self.close()
                raise RuntimeError("decomposition MILP returned no solution (status {s})".format(s=status))
            for i in range(m.SolCount):
                m.Params.SolutionNumber = i
                species, valid = self._read(lambda v: v.Xn)
                if not valid:
                    self.discarded += 1
                    continue
                # keep every valid pool solution: decompose() re-ranks them before truncating
                collect(float(m.PoolObjVal), species)
        elif p.mode == "nogood":
            best = None
            for _ in range(p.num_alternatives):
                m.optimize(cb)
                if m.SolCount == 0:
                    break
                objective = float(m.ObjVal)
                window = abs(best) * p.pool_gap if best is not None else 0.0
                if p.step_penalty > 0:
                    window = max(window, p.ambiguity_window if p.ambiguity_window > 0 else 2.0 * p.species_penalty)
                if best is None:
                    best = objective
                elif objective > best + window + 1e-6:
                    break
                species, valid = self._read(lambda v: v.X)
                if valid:
                    collect(objective, species)
                else:
                    self.discarded += 1
                ones = [self.y[key] for key, var in self.y.items() if var.X > 0.5]
                zeros = [self.y[key] for key, var in self.y.items() if var.X <= 0.5]
                m.addConstr(gp.quicksum(1 - v for v in ones) + gp.quicksum(zeros) >= 1)
        else:
            raise ValueError("unknown decomposition mode {m}".format(m=p.mode))
        if self.discarded:
            log.warning("discarded %d pool solution(s) that violated balance or connectivity", self.discarded)
        self.close()
        solutions.sort(key=lambda t: (round(t[0], 9), len(t[1])))
        if p.step_penalty > 0 and len(solutions) > 1:
            # the widened pool is for the ranking prior; keep what it can prefer: within pool_gap of
            # the best, or cheaper in structure cost than the best by more than the objective lost
            best, best_species = solutions[0]
            best_cost = _step_cost(best_species, p)
            solutions = [t for t in solutions if t[0] <= best + abs(best) * p.pool_gap + 1e-6
                         or _step_cost(t[1], p) + t[0] - best < best_cost]
        return solutions

    def same_fit_alternatives(self, species: Sequence[Species]) -> List[Tuple[float, List[Species]]]:
        """Decompositions that cover every edge exactly as ``species`` do, with no two species on
        the same support, within ``ambiguity_window`` of the best such objective (raw, valid,
        deduplicated; ``species`` itself may be among them)."""
        p = self.params
        if not self.edges:
            return []
        target = {e: sum(int(s.weight_int if s.weight_int is not None else s.weight) * s.counts.get(e, 0) for s in species)
                  for e in self.edges}
        m = self.build(distinct_supports=True)
        for e in self.edges:
            m.addConstr(self.cover[e] == target[e], name="same_fit[{i}]".format(i=self.index[e]))
        m.Params.PoolSearchMode = 2
        m.Params.PoolSolutions = max(1, p.ambiguity_pool)
        m.Params.PoolGapAbs = p.ambiguity_window
        m.optimize(self._callback())
        result, seen = [], set()
        for i in range(m.SolCount):
            m.Params.SolutionNumber = i
            found, valid = self._read(lambda v: v.Xn)
            sig = frozenset(s.signature for s in found)
            if valid and sig not in seen:
                seen.add(sig)
                result.append((float(m.PoolObjVal), found))
        self.close()
        result.sort(key=lambda t: t[0])
        return result

    def close(self) -> None:
        """Free the Gurobi model and close its log file (see OptModelCircular.close)."""
        close_model(self.model)
        self.model = None


def _annotate(graph: ECGraph, species: Species, params: DecompositionParams, paths: Optional[List[ReadPath]]) -> Species:
    orderings = enumerate_circuits(graph, species.counts, limit=params.max_orderings, max_expansions=params.max_expansions)
    if not orderings:
        orderings = [alternating_circuit(graph, species.counts)]
    species.orderings = []
    for steps in orderings:
        score, matched, applicable = path_support(graph, steps, paths or [])
        species.orderings.append(Ordering(steps=steps, path_score=score, paths_matched=matched, paths_applicable=applicable))
    species.orderings.sort(key=lambda o: (-o.path_score, -o.paths_matched))
    species.length_bp = walk_length(graph, species.orderings[0].steps)
    adj_ids = sorted({k[1] for k in species.counts if k[0] == ADJ})
    species.imputed = [a for a in adj_ids if graph.adjacencies[a].kind == IMPUTED]
    species.junctions = [a for a in adj_ids if graph.adjacencies[a].kind == NOVEL]
    return species


def classify(species: Sequence[Species], params: DecompositionParams) -> Tuple[str, Optional[float]]:
    """(call, circular fraction) of a decomposition from its significant species."""
    significant = [s for s in species if not s.low_weight]
    circular = sum(s.weight * s.length_bp for s in significant if s.topology == "circular" and s.weight >= params.min_ecdna_copies)
    linear = sum(s.weight * s.length_bp for s in significant if s.topology == "linear")
    total = circular + linear
    if total <= 0:
        return "none", None
    fraction = circular / float(total)
    if fraction >= params.ecdna_call_fraction:
        return "ecDNA", fraction
    if fraction <= 1.0 - params.ecdna_call_fraction:
        return "linear", fraction
    return "mixed", fraction


def _flag_low_weight(species: List[Species], params: DecompositionParams) -> None:
    if not species:
        return
    threshold = params.min_species_weight
    if threshold is None:  # for integer weights this is the same as ceil(frac x heaviest)
        threshold = max(1.0, params.min_species_weight_frac * max(s.weight for s in species))
    for s in species:
        s.low_weight = len(species) > 1 and s.weight < threshold


def _steps(species: Sequence[Species]) -> int:
    """Segment traversals of the species (each walk's length in segments)."""
    return sum(v for s in species for k, v in s.counts.items() if k[0] == SEG)


def _structures(species: Sequence[Species]) -> FrozenSet:
    """The significant species' walks, ignoring weights."""
    return frozenset(frozenset(s.counts.items()) for s in species if not s.low_weight)


def _step_cost(species: Sequence[Species], params: DecompositionParams) -> float:
    """step_penalty x steps of the significant species of a raw pool solution (low_weight decided on a copy)."""
    flagged = [replace(s) for s in species]
    _flag_low_weight(flagged, params)
    return params.step_penalty * _steps([s for s in flagged if not s.low_weight])


def structure_cost(species: Sequence[Species], params: DecompositionParams) -> float:
    """species_penalty x significant species + step_penalty x their steps (the ranking prior)."""
    significant = [s for s in species if not s.low_weight]
    return params.species_penalty * len(significant) + params.step_penalty * _steps(significant)


def decompose(graph: ECGraph, flow: Dict[EdgeKey, int], params: Optional[DecompositionParams] = None,
              paths: Optional[List[ReadPath]] = None) -> List[Decomposition]:
    """Ranked decompositions of ``flow`` into circular or linear species, each with ranked orderings.

    Decompositions are ranked by (residual, structure cost, number of species not flagged
    ``low_weight``, number of species, objective), then cut to ``num_alternatives``. The
    structure cost is ``species_penalty`` per significant species plus ``step_penalty`` per step
    of a significant species; with ``step_penalty`` 0 it orders like the number of significant
    species. Ranking on significant species means "one molecule plus weight-1 leftovers" beats
    "the molecule split in two" when the solver scores them equally. The top decomposition then
    gets an ``ambiguity`` report (module docstring).
    """
    params = params or DecompositionParams()
    decomposer = CycleDecomposer(graph, flow, params, paths)
    if decomposer.reads:
        log.info("read constraints: %d junction sequences (%.0f reads); %d paths off the flow",
                 len(decomposer.reads), sum(r.support for r in decomposer.reads), decomposer.reads_skipped)
    try:
        raw = decomposer.solve()
    except RuntimeError:
        if not (params.path_fraction > 0 and decomposer.reads):
            raise
        log.warning("no decomposition explains %.0f%% of the read support; solving without --path-fraction",
                    100 * params.path_fraction)
        params = replace(params, path_fraction=0.0)
        decomposer = CycleDecomposer(graph, flow, params, paths)
        raw = decomposer.solve()
    candidates = []
    for objective, species in raw:
        species = sorted(species, key=lambda s: (-s.weight, -sum(s.counts.values())))
        _flag_low_weight(species, params)
        residual = decomposer._residual(species)
        significant = sum(1 for s in species if not s.low_weight)
        candidates.append((round(residual, 6), round(structure_cost(species, params), 6), significant, len(species),
                           round(objective, 9), objective, species))
    candidates.sort(key=lambda c: c[:5])
    result = []
    for rank, (residual, _, _, _, _, objective, species) in enumerate(candidates[:max(1, params.num_alternatives)], start=1):
        annotated = [_annotate(graph, s, params, paths) for s in species]
        call, fraction = classify(annotated, params)
        total, explained = explained_reads(graph, annotated, decomposer.reads)
        result.append(Decomposition(rank=rank, objective=objective, residual=decomposer._residual(annotated), species=annotated,
                                    call=call, circular_fraction=fraction, reads_total=total, reads_explained=explained,
                                    steps=_steps([s for s in annotated if not s.low_weight])))
    if result and result[0].species and params.ambiguity_window > 0:
        result[0].ambiguity = ambiguity_report(graph, flow, params, paths, result)
    return result


def ambiguity_report(graph: ECGraph, flow: Dict[EdgeKey, int], params: DecompositionParams, paths, ranked: List[Decomposition]) -> Dict[str, object]:
    """Alternatives to the top decomposition with exactly the same fit to the flow (module docstring).

    ``gap`` is an alternative's score minus the top decomposition's, where the score is the MILP
    objective plus ``step_penalty`` x the steps of the significant species: positive means the
    priors (and reads) prefer the top one by that much, negative that the top one was chosen by
    the ranking (fewer significant species) over a lower score. ``decided_by`` names what
    separates the top decomposition from its closest alternative."""
    top = ranked[0]
    top_score = top.objective + params.step_penalty * top.steps
    mine = _structures(top.species)
    alts, seen = [], {mine}
    for objective, species in CycleDecomposer(graph, flow, params, paths).same_fit_alternatives(top.species):
        species = sorted(species, key=lambda s: (-s.weight, -sum(s.counts.values())))
        _flag_low_weight(species, params)
        if _structures(species) in seen:  # sorted by objective: keep each set of structures at its best
            continue
        seen.add(_structures(species))
        match = next((d.rank for d in ranked if _structures(d.species) == _structures(species)), None)
        significant = [s for s in species if not s.low_weight]
        score = objective + params.step_penalty * _steps(significant)
        alts.append({"objective": round(objective, 6), "gap": round(score - top_score, 6),
                     "n_species": len(significant), "steps": _steps(significant), "rank": match,
                     "weights": [s.weight for s in significant]})
    report = {"ambiguous": bool(alts), "window": params.ambiguity_window, "n_alternatives": len(alts),
              "n_species": sum(1 for s in top.species if not s.low_weight), "steps": top.steps, "alternatives": alts[:10]}
    if alts:
        closest = min(alts, key=lambda a: abs(a["gap"]))
        report["gap"] = closest["gap"]
        if closest["n_species"] != report["n_species"]:
            report["decided_by"] = "species_penalty"
        elif closest["steps"] != top.steps and params.step_penalty > 0:
            report["decided_by"] = "step_penalty"
        elif closest["steps"] != top.steps:
            report["decided_by"] = "tie-break"
        else:
            report["decided_by"] = "reads or tie-break"
        log.warning("top decomposition is ambiguous: %d alternative(s) fit the flow exactly as well; the closest "
                    "(%d species, %d steps, objective gap %+.4f) was set apart by %s", len(alts), closest["n_species"],
                    closest["steps"], closest["gap"], report["decided_by"])
    return report


def explained_reads(graph: ECGraph, species: Sequence[Species], reads: Sequence[ReadConstraint]) -> Tuple[float, float]:
    """(read support, share found in order in some species' top ordering)."""
    total = sum(r.support for r in reads)
    if total <= 0:
        return 0.0, 0.0
    hit = sum(r.support for r in reads if any(path_matches(graph, s.orderings[0].steps, r.path) for s in species if s.orderings))
    return total, hit / total
