"""Decompose the ecDNA flow into the fewest circular species, with alternatives.

Given the ecDNA-layer flow ``F_e`` on phased edges (from :mod:`rck.ecdna.ilp_circular`), find
up to ``K`` species slots. Slot ``k`` has an integer copy weight ``w_k`` and a small integer
multiplicity ``m_ke`` (0..``max_multiplicity``) on every edge, such that

* ``sum_k w_k * m_ke`` reproduces ``F_e`` (deviations are penalized, not forbidden);
* each slot is balanced with alternation at every phased node
  (segment multiplicity == sum of incident adjacency multiplicities, self-loops twice);
* each slot's support is connected (lazy cuts in a Gurobi callback), so with balance it is
  one alternating closed walk.

The product ``w_k * m_ke`` is linearized exactly through the binary expansion of ``m_ke``
(McCormick on binary x integer). Objective: ``sum_e omega_e |residual_e| + species_penalty *
#species``. Slots are ordered by weight to break symmetry.

Alternatives come from Gurobi's solution pool (``PoolSearchMode=2``) or, with
``mode='nogood'``, from re-solving with a no-good cut on the slot supports. Pool solutions
are re-checked for balance and connectivity and deduplicated by their species multisets.
"""

from __future__ import annotations

import logging
import math
from functools import reduce
from dataclasses import dataclass, field
from typing import Dict, FrozenSet, List, Optional, Sequence, Tuple

from rck.ecdna.gurobi_env import close_model, default_threads
from rck.ecdna.phased import ADJ, SEG, EdgeKey, PhasedNode, endpoints, is_balanced, segment_key_at
from rck.ecdna.subgraph_io import IMPUTED, NOVEL, ECGraph
from rck.ecdna.walk import ReadPath, Step, alternating_circuit, enumerate_circuits, path_support, walk_length

log = logging.getLogger(__name__)


@dataclass
class DecompositionParams:
    max_species: int = 4
    max_multiplicity: int = 3
    species_penalty: float = 1.0
    multiplicity_penalty: float = 1e-4  # tie-break: prefer the walk with fewer edge traversals
    min_species_weight: Optional[int] = None  # species below this weight are flagged low_weight
    min_species_weight_frac: float = 0.10  # default threshold: this fraction of the heaviest species
    adj_residual_weight: float = 1.0
    seg_length_cap: int = 50_000
    min_seg_weight: float = 0.05
    num_alternatives: int = 10
    pool_gap: float = 0.10
    mode: str = "pool"  # pool | nogood
    max_orderings: int = 20
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
    weight: int
    counts: Dict[EdgeKey, int]
    orderings: List[Ordering] = field(default_factory=list)
    length_bp: int = 0
    imputed: List[str] = field(default_factory=list)
    junctions: List[str] = field(default_factory=list)
    low_weight: bool = False  # below the minimum species weight: likely flow residue, not a molecule

    @property
    def status(self) -> str:
        return "low_weight" if self.low_weight else "species"

    @property
    def signature(self) -> Tuple[int, FrozenSet]:
        return self.weight, frozenset(self.counts.items())


@dataclass
class Decomposition:
    rank: int
    objective: float
    residual: float
    species: List[Species]

    @property
    def signature(self) -> FrozenSet:
        return frozenset(s.signature for s in self.species)


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
    def __init__(self, graph: ECGraph, flow: Dict[EdgeKey, int], params: Optional[DecompositionParams] = None):
        self.graph = graph
        self.flow = {k: int(v) for k, v in flow.items() if v > 0}
        self.params = params or DecompositionParams()
        self.edges = sorted(self.flow)
        self.index = {e: i for i, e in enumerate(self.edges)}
        self.model = None
        self.discarded = 0

    def _omega(self, key: EdgeKey) -> float:
        p = self.params
        if key[0] == SEG:
            s = self.graph.segments[key[1]]
            return max(min(s.length / float(p.seg_length_cap), 1.0), p.min_seg_weight)
        return p.adj_residual_weight

    def build(self):
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
        self.beta, self.zeta, self.y, self.mult = {}, {}, {}, {}
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
            # balance with alternation, per phased node
            nodes: Dict[PhasedNode, List[EdgeKey]] = {}
            for e in E:
                if e[0] == ADJ:
                    x, z = endpoints(g, e)
                    nodes.setdefault(x, []).append(e)
                    nodes.setdefault(z, []).append(e)
                else:
                    for x in endpoints(g, e):
                        nodes.setdefault(x, [])
            for j, (node, adjs) in enumerate(sorted(nodes.items())):
                sk = segment_key_at(g, node)
                seg_expr = self.mult[k, sk] if (sk is not None and sk in self.flow) else 0
                m.addConstr(seg_expr == gp.quicksum(self.mult[k, a] for a in adjs), name="bal[{k},{j}]".format(k=k, j=j))
        for k in range(K - 1):
            m.addConstr(self.w[k] >= self.w[k + 1])
            m.addConstr(self.u[k] >= self.u[k + 1])

        self.res = {}
        terms = []
        for e in E:
            cover = gp.quicksum((2 ** b) * self.zeta[k, e, b] for k in range(K) for b in range(nbits))
            r = m.addVar(lb=0, name="r[{i}]".format(i=self.index[e]))
            m.addConstr(r >= cover - self.flow[e])
            m.addConstr(r >= self.flow[e] - cover)
            self.res[e] = r
            terms.append(self._omega(e) * r)
        terms.append(p.species_penalty * gp.quicksum(self.u.values()))
        if p.multiplicity_penalty:
            terms.append(p.multiplicity_penalty * gp.quicksum(self.mult[k, e] for k in range(K) for e in E))
        m.setObjective(gp.quicksum(terms), GRB.MINIMIZE)
        m.update()
        return m

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
            if g > 1:
                counts = {e: v // g for e, v in counts.items()}
                weight *= g
            if is_balanced(self.graph, counts) or len(_components(self.graph, list(counts))) > 1:
                valid = False
            species.append(Species(weight=weight, counts=counts))
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
                if best is None:
                    best = objective
                elif objective > best + abs(best) * p.pool_gap + 1e-6:
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
        return solutions

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


def _flag_low_weight(species: List[Species], params: DecompositionParams) -> None:
    if not species:
        return
    threshold = params.min_species_weight
    if threshold is None:
        threshold = max(1, int(math.ceil(params.min_species_weight_frac * species[0].weight)))
    for s in species:
        s.low_weight = len(species) > 1 and s.weight < threshold


def decompose(graph: ECGraph, flow: Dict[EdgeKey, int], params: Optional[DecompositionParams] = None,
              paths: Optional[List[ReadPath]] = None) -> List[Decomposition]:
    """Ranked decompositions of ``flow`` into circular species, each with ranked orderings.

    Decompositions are ranked by (residual, number of species not flagged ``low_weight``,
    number of species, objective), then cut to ``num_alternatives``. Ranking on significant
    species first means "one molecule plus weight-1 leftovers" beats "the molecule split in
    two" when the solver scores them equally.
    """
    params = params or DecompositionParams()
    decomposer = CycleDecomposer(graph, flow, params)
    raw = decomposer.solve()
    candidates = []
    for objective, species in raw:
        species = sorted(species, key=lambda s: (-s.weight, -sum(s.counts.values())))
        _flag_low_weight(species, params)
        residual = decomposer._residual(species)
        significant = sum(1 for s in species if not s.low_weight)
        candidates.append((round(residual, 6), significant, len(species), round(objective, 9), objective, species))
    candidates.sort(key=lambda c: c[:4])
    result = []
    for rank, (residual, _, _, _, objective, species) in enumerate(candidates[:max(1, params.num_alternatives)], start=1):
        annotated = [_annotate(graph, s, params, paths) for s in species]
        result.append(Decomposition(rank=rank, objective=objective, residual=decomposer._residual(annotated), species=annotated))
    return result
