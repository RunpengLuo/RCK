"""Heaviest-first decomposition of per-clone RCK karyotypes into cycles and telomere-to-telomere paths.

Step 4 of the multi-clone plan (``docs/multiclone_plan.md`` 4.3; ``docs/ecDNA.md`` 11.3). This is
AmpliconArchitect's rule (extract the heaviest walk, subtract it, repeat) in the form CoRAL 3.1.1
implements it as its greedy mode (``coral/cycle_decomposition.py::maximize_weights_greedy``,
``coral/models/concrete.py::get_greedy_objective``): each round solves a one-walk MILP that
maximises the walk's length-weighted copy number, plus a bonus for long-read paths not yet
explained, with the walk's copies capped by the copy number left on every edge.

**Which walk is heaviest** (``heaviest``). ``copies`` (default) is AmpliconArchitect's rule: the
walk with the most copies, sum_k psi_k, wins; ties are broken by CoRAL's length-weighted copy
number. ``length`` is CoRAL's rule alone. They differ when a long low-copy path competes with a
short high-copy cycle. The DIAG has 1 Mb flank segments (CoRAL's graphs hold only the amplified
intervals), so with ``length`` a chromosome path through both flanks can outweigh a short ecDNA
and is extracted first, looping through the amplicon up to ``max_multiplicity`` times (an HSR-like
walk) and taking copies away from the cycle. ``copies`` takes the cycle first (tests/ecdna/test_greedy.py).

Multi-clone version. The input is every clone's integer karyotype from
:class:`rck.core.ilp_gurobi_ecdna.OptModelMultiCloneEcDNA` (RCK's model) on the annotated DIAG:
segment copy numbers c_{k,s,h}, phased adjacency copy numbers N_{k,a,phi}, and the copy-number
excess at telomere positions (where RCK lets chromosomes end). One round finds one walk sigma on the
haplotype-labelled DIAG, shared by the clones, with

* integer multiplicities m_e in 0..``max_multiplicity`` (binary expansion, as in
  :mod:`rck.ecdna.decompose`), alternation and balance at every phased extremity, connectivity
  by lazy cuts;
* either no ends (a **cycle**: an ecDNA candidate) or exactly two, at telomere positions (a
  **path**: a linear chromosome piece), through one virtual end node as in ``decompose.py``;
* integer copies psi_k >= 0 per clone with psi_k * m_e <= R_k(e) on every edge (R_k: what is left
  of clone k's karyotype), the products linearised exactly;
* objective: with ``heaviest='copies'``, first maximise sum_k psi_k; then, with that fixed,
  maximise sum_k psi_k * sum_seg m_e * len_e (CoRAL's length-weighted copy number, summed over the
  clones) + CoRAL's read bonus, minus a tie-break below one bp per edge traversal, so that among
  equally heavy walks the one with fewer traversals wins. With ``heaviest='length'`` only the
  second objective is used.

The walk is subtracted from every clone and the next round starts, until nothing is left (RCK's
karyotype is an exact balanced integer flow, so it can be used up exactly), the remaining
length-weighted copy number falls to ``stop_fraction`` of the total (CoRAL's ``p_total_weight``
rule, off by default), or ``max_walks`` walks were taken.

Deviations from CoRAL's greedy model (also in ``docs/ecDNA.md`` 11.3): the copy-count criterion
first (AA's rule; ``heaviest='copies'``); integer per-clone copies
instead of one fractional weight; walk ends only at RCK's telomere positions (CoRAL: any source
edge); multiplicity capped at ``max_multiplicity``; the traversal tie-break; long-read groups are
RCK's molecule groups (``rck.adjgroups.tsv``), each satisfied when the walk uses all its novel
adjacencies (CoRAL: an ordered subpath).
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from functools import reduce
from typing import Dict, List, Optional, Sequence, Tuple

from rck.core.structures import AdjacencyType, Haplotype, Phasing, AdjacencyGroupType
from rck.ecdna.decompose import DecompositionParams, Decomposition, Species, _annotate, _components, classify
from rck.ecdna.gurobi_env import close_model
from rck.ecdna.phased import ADJ, END, END_NODE, SEG, EdgeKey, adj_key, end_key, endpoints, is_adj_like, is_balanced, n_ends, \
    seg_key, segment_key_at
from rck.ecdna.subgraph_io import ALL_PHASINGS, HEAD, IMPUTED, NOVEL, REFERENCE, TAIL, Adjacency, ECGraph, Node, Segment

log = logging.getLogger("rck.ecdna")

HAPS = ("A", "B")
PHASING_NAMES = tuple(str(ph) for ph in ALL_PHASINGS)
CYCLE, PATH = "cycle", "path"


@dataclass
class GreedyParams:
    heaviest: str = "copies"  # copies: AA's rule (most copies, sum_k psi_k, then CoRAL's weight) | length: CoRAL's weight only
    max_walks: int = 50
    max_multiplicity: int = 3
    read_alpha: float = 0.01  # CoRAL's alpha: weight of unexplained long-read groups (0 = off)
    stop_fraction: float = 0.0  # stop when remaining length-weighted CN <= this share of the total (CoRAL 0.1)
    min_ecdna_copies: int = 2
    min_species_weight_frac: float = 0.10
    max_orderings: int = 20
    time_limit: Optional[float] = None
    threads: int = 2
    log_file: Optional[str] = None
    log_to_console: bool = False


@dataclass
class Walk:
    kind: str  # cycle | path
    counts: Dict[EdgeKey, int]
    psi: Dict[str, int]
    order: int
    objective: float = 0.0
    groups: List[str] = field(default_factory=list)  # molecule groups whose novel adjacencies it uses
    species: Optional[Species] = None


@dataclass
class GreedyResult:
    graph: ECGraph
    walks: List[Walk]
    total_weight: float
    remaining_weight: float
    remaining: Dict[str, Dict[EdgeKey, int]]
    stopped: str


# ----------------------------------------------------------------------------- graph and flow
def _nid(position) -> str:
    return "{c}:{p}_{s}".format(c=position.chromosome, p=position.coordinate, s=str(position.strand))


def diag_graph(inp, name: str = "amplicon") -> Tuple[ECGraph, Dict[str, str]]:
    """The DIAG as an :class:`ECGraph` (so the walk machinery of ``decompose.py``/``walk.py`` applies),
    and RCK's ``stable_id_non_phased`` -> aid. Segment ids are RCK's ``stable_id_non_hap``; node ids
    are ``<chrom>:<coord>_<strand>`` of RCK's positions (``-`` the left, ``+`` the right extremity)."""
    telomeres = {_nid(p) for p in inp.telomeres}
    nodes, segments, adjacencies, aid_of = {}, {}, {}, {}
    for s in inp.segments:
        sid = s.stable_id_non_hap
        tail, head = _nid(s.start_position), _nid(s.end_position)
        nodes[tail] = Node(nid=tail, chrom=s.chromosome, pos=s.start_coordinate, side=TAIL, segment=sid, terminal=tail in telomeres)
        nodes[head] = Node(nid=head, chrom=s.chromosome, pos=s.end_coordinate, side=HEAD, segment=sid, terminal=head in telomeres)
        segments[sid] = Segment(sid=sid, chrom=s.chromosome, start=s.start_coordinate, end=s.end_coordinate, tail=tail, head=head, depth=0.0)
    imputed = set(inp.imputed_aids)
    for a in inp.adjacencies:
        rid = a.stable_id_non_phased
        aid = str(a.extra.get("aid", rid))
        if a.adjacency_type == AdjacencyType.REFERENCE:
            kind = REFERENCE
        else:
            kind = IMPUTED if rid in imputed else NOVEL
        aid_of[rid] = aid
        adjacencies[aid] = Adjacency(aid=aid, u=_nid(a.position1), v=_nid(a.position2), kind=kind, allowed=ALL_PHASINGS, preferred=ALL_PHASINGS)
    return ECGraph(nodes=nodes, segments=segments, adjacencies=adjacencies, name=name), aid_of


def clone_flows(inp, graph: ECGraph, aid_of: Dict[str, str], scnt, acnt) -> Dict[str, Dict[EdgeKey, int]]:
    """Every clone's karyotype as phased edge copy numbers, with the excess at each extremity as a
    linear end (``("T", nid, hap)``); RCK allows excess only at telomere positions (and, with
    ``--mc-linear-ends``, at the penalised extra ends)."""
    flows = {}
    for cid in inp.clone_ids:
        flow: Dict[EdgeKey, int] = {}
        for s in inp.segments:
            for h in (Haplotype.A, Haplotype.B):
                cn = int(scnt[cid].get_cn(sid=s.stable_id_non_hap, haplotype=h, default=0))
                if cn > 0:
                    flow[seg_key(s.stable_id_non_hap, h.value)] = cn
        for a in inp.adjacencies:
            for ph in ALL_PHASINGS:
                cn = int(acnt[cid].get_cn(aid=a.stable_id_non_phased, phasing=ph, default=0))
                if cn > 0:
                    flow[adj_key(aid_of[a.stable_id_non_phased], ph)] = cn
        load: Dict[Tuple[str, str], int] = {}
        for key, value in flow.items():
            if is_adj_like(key):
                for x in endpoints(graph, key):
                    load[x] = load.get(x, 0) + value
        for nid in graph.nodes:
            for h in HAPS:
                sk = segment_key_at(graph, (nid, h))
                excess = flow.get(sk, 0) - load.get((nid, h), 0)
                if excess < 0:
                    raise ValueError("clone {c}: negative excess {e} at {n} {h}".format(c=cid, e=excess, n=nid, h=h))
                if excess > 0:
                    flow[end_key(nid, h)] = excess
        flows[cid] = flow
    return flows


def _seg_len(graph: ECGraph, key: EdgeKey) -> int:
    return graph.segments[key[1]].length if key[0] == SEG else 0


def weighted(graph: ECGraph, flows: Dict[str, Dict[EdgeKey, int]]) -> float:
    return float(sum(v * _seg_len(graph, k) for flow in flows.values() for k, v in flow.items()))


# ------------------------------------------------------------------------------ read groups
def molecule_groups(inp, aid_of: Dict[str, str]) -> Dict[str, List[str]]:
    """RCK molecule groups (long-read paths with >= 2 novel adjacencies) as gid -> aids."""
    groups = {}
    for ag in inp.adjacency_groups:
        if ag.group_type != AdjacencyGroupType.MOLECULE:
            continue
        aids = sorted({aid_of[a.stable_id_non_phased] for a in ag.adjacencies if a.stable_id_non_phased in aid_of})
        if len(aids) >= 2:
            groups[str(ag.gid)] = aids
    return groups


# ------------------------------------------------------------------------------- one round
class HeaviestWalk:
    """The one-walk MILP of one greedy round."""

    def __init__(self, graph: ECGraph, remaining: Dict[str, Dict[EdgeKey, int]], params: GreedyParams,
                 groups: Optional[Dict[str, List[str]]] = None, group_bonus: float = 0.0):
        self.graph, self.params = graph, params
        self.clones = sorted(remaining)
        self.remaining = remaining
        self.edges = sorted({e for flow in remaining.values() for e, v in flow.items() if v > 0})
        self.groups = groups or {}
        self.group_bonus = group_bonus
        self.model = None

    def build(self):
        import gurobipy as gp
        from gurobipy import GRB

        p, g, E = self.params, self.graph, self.edges
        m = gp.Model("rck-ecdna-greedy")
        self.model = m
        m.Params.OutputFlag = 1 if (p.log_to_console or p.log_file) else 0
        m.Params.LogToConsole = 1 if p.log_to_console else 0
        if p.log_file:
            m.Params.LogFile = p.log_file
        if p.time_limit:
            m.Params.TimeLimit = p.time_limit
        m.Params.Threads = p.threads
        m.Params.MIPGap = 0.0
        m.Params.LazyConstraints = 1
        nbits = max(1, int(math.ceil(math.log2(p.max_multiplicity + 1))))
        self.nbits = nbits
        W = {k: max([v for v in self.remaining[k].values()] + [0]) for k in self.clones}
        self.psi = {k: m.addVar(lb=0, ub=W[k], vtype=GRB.INTEGER, name="psi[{k}]".format(k=k)) for k in self.clones}
        self.beta, self.y, self.mult = {}, {}, {}
        cover = {}
        for i, e in enumerate(E):
            bits = [m.addVar(vtype=GRB.BINARY, name="beta[{i},{b}]".format(i=i, b=b)) for b in range(nbits)]
            for b, beta in enumerate(bits):
                self.beta[e, b] = beta
            mult = gp.quicksum((2 ** b) * bits[b] for b in range(nbits))
            self.mult[e] = mult
            y = m.addVar(vtype=GRB.BINARY, name="y[{i}]".format(i=i))
            self.y[e] = y
            m.addConstr(mult <= p.max_multiplicity * y)
            m.addConstr(y <= mult)
            for k in self.clones:
                terms = []
                for b in range(nbits):
                    zeta = m.addVar(lb=0, ub=W[k], name="zeta[{k},{i},{b}]".format(k=k, i=i, b=b))
                    m.addConstr(zeta <= W[k] * bits[b])
                    m.addConstr(zeta <= self.psi[k])
                    m.addConstr(zeta >= self.psi[k] - W[k] * (1 - bits[b]))
                    terms.append((2 ** b) * zeta)
                cover[k, e] = gp.quicksum(terms)
                m.addConstr(cover[k, e] <= self.remaining[k].get(e, 0), name="cap[{k},{i}]".format(k=k, i=i))
        self.cover = cover
        m.addConstr(gp.quicksum(self.y.values()) >= 2)
        m.addConstr(gp.quicksum(self.psi.values()) >= 1)
        ends = [e for e in E if e[0] == END]
        self.linear = m.addVar(vtype=GRB.BINARY, name="linear")
        m.addConstr(gp.quicksum(self.mult[e] for e in ends) == 2 * self.linear)
        nodes: Dict[Tuple[str, str], List[EdgeKey]] = {}
        for e in E:
            if is_adj_like(e):
                x, z = endpoints(g, e)
                nodes.setdefault(x, []).append(e)
                nodes.setdefault(z, []).append(e)
            else:
                for x in endpoints(g, e):
                    nodes.setdefault(x, [])
        nodes.pop(END_NODE, None)
        for j, (node, adjs) in enumerate(sorted(nodes.items())):
            sk = segment_key_at(g, node)
            seg_expr = self.mult[sk] if (sk is not None and sk in self.mult) else 0
            m.addConstr(seg_expr == gp.quicksum(self.mult[a] for a in adjs), name="bal[{j}]".format(j=j))

        weight = gp.quicksum(_seg_len(g, e) * cover[k, e] for k in self.clones for e in E if e[0] == SEG)
        self.weight = weight
        tie = 1.0 / (p.max_multiplicity * max(len(E), 1) + 1)
        terms = [weight, -tie * gp.quicksum(self.mult.values())]
        self.z = {}
        if self.groups and self.group_bonus > 0:
            for gid, aids in sorted(self.groups.items()):
                ors = []
                for aid in aids:
                    keys = [e for e in E if e[0] == ADJ and e[1] == aid]
                    if not keys:
                        ors = None
                        break
                    o = m.addVar(vtype=GRB.BINARY, name="o[{g},{a}]".format(g=gid, a=aid))
                    m.addConstr(o <= gp.quicksum(self.y[e] for e in keys))
                    ors.append(o)
                if not ors:
                    continue
                z = m.addVar(vtype=GRB.BINARY, name="z[{g}]".format(g=gid))
                for o in ors:
                    m.addConstr(z <= o)
                self.z[gid] = z
                terms.append(self.group_bonus * z)
        self.secondary = gp.quicksum(terms)
        self.copies = gp.quicksum(self.psi.values())
        if p.heaviest not in ("copies", "length"):
            raise ValueError("heaviest must be copies or length, not {h}".format(h=p.heaviest))
        m.setObjective(self.copies if p.heaviest == "copies" else self.secondary, GRB.MAXIMIZE)
        m.update()
        return m

    def _callback(self):
        import gurobipy as gp
        from gurobipy import GRB

        E, g = self.edges, self.graph

        def cb(model, where):
            if where != GRB.Callback.MIPSOL:
                return
            values = model.cbGetSolution([self.y[e] for e in E])
            used = [e for e, v in zip(E, values) if v > 0.5]
            comps = _components(g, used) if used else []
            if len(comps) <= 1:
                return
            for i in range(len(comps)):
                side = {x for e in comps[i] for x in endpoints(g, e)}
                cut = [f for f in E if sum(x in side for x in endpoints(g, f)) == 1]
                for j in range(len(comps)):
                    if i != j:
                        model.cbLazy(gp.quicksum(self.y[f] for f in cut) >= self.y[comps[i][0]] + self.y[comps[j][0]] - 1)

        return cb

    def solve(self) -> Optional[Walk]:
        from gurobipy import GRB

        if not self.edges:
            return None
        if self.model is None:
            self.build()
        m = self.model
        try:
            m.optimize(self._callback())
            if m.Status == GRB.INFEASIBLE or m.SolCount == 0:
                return None
            if self.params.heaviest == "copies":  # AA: most copies first; among those, CoRAL's weight, then fewer traversals
                best = int(round(m.ObjVal))
                m.addConstr(self.copies >= best, name="copies")
                m.setObjective(self.secondary, GRB.MAXIMIZE)
                m.optimize(self._callback())
                if m.SolCount == 0:
                    return None
            counts = {}
            for e in self.edges:
                v = int(round(sum((2 ** b) * self.beta[e, b].X for b in range(self.nbits))))
                if v > 0:
                    counts[e] = v
            psi = {k: int(round(self.psi[k].X)) for k in self.clones}
            objective = float(self.weight.getValue())
            groups = sorted(gid for gid, z in self.z.items() if z.X > 0.5)
        finally:
            close_model(m)
            self.model = None
        if not counts or sum(psi.values()) <= 0:
            return None
        g = reduce(math.gcd, counts.values())
        if g > 1 and not n_ends(counts):  # a circle walked g times is that circle at g times the copies
            counts = {e: v // g for e, v in counts.items()}
            psi = {k: v * g for k, v in psi.items()}
        problems = is_balanced(self.graph, counts)
        if problems or len(_components(self.graph, list(counts))) > 1 or n_ends(counts) not in (0, 2):
            raise RuntimeError("greedy walk is not one balanced walk: {p}".format(p=problems[:3]))
        return Walk(kind=PATH if n_ends(counts) else CYCLE, counts=counts, psi=psi, order=0, objective=objective, groups=groups)


# ------------------------------------------------------------------------------- the loop
def heaviest_first(graph: ECGraph, flows: Dict[str, Dict[EdgeKey, int]], params: GreedyParams,
                   groups: Optional[Dict[str, List[str]]] = None) -> GreedyResult:
    remaining = {k: dict(v) for k, v in flows.items()}
    total = weighted(graph, flows)
    unsatisfied = dict(groups or {})
    walks: List[Walk] = []
    stopped = "exhausted"
    while True:
        left = weighted(graph, remaining)
        if all(v == 0 for flow in remaining.values() for v in flow.values()):
            stopped = "exhausted"
            break
        if params.stop_fraction > 0 and left <= params.stop_fraction * total and not unsatisfied:
            stopped = "stop_fraction"
            break
        if len(walks) >= params.max_walks:
            stopped = "max_walks"
            break
        bonus = 0.0
        if unsatisfied and params.read_alpha > 0:  # CoRAL: pp = alpha x remaining weight / unsatisfied, at least 1
            bonus = max(params.read_alpha * left / len(unsatisfied), 1.0)
        walk = HeaviestWalk(graph, remaining, params, unsatisfied, bonus).solve()
        if walk is None:
            stopped = "no_walk"
            break
        walk.order = len(walks) + 1
        for k, psi in walk.psi.items():
            for e, mult in walk.counts.items():
                remaining[k][e] = remaining[k].get(e, 0) - psi * mult
                if remaining[k][e] < 0:
                    raise RuntimeError("greedy walk over-uses {e} in clone {k}".format(e=e, k=k))
        for gid in walk.groups:
            unsatisfied.pop(gid, None)
        walks.append(walk)
        log.info("greedy walk %d: %s, psi %s, %d edges, weight %.0f", walk.order, walk.kind,
                 ",".join("{k}={v}".format(k=k, v=v) for k, v in sorted(walk.psi.items())), len(walk.counts), walk.objective)
    remaining = {k: {e: v for e, v in flow.items() if v} for k, flow in remaining.items()}
    return GreedyResult(graph=graph, walks=walks, total_weight=total, remaining_weight=weighted(graph, remaining),
                        remaining=remaining, stopped=stopped)


# ----------------------------------------------------------------------------- reporting
def as_decomposition(result: GreedyResult, proportions: Dict[str, float], params: GreedyParams) -> Decomposition:
    """The walks as one :class:`Decomposition` in extraction order. A species' ``weight`` is its bulk
    copies per cell, sum_k u_k psi_k (comparable to the truth's bulk weights); paths are
    ``chromosomal`` and cycles lighter than ``min_species_weight_frac`` of the heaviest cycle
    ``low_weight``, as in the bulk route."""
    dparams = DecompositionParams(max_orderings=params.max_orderings, min_ecdna_copies=params.min_ecdna_copies,
                                  min_species_weight_frac=params.min_species_weight_frac)
    species = []
    for w in result.walks:
        bulk = sum(proportions.get(k, 0.0) * v for k, v in w.psi.items())
        sp = Species(weight=bulk, counts=w.counts, weight_int=None)
        _annotate(result.graph, sp, dparams, None)
        w.species = sp
        species.append(sp)
    cycles = [w for w in result.walks if w.kind == CYCLE]
    heaviest = max([w.species.weight for w in cycles] + [0.0])
    for w in cycles:
        w.species.low_weight = len(cycles) > 1 and w.species.weight < max(1e-9, params.min_species_weight_frac * heaviest)
    call, fraction = classify([w.species for w in cycles], dparams)
    return Decomposition(rank=1, objective=result.total_weight - result.remaining_weight, residual=result.remaining_weight,
                         species=species, call=call, circular_fraction=fraction)


def decomposition_frame(result: GreedyResult, decomposition: Decomposition, clones: Sequence[str], proportions: Dict[str, float]):
    """``decompositions.tsv``: the bulk route's columns (``output.decomposition_table``), one species
    per walk in extraction order, plus ``kind``, ``order``, ``psi_<clone>``, ``proportion_<clone>``
    (psi over the clone's cycles) and ``proportion_bulk`` (sum_k u_k psi_k over the cycles). Paths
    get ``status`` ``chromosomal`` and no proportions."""
    from rck.ecdna.output import decomposition_table

    frame = decomposition_table(result.graph, [decomposition])
    names, n_cycle, n_path = {}, 0, 0
    for i, w in enumerate(result.walks, start=1):
        if w.kind == CYCLE:
            n_cycle += 1
            names["ec{i}".format(i=i)] = "ec{n}".format(n=n_cycle)
        else:
            n_path += 1
            names["ec{i}".format(i=i)] = "path{n}".format(n=n_path)
    by_name = {"ec{i}".format(i=i): w for i, w in enumerate(result.walks, start=1)}
    cycles = [w for w in result.walks if w.kind == CYCLE]
    clone_total = {k: sum(w.psi.get(k, 0) for w in cycles) for k in clones}
    bulk_total = sum(proportions.get(k, 0.0) * w.psi.get(k, 0) for w in cycles for k in clones)
    extra = {"kind": [], "order": []}
    for k in clones:
        extra["psi_" + k] = []
        extra["proportion_" + k] = []
    extra["proportion_bulk"] = []
    status = []
    for name in frame["species"]:
        w = by_name[name]
        extra["kind"].append(w.kind)
        extra["order"].append(w.order)
        for k in clones:
            extra["psi_" + k].append(w.psi.get(k, 0))
            extra["proportion_" + k].append(w.psi.get(k, 0) / clone_total[k] if (w.kind == CYCLE and clone_total[k]) else None)
        bulk = sum(proportions.get(k, 0.0) * w.psi.get(k, 0) for k in clones)
        extra["proportion_bulk"].append(bulk / bulk_total if (w.kind == CYCLE and bulk_total) else None)
        status.append("chromosomal" if w.kind == PATH else w.species.status)
    frame["status"] = status
    frame["species"] = [names[n] for n in frame["species"]]
    for column, values in extra.items():
        frame[column] = values
    return frame
