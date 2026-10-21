"""Two-layer circular flow MILP: haplotype-specific ecDNA flow on top of a chromosomal background.

This adapts RCK's ``OptModelMultiClone`` (``rck/core/ilp_gurobi.py``) to circular DNA. The
pieces kept from RCK are its variable layout (integer copy numbers per segment haplotype and
per adjacency labeling, binary presence indicators tied to them with big-M constraints, at
most one labeling per novel adjacency, AA/BB-only reference adjacencies) and its
haplotype-aware balance constraint (``define_constraints_on_nodes`` /
``get_lin_expr_for_node_haplotype``, self-loops counted twice). What changes:

* **No telomeres in the ecDNA layer.** Every extremity is strictly balanced
  (segment copy number == sum of incident adjacency copy numbers, per haplotype), so by the
  balance theorem the ecDNA layer decomposes into alternating closed walks.
* **A background layer.** Observed copy number = ecDNA + chromosomal copies. The background
  layer is balanced too, but may start / end at terminal extremities (region flanks), may
  not use imputed edges, and is softly tied to the copy number measured just outside the
  region. Without it the chromosomal copies under an amplicon would either be forced into
  cycles or charged to the deviation term.
* **Imputed candidate edges** carry a usage penalty and never count as observed junctions.
* **Objective.** Length-capped, sqrt-normalized deviation between observed and inferred
  allele-specific segment copy numbers, plus junction read support, plus penalties for
  imputed edges, novel-adjacency use and (optionally) disfavoured junction labelings.
  RCK's novel-adjacency utilization constraint is dropped: chromosomal SVs that touch the
  amplicon should be free to stay unused.
* Optional per-haplotype extremity exclusivity (``enforce_is``), RCK's infinite-sites
  constraint, off by default because chromothripsis-derived ecDNA violates it.
* **Linear ends (optional).** With ``linear_ends`` set to ``boundaries`` (default) or ``all``,
  amplified copies may also start or stop at selected extremities, at ``linear_end_penalty``
  per copy and end. ``boundaries`` allows ends only where copy number steps between a segment
  and its reference neighbour, at region flanks, and at segment ends without a reference join
  (the places AmpliconArchitect and CoRAL connect to their source / sink). This lets the model
  explain an amplicon that has no closing junction as linear instead of forcing a circle or an
  imputed join; ``none`` restores the strictly circular model. Penalty trade-off: a circle
  through observed junctions is always cheapest; an imputed closing edge costs
  ``imputed_penalty`` once, a linear walk ``2 x linear_end_penalty`` per copy, so imputation
  wins for high-copy amplicons and linear ends for low-copy ones.
"""

from __future__ import annotations

import logging
import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import pandas as pd

from rck.core.structures import Haplotype, Phasing
from rck.ecdna.gurobi_env import close_model, default_threads
from rck.ecdna.phased import adj_key, end_key, seg_key
from rck.ecdna.scale import Observations
from rck.ecdna.subgraph_io import HAPLOTYPES, IMPUTED, REFERENCE, ECGraph

log = logging.getLogger(__name__)


@dataclass
class CircularParams:
    ub_factor: float = 1.5  # segment copy-number upper bound = ceil(observed * ub_factor) + ub_add
    ub_add: int = 3
    bg_max: int = 4  # per-haplotype background copy-number cap
    alpha_novel: float = 0.5  # weight of the junction read-support term, novel adjacencies
    alpha_ref: float = 0.25  # same for reference adjacencies (#REF)
    imputed_penalty: float = 2.0  # per imputed edge used
    novel_penalty: float = 0.05  # per novel adjacency used (parsimony)
    soft_junction_penalty: float = 1.0  # per junction labeling outside the read-based label (soft mode)
    bg_prior_weight: float = 0.1  # pull of the background layer towards the flank copy number
    seg_length_cap: int = 50_000  # segment weight = min(length / cap, 1)
    min_seg_weight: float = 0.02
    enforce_is: bool = False
    allow_background_novel: bool = True  # chromosomal SVs may carry background flow
    linear_ends: str = "boundaries"  # boundaries | all | none: where ecDNA-layer copies may start / stop
    linear_end_penalty: float = 0.25  # per copy per linear end
    bg_balance_weight: float = 0.01  # tie-break: background copies split evenly between haplotypes
    boundary_min_copies: float = 2.0  # copy-number step that makes an extremity a boundary ...
    boundary_frac: float = 0.25  # ... or this fraction of the larger copy number, whichever is larger
    time_limit: Optional[float] = None
    threads: Optional[int] = None
    mip_gap: float = 1e-4
    log_file: Optional[str] = None
    log_to_console: bool = False


@dataclass
class FlowSolution:
    status: str
    objective: float
    mip_gap: float
    seg_ec: Dict[Tuple[str, Haplotype], int]
    seg_bg: Dict[Tuple[str, Haplotype], int]
    adj_ec: Dict[Tuple[str, Phasing], int]
    adj_bg: Dict[Tuple[str, Phasing], int]
    terminal_bg: Dict[Tuple[str, Haplotype], int]
    imputed_used: List[str] = field(default_factory=list)
    linear_ends: Dict[Tuple[str, Haplotype], int] = field(default_factory=dict)

    def ecdna_edges(self) -> Dict[Tuple[str, str, str], int]:
        """Positive ecDNA-layer flow on phased edge keys (see :mod:`rck.ecdna.phased`)."""
        result = {}
        for (sid, h), v in self.seg_ec.items():
            if v > 0:
                result[seg_key(sid, h)] = v
        for (aid, ph), v in self.adj_ec.items():
            if v > 0:
                result[adj_key(aid, ph)] = v
        for (nid, h), v in self.linear_ends.items():
            if v > 0:
                result[end_key(nid, h)] = v
        return result

    def segment_table(self, graph: ECGraph, obs: Observations) -> pd.DataFrame:
        rows = []
        for sid, s in sorted(graph.segments.items(), key=lambda kv: (kv[1].chrom, kv[1].start)):
            so = obs.segments[sid]
            rows.append({
                "sid": sid, "chrom": s.chrom, "start": s.start, "end": s.end, "length": s.length,
                "obs_total": so.total, "obs_A": so.a, "obs_B": so.b, "obs_source": so.source,
                "ecdna_A": self.seg_ec[(sid, Haplotype.A)], "ecdna_B": self.seg_ec[(sid, Haplotype.B)],
                "bg_A": self.seg_bg[(sid, Haplotype.A)], "bg_B": self.seg_bg[(sid, Haplotype.B)],
            })
        frame = pd.DataFrame(rows)
        if len(frame):
            frame["inferred_total"] = frame[["ecdna_A", "ecdna_B", "bg_A", "bg_B"]].sum(axis=1)
        return frame

    def adjacency_table(self, graph: ECGraph, obs: Observations) -> pd.DataFrame:
        rows = []
        for aid, a in sorted(graph.adjacencies.items()):
            for ph in a.allowed:
                ec = self.adj_ec.get((aid, ph), 0)
                bg = self.adj_bg.get((aid, ph), 0)
                rows.append({"aid": aid, "kind": a.kind, "u": a.u, "v": a.v, "phasing": str(ph), "reads": a.reads,
                             "obs_cn": obs.adjacencies.get(aid), "junction_hap": a.junction_hap, "ecdna": ec, "bg": bg})
        for (nid, h), v in sorted(self.linear_ends.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
            if v > 0:
                rows.append({"aid": "end:{n}".format(n=nid), "kind": "linear_end", "u": nid, "v": None, "phasing": str(h),
                             "reads": None, "obs_cn": None, "junction_hap": None, "ecdna": v, "bg": 0})
        return pd.DataFrame(rows)


class OptModelCircular:
    """Build and solve the two-layer flow MILP for one amplicon graph."""

    def __init__(self, graph: ECGraph, obs: Observations, params: Optional[CircularParams] = None):
        self.graph = graph
        self.obs = obs
        self.params = params or CircularParams()
        self.model = None
        self.xe, self.xb, self.ye, self.yb, self.pv, self.t, self.te = {}, {}, {}, {}, {}, {}, {}
        self.end_nodes: List[str] = []
        self.ub_seg: Dict[str, int] = {}

    # ---- helpers ---- #
    def _seg_weight(self, sid: str) -> float:
        s = self.graph.segments[sid]
        return max(min(s.length / float(self.params.seg_length_cap), 1.0), self.params.min_seg_weight)

    def _big_m(self, aid: str) -> int:
        a = self.graph.adjacencies[aid]
        bounds = []
        for nid in (a.u, a.v):
            sid = self.graph.nodes[nid].segment
            if sid is not None:
                bounds.append(self.ub_seg[sid])
        m = max(bounds) if bounds else self.params.bg_max
        if a.is_self_loop:
            m = max(m // 2, 1)
        return m

    def linear_end_nodes(self) -> List[str]:
        """Extremities where ecDNA-layer copies may start or stop (see module docstring)."""
        p, g, obs = self.params, self.graph, self.obs
        mode = (p.linear_ends or "none").lower()
        if mode == "none":
            return []
        candidates = sorted(nid for nid, n in g.nodes.items() if n.segment is not None)
        if mode == "all":
            return candidates
        if mode != "boundaries":
            raise ValueError("linear_ends must be boundaries, all or none, not {m}".format(m=p.linear_ends))
        result = []
        for nid in candidates:
            c = obs.segments[g.nodes[nid].segment].total
            refs = [a for a in g.adjacencies_at(nid) if a.kind == REFERENCE]
            if not refs:
                result.append(nid)
                continue
            for a in refs:
                other_sid = g.nodes[a.other_end(nid)].segment
                if other_sid is None:  # region flank
                    result.append(nid)
                    break
                oc = obs.segments[other_sid].total
                if abs(c - oc) >= max(p.boundary_min_copies, p.boundary_frac * max(c, oc)):
                    result.append(nid)
                    break
        return result

    def _uses_background(self, aid: str) -> bool:
        a = self.graph.adjacencies[aid]
        if a.kind == IMPUTED:
            return False
        return a.kind == REFERENCE or self.params.allow_background_novel

    # ---- model ---- #
    def build(self):
        import gurobipy as gp
        from gurobipy import GRB

        p = self.params
        g = self.graph
        m = gp.Model("rck-ecdna-flow")
        self.model = m
        m.Params.OutputFlag = 1 if (p.log_to_console or p.log_file) else 0
        m.Params.LogToConsole = 1 if p.log_to_console else 0
        if p.log_file:
            m.Params.LogFile = p.log_file
        if p.time_limit:
            m.Params.TimeLimit = p.time_limit
        m.Params.Threads = default_threads(p.threads)
        m.Params.MIPGap = p.mip_gap

        for sid in sorted(g.segments):
            so = self.obs.segments[sid]
            ub = int(math.ceil(max(so.total, 1.0) * p.ub_factor)) + p.ub_add
            self.ub_seg[sid] = ub
            for h in HAPLOTYPES:
                self.xe[sid, h] = m.addVar(lb=0, ub=ub, vtype=GRB.INTEGER, name="xe[{s},{h}]".format(s=sid, h=h))
                self.xb[sid, h] = m.addVar(lb=0, ub=p.bg_max, vtype=GRB.INTEGER, name="xb[{s},{h}]".format(s=sid, h=h))

        for aid in sorted(g.adjacencies):
            a = g.adjacencies[aid]
            big_m = self._big_m(aid)
            for ph in a.allowed:
                self.ye[aid, ph] = m.addVar(lb=0, ub=big_m, vtype=GRB.INTEGER, name="ye[{a},{p}]".format(a=aid, p=ph))
                if self._uses_background(aid):
                    self.yb[aid, ph] = m.addVar(lb=0, ub=p.bg_max, vtype=GRB.INTEGER, name="yb[{a},{p}]".format(a=aid, p=ph))
                if a.is_novel_like:
                    pv = m.addVar(vtype=GRB.BINARY, name="p[{a},{p}]".format(a=aid, p=ph))
                    self.pv[aid, ph] = pv
                    total = self.ye[aid, ph] + (self.yb[aid, ph] if (aid, ph) in self.yb else 0)
                    m.addConstr(total <= (big_m + p.bg_max) * pv, name="bigM[{a},{p}]".format(a=aid, p=ph))
                    m.addConstr(total >= pv, name="presence[{a},{p}]".format(a=aid, p=ph))
            if a.is_novel_like:
                m.addConstr(gp.quicksum(self.pv[aid, ph] for ph in a.allowed) <= 1, name="one-labeling[{a}]".format(a=aid))

        for nid in g.terminals():
            for h in HAPLOTYPES:
                self.t[nid, h] = m.addVar(lb=0, ub=4 * p.bg_max, vtype=GRB.INTEGER, name="t[{n},{h}]".format(n=nid, h=h))

        self.end_nodes = self.linear_end_nodes()
        for nid in self.end_nodes:
            ub = self.ub_seg[g.nodes[nid].segment]
            for h in HAPLOTYPES:
                self.te[nid, h] = m.addVar(lb=0, ub=ub, vtype=GRB.INTEGER, name="te[{n},{h}]".format(n=nid, h=h))

        self._balance(gp)
        self._objective(gp, GRB)
        if p.enforce_is:
            self._infinite_sites(gp)
        m.update()
        return m

    def _balance(self, gp):
        g, m = self.graph, self.model
        for nid in sorted(g.nodes):
            node = g.nodes[nid]
            for h in HAPLOTYPES:
                ec_terms, bg_terms = [], []
                for a in g.adjacencies_at(nid):
                    for ph in a.allowed:
                        k = a.hap_at(ph, nid).count(h)
                        if k == 0:
                            continue
                        ec_terms.append(k * self.ye[a.aid, ph])
                        if (a.aid, ph) in self.yb:
                            bg_terms.append(k * self.yb[a.aid, ph])
                ec = gp.quicksum(ec_terms)
                bg = gp.quicksum(bg_terms)
                slack = self.t.get((nid, h), 0)
                if node.segment is not None:
                    m.addConstr(self.xe[node.segment, h] == ec + self.te.get((nid, h), 0), name="ec-balance[{n},{h}]".format(n=nid, h=h))
                    m.addConstr(self.xb[node.segment, h] == bg + slack, name="bg-balance[{n},{h}]".format(n=nid, h=h))
                else:
                    m.addConstr(ec == 0, name="ec-nosegment[{n},{h}]".format(n=nid, h=h))
                    m.addConstr(bg == slack, name="bg-terminal[{n},{h}]".format(n=nid, h=h))

    def _abs_dev(self, gp, expr, target, name):
        d = self.model.addVar(lb=0, name=name)
        self.model.addConstr(d >= expr - target)
        self.model.addConstr(d >= target - expr)
        return d

    def _objective(self, gp, GRB):
        p, g, obs, m = self.params, self.graph, self.obs, self.model
        terms = []
        for sid in sorted(g.segments):
            so = obs.segments[sid]
            w = self._seg_weight(sid)
            norm = 1.0 / math.sqrt(max(so.total, 1.0))
            if so.a is not None and so.b is not None:
                for h, value in ((Haplotype.A, so.a), (Haplotype.B, so.b)):
                    d = self._abs_dev(gp, self.xe[sid, h] + self.xb[sid, h], value, "dseg[{s},{h}]".format(s=sid, h=h))
                    terms.append(w * norm * d)
            else:
                total = gp.quicksum(self.xe[sid, h] + self.xb[sid, h] for h in HAPLOTYPES)
                d = self._abs_dev(gp, total, so.total, "dseg[{s},T]".format(s=sid))
                terms.append(w * norm * d)
            bg_total = gp.quicksum(self.xb[sid, h] for h in HAPLOTYPES)
            d = self._abs_dev(gp, bg_total, obs.background.get(sid, 2.0), "dbg[{s}]".format(s=sid))
            terms.append(p.bg_prior_weight * w * d)
            if p.bg_balance_weight > 0 and so.a is not None and so.b is not None:
                d = self._abs_dev(gp, self.xb[sid, Haplotype.A] - self.xb[sid, Haplotype.B], 0.0, "dbal[{s}]".format(s=sid))
                terms.append(p.bg_balance_weight * d)

        for aid in sorted(g.adjacencies):
            a = g.adjacencies[aid]
            value = obs.adjacencies.get(aid)
            if value is not None:
                flow = gp.quicksum(self.ye[aid, ph] + (self.yb[aid, ph] if (aid, ph) in self.yb else 0) for ph in a.allowed)
                alpha = p.alpha_ref if a.kind == REFERENCE else p.alpha_novel
                d = self._abs_dev(gp, flow, value, "dadj[{a}]".format(a=aid))
                terms.append(alpha / math.sqrt(max(value, 1.0)) * d)
            if a.is_novel_like:
                cost = p.imputed_penalty if a.kind == IMPUTED else p.novel_penalty
                for ph in a.allowed:
                    extra = p.soft_junction_penalty if (a.kind != IMPUTED and ph not in a.preferred) else 0.0
                    terms.append((cost + extra) * self.pv[aid, ph])
        if self.te:
            terms.append(p.linear_end_penalty * gp.quicksum(self.te.values()))
        m.setObjective(gp.quicksum(terms), GRB.MINIMIZE)

    def _infinite_sites(self, gp):
        g, m = self.graph, self.model
        for nid in sorted(g.nodes):
            for h in HAPLOTYPES:
                terms = []
                for a in g.adjacencies_at(nid):
                    if not a.is_novel_like:
                        continue
                    for ph in a.allowed:
                        if h in a.hap_at(ph, nid):
                            terms.append(self.pv[a.aid, ph])
                if len(terms) > 1:
                    m.addConstr(gp.quicksum(terms) <= 1, name="is[{n},{h}]".format(n=nid, h=h))

    # ---- solve ---- #
    def solve(self) -> FlowSolution:
        from gurobipy import GRB

        if self.model is None:
            self.build()
        m = self.model
        m.optimize()
        names = {GRB.OPTIMAL: "optimal", GRB.TIME_LIMIT: "time_limit", GRB.SUBOPTIMAL: "suboptimal",
                 GRB.INFEASIBLE: "infeasible", GRB.INF_OR_UNBD: "inf_or_unbd", GRB.INTERRUPTED: "interrupted"}
        status = names.get(m.Status, str(m.Status))
        if m.SolCount == 0:
            self.close()
            raise RuntimeError("flow MILP returned no solution (status {s})".format(s=status))

        def val(v):
            return int(round(v.X))

        seg_ec = {k: val(v) for k, v in self.xe.items()}
        seg_bg = {k: val(v) for k, v in self.xb.items()}
        adj_ec = {k: val(v) for k, v in self.ye.items()}
        adj_bg = {k: val(v) for k, v in self.yb.items()}
        terminal = {k: val(v) for k, v in self.t.items()}
        linear_ends = {k: val(v) for k, v in self.te.items()}
        imputed = sorted({aid for (aid, ph), v in self.pv.items() if v.X > 0.5 and self.graph.adjacencies[aid].kind == IMPUTED})
        gap = m.MIPGap if m.IsMIP else 0.0
        solution = FlowSolution(status=status, objective=float(m.ObjVal), mip_gap=float(gap), seg_ec=seg_ec, seg_bg=seg_bg,
                                adj_ec=adj_ec, adj_bg=adj_bg, terminal_bg=terminal, imputed_used=imputed, linear_ends=linear_ends)
        self.close()
        return solution

    def close(self) -> None:
        """Free the Gurobi model and close its log file.

        An open log file cannot be removed cleanly on NFS (it leaves a ``.nfsXXXX`` file behind),
        so the model is disposed as soon as its values have been read.
        """
        close_model(self.model)
        self.model = None


def solve_flow(graph: ECGraph, obs: Observations, params: Optional[CircularParams] = None) -> FlowSolution:
    model = OptModelCircular(graph, obs, params)
    model.build()
    return model.solve()
