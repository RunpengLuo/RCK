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
* **Depth bias (optional, off by default).** With ``depth_bias_window`` > 0, segments are grouped
  into windows of about that many bp of nearby sequence (:func:`depth_bias_windows`), and each
  window k gets a continuous bias beta_k in [-``depth_bias_max``, ``depth_bias_max``]. The
  segment (and phase-block piece) targets become obs x (1 - beta_k), which stays linear because obs
  is data; a flippable piece's target uses one exact product variable (flip x beta). The objective
  adds ``depth_bias_cost`` x |beta_k| and ``depth_bias_smooth`` x |beta_k - beta_j| over neighbouring
  windows. This absorbs a smooth multiplicative coverage bias (GC, mappability) shared by nearby
  segments, so the integer flow need not follow a few-copy dip that no junction supports; junction
  reads, which the depth bias does not reach, are left as they are.
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
    depth_bias_window: int = 0  # bp; > 0 fits one multiplicative depth bias per window of nearby segments; 0 = off
    depth_bias_cost: float = 1.0  # per unit |bias| (a 10 % shift of one window costs 0.1 x this)
    depth_bias_smooth: float = 1.0  # per unit |bias_k - bias_j| between neighbouring windows
    depth_bias_max: float = 0.3  # bound on |bias|
    enforce_is: bool = False
    allow_background_novel: bool = True  # chromosomal SVs may carry background flow
    linear_ends: str = "boundaries"  # boundaries | all | none: where ecDNA-layer copies may start / stop
    linear_end_penalty: float = 0.25  # per copy per linear end
    bg_balance_weight: float = 0.01  # tie-break: background copies split evenly between haplotypes
    boundary_min_copies: float = 2.0  # copy-number step that makes an extremity a boundary ...
    boundary_frac: float = 0.25  # ... or this fraction of the larger copy number, whichever is larger
    flip_penalty: float = 1e-6  # tie-break per flipped unoriented phase block (--block-flip): keep the observed frame
    flip_oriented_penalty: Optional[float] = None  # cost of flipping a block the subgraph oriented; None = never
    flip_tol: float = 0.01  # a block whose reversed flip changes the objective by at most this is "unresolved"
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
    block_flips: Dict[str, int] = field(default_factory=dict)  # phase block -> 1 if its A / B were swapped (--block-flip)
    blocks: List[Dict[str, object]] = field(default_factory=list)  # one block_flips.tsv row per phase block
    depth_bias: Dict[str, float] = field(default_factory=dict)  # window -> fitted bias (depth_bias_window > 0 only)
    seg_bias: Dict[str, float] = field(default_factory=dict)  # segment -> its window's bias (empty when off)

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
            obs_a, obs_b = obs.segment_alleles(sid, self.block_flips)
            rows.append({
                "sid": sid, "chrom": s.chrom, "start": s.start, "end": s.end, "length": s.length,
                "obs_total": so.total, "obs_A": obs_a, "obs_B": obs_b, "obs_source": so.source,
                "ecdna_A": self.seg_ec[(sid, Haplotype.A)], "ecdna_B": self.seg_ec[(sid, Haplotype.B)],
                "bg_A": self.seg_bg[(sid, Haplotype.A)], "bg_B": self.seg_bg[(sid, Haplotype.B)],
            })
        frame = pd.DataFrame(rows)
        if len(frame):
            frame["inferred_total"] = frame[["ecdna_A", "ecdna_B", "bg_A", "bg_B"]].sum(axis=1)
            if self.seg_bias:  # only with --depth-bias-window, so default outputs keep their columns
                frame["depth_bias"] = [round(self.seg_bias.get(sid, 0.0), 6) for sid in frame["sid"]]
                frame["obs_corrected"] = frame["obs_total"] * (1.0 - frame["depth_bias"])
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


def depth_bias_windows(graph: ECGraph, window: int) -> Tuple[Dict[str, str], List[str], List[Tuple[str, str]]]:
    """Group segments into depth-bias windows: (segment -> window, windows in order, neighbouring
    window pairs).

    Per chromosome, segments in coordinate order are added to the current window until the next one
    starts ``window`` bp or more after the window's first segment, or a gap of more than ``window``
    bp separates them; a segment longer than ``window`` is a window of its own. Consecutive windows
    are neighbours (smoothness term) unless such a gap separates them: distant parts of an amplicon do
    not share a coverage bias. Windows are named ``<chrom>:<n>``."""
    if window <= 0:
        return {}, [], []
    by_chrom: Dict[str, List] = {}
    for sid, s in graph.segments.items():
        by_chrom.setdefault(s.chrom, []).append((s.start, s.end, sid))
    seg_window, windows, pairs = {}, [], []
    for chrom in sorted(by_chrom):
        current, first_start, prev_end = None, None, None
        for start, end, sid in sorted(by_chrom[chrom]):
            gap = None if prev_end is None else start - prev_end - 1
            if current is None or gap > window or start - first_start >= window or end - first_start + 1 > window:
                name = "{c}:{n}".format(c=chrom, n=sum(1 for w in windows if w.startswith(chrom + ":")))
                if current is not None and gap <= window:
                    pairs.append((current, name))
                current, first_start = name, start
                windows.append(name)
            seg_window[sid] = current
            prev_end = end if prev_end is None else max(prev_end, end)
    return seg_window, windows, pairs


class OptModelCircular:
    """Build and solve the two-layer flow MILP for one amplicon graph."""

    def __init__(self, graph: ECGraph, obs: Observations, params: Optional[CircularParams] = None):
        self.graph = graph
        self.obs = obs
        self.params = params or CircularParams()
        self.model = None
        self.xe, self.xb, self.ye, self.yb, self.pv, self.t, self.te = {}, {}, {}, {}, {}, {}, {}
        self.flip = {}  # phase block -> binary A / B swap (--block-flip)
        self.beta = {}  # depth-bias window -> continuous bias (--depth-bias-window)
        self.flip_beta = {}  # (phase block, window) -> flip x bias (exact linearisation)
        self.seg_window: Dict[str, str] = {}
        self.bias_pairs: List[Tuple[str, str]] = []
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

    def flippable_blocks(self) -> set:
        """Phase blocks that get a flip variable: those the subgraph left unoriented, plus the
        oriented ones when ``flip_oriented_penalty`` is set; only blocks with allele data."""
        return {pc.block for pieces in self.obs.pieces.values() for pc in pieces
                if pc.has_alleles and (not pc.fixed or self.params.flip_oriented_penalty is not None)}

    def _flip_cost(self, block: str) -> float:
        fixed = any(pc.fixed for ps in self.obs.pieces.values() for pc in ps if pc.block == block)
        return self.params.flip_oriented_penalty if fixed else self.params.flip_penalty

    def _piece_terms(self, gp, sid: str, w: float, norm: float) -> List:
        """Segment data term split over its phase-block pieces by length. Allele targets are
        linear in the block's flip f: t_A = a + f (b - a), t_B = b + f (a - b)."""
        terms = []
        beta = self._beta(sid)
        for i, pc in enumerate(self.obs.pieces[sid]):
            weight = w * pc.share * norm
            if pc.has_alleles:
                f = self.flip.get(pc.block)
                for h, (own, other) in ((Haplotype.A, (pc.a, pc.b)), (Haplotype.B, (pc.b, pc.a))):
                    if f is None:
                        target = self._biased(own, beta)
                    elif beta is None:
                        target = own + (other - own) * f
                    else:  # (own + (other - own) f)(1 - beta) = own (1 - beta) + (other - own)(f - f beta)
                        target = own - own * beta + (other - own) * (f - self._flip_beta(gp, pc.block, self.seg_window[sid]))
                    d = self._abs_dev(gp, self.xe[sid, h] + self.xb[sid, h], target, "dpc[{s},{i},{h}]".format(s=sid, i=i, h=h))
                    terms.append(weight * d)
            else:
                total = gp.quicksum(self.xe[sid, h] + self.xb[sid, h] for h in HAPLOTYPES)
                d = self._abs_dev(gp, total, self._biased(pc.total, beta), "dpc[{s},{i},T]".format(s=sid, i=i))
                terms.append(weight * d)
        return terms

    # ---- depth bias (--depth-bias-window) ---- #
    def _beta(self, sid: str):
        """The bias variable of the segment's window, or None when the option is off."""
        wid = self.seg_window.get(sid)
        return self.beta.get(wid) if wid is not None else None

    @staticmethod
    def _biased(value: float, beta):
        """value x (1 - beta): the observation with the window's bias taken out (value when off)."""
        return value if beta is None else value - value * beta

    def _flip_beta(self, gp, block: str, wid: str):
        """z = flip x beta, exact for a binary flip and |beta| <= depth_bias_max (McCormick)."""
        key = (block, wid)
        if key not in self.flip_beta:
            m, f, b, bound = self.model, self.flip[block], self.beta[wid], self.params.depth_bias_max
            z = m.addVar(lb=-bound, ub=bound, name="fbeta[{b},{w}]".format(b=block, w=wid))
            m.addConstr(z <= bound * f)
            m.addConstr(z >= -bound * f)
            m.addConstr(z <= b + bound * (1 - f))
            m.addConstr(z >= b - bound * (1 - f))
            self.flip_beta[key] = z
        return self.flip_beta[key]

    def _bias_terms(self, gp) -> List:
        p = self.params
        terms = []
        for wid, b in sorted(self.beta.items()):
            terms.append(p.depth_bias_cost * self._abs_dev(gp, b, 0.0, "dbeta[{w}]".format(w=wid)))
        for u, v in self.bias_pairs:
            terms.append(p.depth_bias_smooth * self._abs_dev(gp, self.beta[u] - self.beta[v], 0.0, "dsmooth[{u},{v}]".format(u=u, v=v)))
        return terms

    def block_report(self, seg_cn: Dict[Tuple[str, Haplotype], int], flips: Dict[str, int]) -> Tuple[Dict[str, int], List[Dict[str, object]]]:
        """Final flips and one row per phase block, from the solved copy numbers.

        ``delta_if_reversed`` is the objective change if the block's flip were reversed with every
        copy number fixed. A not-oriented block within ``flip_tol`` of that is ``unresolved`` and
        is kept as observed (its flip set back to 0)."""
        p, obs = self.params, self.obs
        rows: Dict[str, Dict[str, object]] = {}
        misfit: Dict[str, List[float]] = {}  # block -> [data cost unflipped, flipped]
        for sid, pieces in obs.pieces.items():
            norm = 1.0 / math.sqrt(max(obs.segments[sid].total, 1.0))
            w = self._seg_weight(sid)
            beta = self._beta(sid)
            scale = 1.0 - (float(beta.X) if beta is not None else 0.0)  # the targets the flow fitted (--depth-bias-window)
            for pc in pieces:
                row = rows.setdefault(pc.block, {"phase_block": pc.block, "orientation": pc.orientation, "n_pieces": 0, "n_snps": 0})
                row["n_pieces"] += 1
                row["n_snps"] += pc.n_snps
                cost = misfit.setdefault(pc.block, [0.0, 0.0])
                if not pc.has_alleles:
                    continue
                for k, flip in enumerate((0, 1)):
                    ta, tb = (scale * x for x in pc.alleles(flip))
                    cost[k] += w * pc.share * norm * (abs(seg_cn[(sid, Haplotype.A)] - ta) + abs(seg_cn[(sid, Haplotype.B)] - tb))
        final = {}
        for block in sorted(rows):
            row, (c0, c1) = rows[block], misfit[block]
            flip = int(flips.get(block, 0))
            variable = block in self.flip
            penalty = self._flip_cost(block) if variable else 0.0
            delta = (c0 - c1 - penalty) if flip else (c1 + penalty - c0)
            if row["orientation"] != 0 and not flip:
                status = "fixed"
            elif row["orientation"] == 0 and (not variable or delta <= p.flip_tol):
                status = "unresolved"
                if flip:
                    log.info("phase block %s: flip changes the objective by only %.4g; kept as observed", block, delta)
                    flip, delta = 0, -delta
            else:
                status = "flipped" if flip else "kept"
            row.update({"flip": flip, "status": status, "delta_if_reversed": round(delta, 6)})
            final[block] = flip
        return final, [rows[b] for b in sorted(rows, key=_block_order)]

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

        for block in sorted(self.flippable_blocks()):
            self.flip[block] = m.addVar(vtype=GRB.BINARY, name="flip[{b}]".format(b=block))

        if p.depth_bias_window > 0:
            self.seg_window, windows, self.bias_pairs = depth_bias_windows(g, p.depth_bias_window)
            for wid in windows:
                self.beta[wid] = m.addVar(lb=-p.depth_bias_max, ub=p.depth_bias_max, name="beta[{w}]".format(w=wid))

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
            if obs.pieces.get(sid):
                terms.extend(self._piece_terms(gp, sid, w, norm))
            elif so.a is not None and so.b is not None:
                for h, value in ((Haplotype.A, so.a), (Haplotype.B, so.b)):
                    d = self._abs_dev(gp, self.xe[sid, h] + self.xb[sid, h], self._biased(value, self._beta(sid)),
                                      "dseg[{s},{h}]".format(s=sid, h=h))
                    terms.append(w * norm * d)
            else:
                total = gp.quicksum(self.xe[sid, h] + self.xb[sid, h] for h in HAPLOTYPES)
                d = self._abs_dev(gp, total, self._biased(so.total, self._beta(sid)), "dseg[{s},T]".format(s=sid))
                terms.append(w * norm * d)
            bg_total = gp.quicksum(self.xb[sid, h] for h in HAPLOTYPES)
            d = self._abs_dev(gp, bg_total, obs.background.get(sid, 2.0), "dbg[{s}]".format(s=sid))
            terms.append(p.bg_prior_weight * w * d)
            if p.bg_balance_weight > 0 and obs.has_alleles(sid):
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
        if self.beta:
            terms.extend(self._bias_terms(gp))
        for block, f in sorted(self.flip.items()):
            terms.append(self._flip_cost(block) * f)
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
        flips, blocks = {}, []
        if self.obs.pieces:
            seg_cn = {k: seg_ec[k] + seg_bg[k] for k in seg_ec}
            flips, blocks = self.block_report(seg_cn, {b: val(v) for b, v in self.flip.items()})
        depth_bias = {wid: float(b.X) for wid, b in self.beta.items()}
        seg_bias = {sid: depth_bias[wid] for sid, wid in self.seg_window.items()} if depth_bias else {}
        solution = FlowSolution(status=status, objective=float(m.ObjVal), mip_gap=float(gap), seg_ec=seg_ec, seg_bg=seg_bg,
                                adj_ec=adj_ec, adj_bg=adj_bg, terminal_bg=terminal, imputed_used=imputed, linear_ends=linear_ends,
                                block_flips=flips, blocks=blocks, depth_bias=depth_bias, seg_bias=seg_bias)
        self.close()
        return solution

    def close(self) -> None:
        """Free the Gurobi model and close its log file.

        An open log file cannot be removed cleanly on NFS (it leaves a ``.nfsXXXX`` file behind),
        so the model is disposed as soon as its values have been read.
        """
        close_model(self.model)
        self.model = None


def _block_order(block: str):
    chrom, _, idx = block.partition(":")
    return chrom, int(idx) if idx.isdigit() else -1


def solve_flow(graph: ECGraph, obs: Observations, params: Optional[CircularParams] = None) -> FlowSolution:
    model = OptModelCircular(graph, obs, params)
    model.build()
    return model.solve()
