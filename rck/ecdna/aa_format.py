"""Read AmpliconArchitect (AA) graph and cycles files, from any method.

The AA formats (as ecSimulator's ``ecSim_IO`` writes them, and AA, CoRAL and ``rck-ecdna`` too):

``*_graph.txt``
    ``sequence  chr:start-  chr:end+  PredictedCopyCount  AverageCoverage  Size  NumberReadsMapped``
    and ``concordant`` / ``discordant`` / ``source`` lines
    ``kind  chr:pos±->chr:pos±  PredictedCopyCount  NumberOfReadPairs  ...``.
``*_cycles.txt``
    ``Interval  i  chrom  start  end``, ``Segment  id  chrom  start  end`` and
    ``Cycle=k;Copy_count=w;Segments=1+,2-,...`` lines. A path that is not a cycle starts and ends
    with segment ``0`` (AA's source: ``0+,...,0-``; ecSimulator writes ``0+,...,0+``).

Coordinates are 1-based and inclusive. Cycles become :class:`rck.ecdna.evaluate.Cycle` objects
with 0-based inclusive steps (``start - 1``, ``end - 1``), the convention under which
:func:`rck.ecdna.evaluate.read_aa_cycles` reads the simulator's truth, so truth and inferred
cycles read from AA files are compared in one frame.

Optional keys on a ``Cycle=`` line, which AA and CoRAL ignore:

* ``Haplotypes=A,B,...``: one haplotype per step (``rck-ecdna``; ``.`` for the source). Without
  it, the cycle is unphased and its steps carry haplotype ``.``;
* ``Status=low_weight``: the method itself does not call this cycle a species (``rck-ecdna``);
* ``Alt_segments=...|...`` and ``Alt_haplotypes=...|...``: other orderings of the same species
  (``rck-ecdna`` reports every distinct Eulerian ordering of a species' edges);
* ``Species=name``: a name for the cycle (the truth mixture's cycles files).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from rck.ecdna.evaluate import CStep, Cycle

UNPHASED = "."

Breakend = Tuple[str, int, str]  # (chrom, AA position, strand)


@dataclass
class AAEdge:
    kind: str  # concordant, discordant or source
    a: Breakend
    b: Breakend
    cn: float
    reads: Optional[float] = None


@dataclass
class AASequence:
    chrom: str
    start: int  # 1-based inclusive
    end: int
    cn: float


@dataclass
class AAGraph:
    path: str
    sequences: List[AASequence] = field(default_factory=list)
    edges: List[AAEdge] = field(default_factory=list)

    @property
    def discordant(self) -> List[AAEdge]:
        return [e for e in self.edges if e.kind == "discordant"]


@dataclass
class AACycleRecord:
    cycle_id: str
    copy_count: float
    tokens: List[Tuple[str, str]]  # (segment id, strand), the source "0" left out
    linear: bool
    haplotypes: Optional[List[str]]
    fields: Dict[str, str]
    alternatives: List[Tuple[List[Tuple[str, str]], Optional[List[str]]]] = field(default_factory=list)


@dataclass
class AACycles:
    path: str
    intervals: List[Tuple[str, int, int]]  # 1-based inclusive
    segments: Dict[str, Tuple[str, int, int]]
    cycles: List[AACycleRecord]

    @property
    def phased(self) -> bool:
        return bool(self.cycles) and all(c.haplotypes is not None for c in self.cycles)


def parse_breakend(text: str) -> Breakend:
    chrom, rest = text.strip().rsplit(":", 1)
    return chrom, int(rest[:-1]), rest[-1]


def read_graph(path: str) -> AAGraph:
    graph = AAGraph(path)
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] == "sequence" and len(f) >= 4:
                a, b = parse_breakend(f[1]), parse_breakend(f[2])
                graph.sequences.append(AASequence(a[0], a[1], b[1], float(f[3])))
            elif f[0] in ("concordant", "discordant", "source") and len(f) >= 3:
                left, right = f[1].split("->")
                reads = None
                if len(f) > 3:
                    try:
                        reads = float(f[3])
                    except ValueError:
                        reads = None
                graph.edges.append(AAEdge(f[0], parse_breakend(left), parse_breakend(right), float(f[2]), reads))
    return graph


def _tokens(text: str) -> Tuple[List[Tuple[str, str]], bool, List[int]]:
    """(steps, linear, positions of the kept steps among all tokens) of a ``Segments=`` value."""
    raw = [t.strip() for t in text.split(",") if t.strip()]
    steps, kept, linear = [], [], False
    for i, t in enumerate(raw):
        sid, strand = t[:-1], t[-1]
        if strand not in "+-":
            raise ValueError("cycle step {t} has no strand".format(t=t))
        if sid == "0":
            linear = True
            continue
        steps.append((sid, strand))
        kept.append(i)
    return steps, linear, kept


def _haplotypes(text: Optional[str], kept: Sequence[int]) -> Optional[List[str]]:
    if text is None:
        return None
    raw = [h.strip() for h in text.split(",")]
    return [raw[i] if i < len(raw) else UNPHASED for i in kept]


def read_cycles(path: str) -> AACycles:
    intervals, segments, cycles = [], {}, []
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] == "Interval" and len(f) >= 5:
                intervals.append((f[2], int(f[3]), int(f[4])))
            elif f[0] == "Segment" and len(f) >= 5:
                segments[f[1]] = (f[2], int(f[3]), int(f[4]))
            elif line.startswith("Cycle="):
                parts = dict(p.split("=", 1) for p in line.strip().split(";") if "=" in p)
                steps, linear, kept = _tokens(parts.get("Segments", ""))
                haps = _haplotypes(parts.get("Haplotypes"), kept)
                record = AACycleRecord(cycle_id=parts["Cycle"], copy_count=float(parts.get("Copy_count", "nan")), tokens=steps,
                                       linear=linear, haplotypes=haps, fields=parts)
                if parts.get("Alt_segments"):
                    alt_haps = parts.get("Alt_haplotypes", "").split("|") if parts.get("Alt_haplotypes") else []
                    for k, text in enumerate(parts["Alt_segments"].split("|")):
                        alt_steps, _, alt_kept = _tokens(text)
                        record.alternatives.append((alt_steps, _haplotypes(alt_haps[k], alt_kept) if k < len(alt_haps) else None))
                for sid, _ in steps:
                    if sid not in segments:
                        raise ValueError("{p}: cycle {c} uses segment {s}, which is not listed".format(p=path, c=record.cycle_id, s=sid))
                cycles.append(record)
    return AACycles(path, intervals, segments, cycles)


def _steps(cycles: AACycles, tokens, haps) -> List[CStep]:
    out = []
    for i, (sid, strand) in enumerate(tokens):
        chrom, start, end = cycles.segments[sid]
        out.append(CStep(chrom, start - 1, end - 1, haps[i] if haps is not None else UNPHASED, strand))
    return out


def to_cycles(cycles: AACycles, prefix: str = "") -> List[Cycle]:
    """Every cycle (and path) of a cycles file as a scoring :class:`Cycle`, named ``<prefix>``
    + ``Species=`` when given, else ``<prefix>cycle<k>``."""
    out = []
    for rec in cycles.cycles:
        if not rec.tokens:
            continue
        orderings = [_steps(cycles, rec.tokens, rec.haplotypes)]
        orderings += [_steps(cycles, t, h if h is not None else rec.haplotypes) for t, h in rec.alternatives]
        name = prefix + (rec.fields.get("Species") or "cycle{k}".format(k=rec.cycle_id))
        out.append(Cycle(name=name, weight=rec.copy_count, orderings=orderings, circular=not rec.linear))
    return out


def raw_walks(cycles: AACycles) -> List[Tuple[List[dict], List[Tuple[int, str]]]]:
    """[(segment rows, walk)] in the layout of :func:`rck.ecdna.evaluate.read_aa_cycles` (which
    reads the simulator's truth this way): rows ``{#CHR, START = start - 1, END = end}`` in order
    of first use, walk as (1-based row index, strand)."""
    result = []
    for rec in cycles.cycles:
        if rec.linear:
            raise ValueError("{p}: cycle {c} is a path (segment 0); truth species must be cycles".format(p=cycles.path, c=rec.cycle_id))
        rows, walk, index = [], [], {}
        for sid, strand in rec.tokens:
            if sid not in index:
                chrom, start, end = cycles.segments[sid]
                rows.append({"#CHR": chrom, "START": start - 1, "END": end})
                index[sid] = len(rows)
            walk.append((index[sid], strand))
        if walk:
            result.append((rows, walk))
    return result
