# Graph route (Task 2): every decomposer on the same breakpoint graph

Every decomposer (CoRAL 3.1.1 `cycle`, CycleExtractor 1.0.0, original RCK and RCK-ecDNA) gets the
same AA-format breakpoint graph, so differences come from the decomposition alone. Two graphs
per case: the **truth mixture graph** (oracle) and our **inferred union subgraph** (realistic),
for demo5–demo8 and the 32 bench cases. The plan and the state of the work are in
[NEXT_STEPS.md](NEXT_STEPS.md); environments in [envs/README.md](envs/README.md).

| file | what |
|---|---|
| `envs/` | lock files of `/n/fs/yubi-iw/conda/envs/coral-3.1.1` (CoRAL and CycleExtractor both run with it) |
| `fixtures/` | graphs written by CoRAL 3.1.1's own writer (`coral_rewrite.py`): the reference layout |
| `coral_rewrite.py` | parse a graph with CoRAL and write it back with CoRAL's writer |
| `check_coral_format.py` | check AA graphs against CoRAL 3.1.1's parser, sorter and writer; `--golden`; `--ce` (CycleExtractor, HiGHS) |
| `roundtrip.py`, `roundtrip.sbatch`, `roundtrip.list` | the round-trip test (subgraph → AA → subgraph → RCK-ecDNA total-only) |
| `roundtrip_prepare.tsv` | its offline part: per case, the graph comparison of each arm |
| `graph_route.py`, `graph_route.sbatch`, `graph_route.list` | the head-to-head runs (priority 1 a): CoRAL, RCK-ecDNA (`aa_import --imputed none`) and CycleExtractor on the truth `graph_reads_paths` and inferred graphs, one run per line, `run.json` per run |
| `graph_route_prepare.tsv` | per graph: input file and the import's counts and fitted scale |
| `$G/<case>/runs/<graph>/coral_input/` | the graph CoRAL runs on: the input graph with 10 kb flanks at every interval end (`rck.ecdna.coral_flank`; CoRAL cannot use a breakpoint on an interval's end node) |
| `native_aa_metrics.tsv` | the native arm: RCK-ecDNA's bench runs at λ = 0.2, scored by aa_eval (rows of `results/bench_nested/aa_metrics.tsv`) |
| `score_h2h.py` | scores every run and the native arm with aa_eval (offline): `h2h_table.tsv` (case x graph x method x mode), `h2h_scenario.tsv`, `h2h_species.tsv` (per truth species), `h2h_graph_metrics.tsv` (once per input graph); `--drop-paths` writes the `_cycles_only` sensitivity tables |

## Head-to-head results (2026-11-06, single seed)

**RCK-ecDNA's numbers are single-seed.** In the round trip, its decomposition changed in 29 of 36
cases with an unchanged flow objective; seed repeats are pending (NEXT_STEPS step e). Its per-case
scores swing more than the other methods' between truth and inferred graphs. It scores 0.5–0.6
*lower* on the cleaner truth graph in nest_cross_1to1_s509, nest_cross_2to1_s103 and
nest_samehap_3to1_s304, which fits tie-breaking. Read its case-level numbers with that in mind.

All 225 runs finished (`graph_route.py status`: 225 done, none failed). Every CoRAL run has
≥ 1 cycle. Scored by `score_h2h.py` with aa_eval's defaults (strict: CoRAL / CE `cn-frac:0.1`,
RCK-ecDNA `status`; best-of: every cycle, plus RCK-ecDNA's 2 same-fit alternatives), on the 67
truth species of the 36 cases, all circular. Values are the mean composite-H0 over cases,
strict / best-of. RCK-ecDNA* = total-only on the shared graph. Native = its bench run, which has
allele data, its own subgraph and imputed edges; the re-score reproduces `native_aa_metrics.tsv`.

| scenario (cases, species) | truth CoRAL | truth CE | truth RCK-ecDNA* | inferred CoRAL | inferred CE | inferred RCK-ecDNA* | native |
|---|---|---|---|---|---|---|---|
| demo (4, 7) | 0.99 / 0.99 | 0.79 / 0.85 | 0.93 / 0.93 | 0.98 / 0.98 | 0.59 / 0.69 | 0.98 / 0.98 | 0.99 / 0.99 |
| ctrl_dup (6, 6) | 1.00 / 1.00 | 1.00 / 1.00 | 0.84 / 0.84 | 0.96 / 0.96 | 0.96 / 0.96 | 0.82 / 0.82 | 0.82 / 0.87 |
| ctrl_separate (2, 4) | 1.00 / 1.00 | 1.00 / 1.00 | 1.00 / 1.00 | 0.95 / 0.95 | 0.82 / 0.95 | 0.97 / 0.97 | 0.97 / 0.97 |
| nest_cross_1to1 (6, 12) | 0.35 / 0.52 | 0.57 / 0.60 | 0.63 / 0.69 | 0.36 / 0.49 | 0.39 / 0.48 | 0.48 / 0.56 | 0.86 / 0.86 |
| nest_cross_2to1 (2, 4) | 0.65 / 0.71 | 0.67 / 0.71 | 0.71 / 0.72 | 0.66 / 0.67 | 0.66 / 0.67 | 0.83 / 0.83 | 0.96 / 0.96 |
| nest_cross_3sp (2, 6) | 0.26 / 0.37 | 0.30 / 0.43 | 0.37 / 0.62 | 0.31 / 0.36 | 0.29 / 0.37 | 0.36 / 0.50 | 0.70 / 0.74 |
| nest_samehap_2to1 (2, 4) | 1.00 / 1.00 | 0.60 / 0.71 | 1.00 / 1.00 | 0.66 / 0.66 | 0.29 / 0.39 | 0.69 / 0.84 | 0.99 / 0.99 |
| nest_samehap_3to1 (6, 12) | 0.82 / 0.82 | 0.90 / 0.90 | 0.79 / 0.87 | 0.62 / 0.64 | 0.55 / 0.60 | 0.65 / 0.75 | 0.82 / 0.96 |
| nest_samehap_6to1 (6, 12) | 0.65 / 0.76 | 0.92 / 0.92 | 0.67 / 0.80 | 0.55 / 0.67 | 0.63 / 0.71 | 0.60 / 0.75 | 0.55 / 0.86 |
| **all (36, 67)** | **0.74 / 0.80** | **0.80 / 0.82** | **0.76 / 0.82** | **0.67 / 0.71** | **0.60 / 0.67** | **0.69 / 0.76** | **0.82 / 0.91** |

Species counts over the 67 (strict; exact = AA/CoRAL exact reconstruction at 100 bp, haplotype-blind):

| graph | method | exact (strict / best-of) | missed | merged | missed or merged | false (circular + linear) |
|---|---|---|---|---|---|---|
| truth | CoRAL | 39 / 39 | 5 | 20 | 24 | 2 + 10 |
| truth | CE | 43 / 43 | 5 | 24 | 25 | 1 + 11 |
| truth | RCK-ecDNA* | 37 / 41 | 3 | 30 | 32 | 7 + 0 |
| inferred | CoRAL | 32 / 32 | 5 | 28 | 32 | 10 + 8 |
| inferred | CE | 24 / 25 | 4 | 33 | 36 | 17 + 72 |
| inferred | RCK-ecDNA* | 33 / 36 | 4 | 31 | 35 | 16 + 3 |
| native | RCK-ecDNA | 50 / 56 | 6 | 11 | 17 | 8 + 2 |

**Main results.** On the shared total-only graph, no method is clearly ahead. On truth graphs, CE
leads on strict composite-H0 (0.80, 43 exact) and RCK-ecDNA and CE tie on best-of (0.82). On the
inferred graphs RCK-ecDNA leads (0.69 / 0.76, against CoRAL 0.67 / 0.71 and CE 0.60 / 0.67).
These margins are mostly within RCK-ecDNA's tie-breaking noise. What clearly separates is allele
data: the native arm (0.82 / 0.91, 50 exact) beats every total-only arm, most in the nest_cross
scenarios (0.86 against at most 0.63), where nested species on different haplotypes cannot be
separated from total copy number alone. By scenario:
- CoRAL and CE solve the fold-back controls (ctrl_dup, 1.00 on truth). There RCK-ecDNA splits the
  truth walk into two cycles at an identical copy number (s708, s710; the same in the native arm).
- CE wins nest_samehap_3to1 and 6to1 on truth graphs (0.90, 0.92), where its first greedy cycle
  satisfies the path constraints. The same greedy step makes one long walk that merges both demo8
  species at the minor CN, and on inferred graphs it emits many linear walks (72 false). CE loses
  the most from truth to inferred graph (0.80 → 0.60).
- CoRAL is weakest on nest_cross (0.35 strict on 1to1) and, on truth graphs, on samehap_6to1.

Every arm loses 0.05–0.20 from truth to inferred graph. The inferred graphs keep every breakpoint
(`h2h_graph_metrics.tsv`: `bp_accuracy` 1.0, no false edge, 38 of 39 `error_free_merged`), so the
loss comes from segmentation and copy number, not missing junctions.

**What depends on the CoRAL flank fix and on CoRAL's multiplicity bound.**
- Every CoRAL number depends on `coral_flank`: without it, CoRAL returned 0 cycles.
- No CoRAL *cycle* uses a flank segment. Flanks appear only in CoRAL's linear s/t paths (22 strict
  paths on truth graphs, 15 on inferred). Of these, 12 and 7 are the one-to-one match of a truth
  species (partial credit, never exact).
- `aa_format.read_cycles` skips CoRAL/CE `Path=` lines, so the headline reads them as AA-style
  linear walks, as RCK-ecDNA's are scored. Dropping every linear walk (`h2h_*_cycles_only.tsv`)
  gives CoRAL truth strict 0.71 (missed 5 → 15) and CE 0.79; RCK-ecDNA is unchanged. No ranking
  changes.
- nest_samehap_6to1_s804: the CN-40.8 species exceeds CoRAL's read-clustered discordant
  multiplicity bound, and CoRAL's heaviest cycle is CN 6.8. CoRAL scores 0.39 truth / 0.35
  inferred (CE 0.54 / 0.52, RCK-ecDNA 0.72 / 0.71). Leaving s804 out raises CoRAL's overall
  truth strict mean by 0.01 (0.742 → 0.752).

**Not computed.** `subwalks_satisfied`: `aa_eval.read_paths` needs `chr:pos` ids, and the inferred
graphs' `paths.tsv` carry the subgraph's `imputed_ext_*` ids. CoRAL's own
`Path_constraints_satisfied` is in its cycles files.

Outputs go to `/n/fs/yubi-iw/runs/graph_route/` (project space, not purged; rebuild with the commands below).

## Commands

```bash
cd ~/RCK
PY=/n/fs/yubi-iw/conda/envs/gurobi-env/bin/python           # RCK
CPY=/n/fs/yubi-iw/conda/envs/coral-3.1.1/bin/python         # CoRAL / CycleExtractor
G=/n/fs/yubi-iw/runs/graph_route

# inferred subgraph -> AA graph (+ _segment_cn.tsv, _export.json)
$PY -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR $G/CASE/inferred/graph --name CASE [--imputed comment|drop]

# truth mixture graph, optionally with read counts and long-read path constraints (amplicon-allele-analysis)
$PY ~/amplicon-allele-analysis/simulations/truth_mixture.py RUN_DIR --outdir $G/CASE/truth/graph \
    [--reads-per-copy 6.5 --depth-per-copy D] [--path-constraints]

# AA graph -> subgraph folder for rck-ecdna --inputs subgraph
$PY -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate] [--cn-profile F] [--background-cn X] [--compare ORIGINAL]

# check graphs with CoRAL itself (no solver), and CycleExtractor with HiGHS (no licence)
PYTHONPATH=~/RCK $CPY results/graph_route/check_coral_format.py GRAPH... [--paths-from-export]
PYTHONPATH=~/RCK $CPY results/graph_route/check_coral_format.py --ce OUTDIR GRAPH...
```

## What the AA graph holds

`aa_export subgraph` writes the graph RCK-ecDNA solves on (`read_subgraph`: phantoms contracted,
parallel imputed edges merged) with the copy numbers it uses: segment CN = `CN_TOTAL`, junction CN
= reads / junction reads per copy (RCK-ecDNA observes junctions this way, not through the rounded
`CN` column of edges.tsv). Positions are subgraph coordinates + 1, as for exported answers.

| line | content | seen by CoRAL / CE |
|---|---|---|
| `sequence` | one per segment, CoRAL's sort order; AverageCoverage = depth | yes |
| `concordant` | reference joins between segments; an unobserved join (imputed grid point) gets the smaller segment CN | yes |
| `discordant` | observed junctions; fold-backs (both ends on one node) as self-loops | yes |
| `path_constraint` | read paths of paths.tsv over segments, reference joins and observed junctions, CoRAL 3.1.1 tokens, support = `#NOV` | yes |
| `interval` | runs of contiguous segments (what CoRAL infers itself) | ignored |
| `#rck_imputed` | RCK-ecDNA's imputed candidates (no read evidence) | no |
| `#rck_outside` | joins to a region's outer flank (a node without a segment); CoRAL and CE add their own s/t ends instead | no |
| `#rck_node`, `#rck_seq`, `#rck_edge`, `#rck_pc`, `#rck_path` | the subgraph's node / segment / adjacency / path ids and orders, so the import is exact | no |

`aa_import` inverts it. `--imputed comment` restores the `#rck_imputed` lines, `none` leaves them out
(what CoRAL and CE see), `regenerate` proposes RCK-ecDNA's own closing candidates on the graph
(subgraph_construct.py's rule adapted to sequence edges; NEXT_STEPS.md). `#rck_outside` lines are
always restored. `scale.json` is fitted from the graph (junction scale = median reads / CN; a graph
without reads gets WEIGHT = CN and scale 1; no depth fit from a CoRAL-rewritten graph).

## Validation done offline (2026-11-03, no Gurobi)

* **Round trip, 36/36 cases** (demo5–demo8 + 32 bench; `roundtrip_prepare.tsv`): `--imputed comment`
  with the original cn_profile gives back the normalised graph exactly (node, segment and adjacency
  ids and orders, reads, observations, junction scale to the bit, read paths); `--imputed none` is
  exact apart from the imputed candidates; graph-only `none` differs only in the flank background
  (2 instead of the cn_profile). Unit tests: `tests/ecdna/test_aa_import.py`.
* **CoRAL 3.1.1, 36 exported + 78 truth graphs** (`check_coral_format.py`): parse counts equal ours,
  `infer_amplicon_intervals` gives our intervals and keeps the edge order, CoRAL's writer round trip
  keeps edges, CN and constraints; all 418 path constraints of the exported graphs are rebuilt
  independently on CoRAL's own adjacencies and written by its `path_to_str` with our exact tokens;
  every truth constraint crosses exactly the junctions of its read.
* **CycleExtractor (HiGHS)** on demo5, demo7, ctrl_dup_s702 (exported) and ctrl_dup_s710 (truth
  with reads and constraints): parses, writes cycles `aa_format.read_cycles` reads; on the
  ctrl_dup_s710 truth graph it returns the truth cycle exactly, both fold-backs included.

## Fold-back junctions

A fold-back (inverted duplication whose two breakends are the same node) is a discordant self-loop.
14 of the 32 bench graphs have one: ctrl_dup_s702, s703, s708, s709, s710; nest_cross_1to1_s505,
s506, s508; nest_cross_3sp_s402, s403; nest_samehap_2to1_s201; nest_samehap_6to1_s803, s804, s805
(the ids are in each `_export.json` "foldbacks").

* **CycleExtractor 1.0.0 does not drop them** (correcting the earlier note in NEXT_STEPS.md):
  `Create_graph` replaces each self-loop by a gadget, `node -> i`, a 1 bp sequence edge `i -> i+1`,
  `i+1 -> node` (the `if start != end` branch only routes loops to it). Seen working: on
  ctrl_dup_s702 its cycle turns at the fold-back (`...5+,6+` then `6-,5-...`) and satisfies the
  fold-back path constraint `e6+:2,d4+:1,e6-:2`; on the ctrl_dup_s710 truth graph it recovers the
  truth cycle with both fold-backs. What remains to watch: the gadget's two half-edges each carry
  the junction's read count into CE's multiplicity bound.
* CoRAL models self-loops directly: `add_discordant_edge` lists a self-loop twice at its node, so
  the cycle model's Eulerian balance counts one traversal twice there (and its CN fit,
  `DiscordantEdge.is_self_loop`, treats it apart).
* RCK-ecDNA: a self-loop's phasing is same-haplotype only, as in RCK.

## Known limitations (recorded, not worked around)

* **CoRAL ignores breakpoints on interval end nodes** (worked around for CoRAL only). The two outer
  nodes of each run of contiguous sequence edges are CoRAL's end nodes, balanced only by its s/t
  edges (`get_eulerian_node_constraints`), so a cycle through a breakpoint there is infeasible.
  Our graphs end at breakpoints (69 of 75), so `graph_route.py` gives CoRAL the graph with 10 kb
  flanks (`python -m rck.ecdna.coral_flank GRAPH OUT [--flank 10000]`); NEXT_STEPS.md has the checks.

* **CoRAL re-sorts without renumbering.** `coral cycle` calls `infer_amplicon_intervals`, which sorts
  sequence and concordant edges (natural chromosome order) after the path constraints were parsed by
  file index. The exporter and `truth_mixture.py --path-constraints` write CoRAL's order; a graph
  from elsewhere with constraints in another order is silently corrupted by CoRAL
  (ctrl_separate_s603, chr15 + chr7, would have been).
* **CoRAL's writer** puts `lr_nc / size` in AverageCoverage and truncates read counts to integers;
  never fit a depth scale from a CoRAL-rewritten graph (`aa_import` refuses to).
* **Direction of a path constraint**: CoRAL writes a path in canonical direction; RCK-ecDNA matches
  paths in both directions, and `#rck_pc` restores the read's own direction anyway.
* **What CoRAL / CE cannot see**: allele copy numbers, phase blocks and junction haplotype labels
  (the import is total-only; `_segment_cn.tsv` keeps alleles for scoring); RCK-ecDNA's flank joins
  (they add their own s/t ends; the import restores the joins even with `--imputed none`, as they
  are RCK-ecDNA's counterpart of those ends).
* **Unobserved concordant joins** (at imputed grid points, 10 in demo7) get the smaller neighbouring
  segment CN in the AA graph; RCK-ecDNA has no observation there.
* **Truth vs inferred read paths**: the truth graph gets every read of `obs/path.tsv` whose junctions
  are mixture edges (demo7: 12); the subgraph's paths.tsv was filtered by `rck_inputs.clip_path`
  (demo7: 3). Synthetic truth reads are `round(CN x R)`, so the fitted junction scale is close to,
  not exactly, R (demo7: 6.497 for 6.5).
