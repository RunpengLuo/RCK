# Graph route (Task 2): every decomposer on the same breakpoint graph

The graph route gives every decomposer the same AmpliconArchitect-format (AA) breakpoint graph,
so that differences in the answers come from the decomposition alone and not from how each tool
builds its graph. The decomposers are:

- CoRAL 3.1.1 `cycle` [Zhu et al. 2024](../../docs/pipeline/references.md#zhu2024);
- CycleExtractor 1.0.0 [Faizrahnemoon et al. 2026](../../docs/pipeline/references.md#cycleextractor);
- RCK-ecDNA (this pipeline's `rck-ecdna`);
- later, the original RCK [Aganezov & Raphael 2020](../../docs/pipeline/references.md#aganezov2020).

Each case gets two graphs:

- the **truth mixture graph** (an oracle input);
- our **inferred union subgraph** (a realistic input).

There are 36 cases: demo5–demo8 and the 32 bench_nested cases. The pipeline is described in
[docs/pipeline/README.md](../../docs/pipeline/README.md). The scorer and every metric used here
are defined in [docs/pipeline/evaluation.md](../../docs/pipeline/evaluation.md): [§4, method-agnostic
scorer](../../docs/pipeline/evaluation.md#4-the-method-agnostic-scorer-aa_eval) for the metrics; [§5, translators](../../docs/pipeline/evaluation.md#5-translators-aa_export-aa_import-coral_flank) and [§5.4, CoRAL flank fix](../../docs/pipeline/evaluation.md#54-the-coral-flank-fix-coral_flank) for the
translators. The plan and the state of the work are in [NEXT_STEPS.md](NEXT_STEPS.md), and the
environments in [envs/README.md](envs/README.md).

> **Provenance:** (c) new, RCK yubi_dev@c2ecf50 (environments), @2f921a5 / @1e21c8e (export,
> import, round trip), @1c0f3b4 (head-to-head runs), @bb15c21 (CoRAL flank fix), @b4ff12d
> (scoring). **Code:** `results/graph_route/{graph_route.py,score_h2h.py,roundtrip.py,check_coral_format.py}`,
> `rck/ecdna/{aa_export.py,aa_import.py,aa_format.py,coral_flank.py,aa_eval.py}`.

## 1. Design

### 1.1 Graph arms

> **Code:** `graph_route.py::graphs_of`, `graph_route.py::prepare`.

| arm | graph | role |
|---|---|---|
| `truth_rp` | the truth mixture graph of amplicon-allele-analysis `simulations/truth_mixture.py`, with synthetic reads and long-read path constraints (`$G/<case>/truth/graph_reads_paths/mixture_amplicon<N>_graph.txt`), one run per truth amplicon: 39 graphs | oracle headline |
| `inferred` | our inferred union subgraph as an AA graph, written by `python -m rck.ecdna.aa_export subgraph` (`$G/<case>/inferred/graph/`): 36 graphs | realistic input |
| `truth_plain` | the truth mixture graph without reads or path constraints | secondary arm (priority 2 g), not in the default LIST |

How the truth graph's evidence is generated:

- synthetic reads: $\mathrm{NumberOfReadPairs} = \operatorname{round}(\mathrm{CN} \cdot R)$ on
  breakpoint edges, with $R = 6.5$;
- coverage: $\mathrm{AverageCoverage} = \mathrm{CN} \cdot D$, with $D$ from the case's
  `scale.json`;
- path constraints: the run's long-read paths whose junctions are mixture discordant edges.

The options are `truth_mixture.py --reads-per-copy`, `--depth-per-copy` and `--path-constraints`.

### 1.2 Methods and what each is given

> **Code:** `graph_route.py::run_coral`, `graph_route.py::run_ce`, `graph_route.py::rck_argv`.

| method | version | command and settings | input graph |
|---|---|---|---|
| CoRAL — (a) existing [Zhu et al. 2024](../../docs/pipeline/references.md#zhu2024), [CoRAL 3.1.1 code](../../docs/pipeline/references.md#coral_code) | 3.1.1 (commit 7f3bf317) | `coral cycle --solver gurobi_direct --solver-time-limit 1800 --threads T`; other options at CoRAL's defaults (`--cycle-decomp-mode max_weight`, alpha 0.01); `--global-time-limit` default 21600 | the graph with 10 kb flanks at every interval end (`rck.ecdna.coral_flank`; §1.3) |
| CycleExtractor — (a) existing [Faizrahnemoon et al. 2026](../../docs/pipeline/references.md#cycleextractor) | 1.0.0 (commit 1b29b33) | `CE.py --solver gurobi` at its defaults; its time limit is hard-coded at 7200 s per MILP; `CE_SOLVER=highs` runs it with HiGHS [Huangfu & Hall 2018](../../docs/pipeline/references.md#huangfu2018), without the licence | the graph as is |
| RCK-ecDNA — (c) new | RCK yubi_dev | `rck-ecdna --inputs subgraph` on `aa_import GRAPH --imputed none`; total-only (`--min-snps-allele 100000000 --soft-junction-hap --soft-junction-penalty 0`); the bench's λ procedure (`--step-penalty 0.1 --keep-pool --pool-lambdas 0 0.02 0.05 0.1 0.2 0.5`, answer re-ranked at λ = 0.2 when scoring); `--time-limit 1800` per MILP; `--threads T`; single seed (Gurobi's default) | the graph only, as CoRAL and CE see it: flank background 2, no `cn_profile`, no allele data, so `--block-flip` does not apply |
| native — (c) new | RCK yubi_dev | RCK-ecDNA's bench_nested runs at λ = 0.2 (its own subgraph, allele data, its own imputed edges), re-scored from the AA exports of `results/bench_nested/diag/aa_regression.py` | its own subgraph |
| original RCK — (a) existing [Aganezov & Raphael 2020](../../docs/pipeline/references.md#aganezov2020) | upstream 4c6c288 + `rck/version-fix` | pending (NEXT_STEPS priority 2 f): a converter to RCK's input format, then graph-level scoring | – |

A truth graph has no `paths.tsv`, so for the RCK-ecDNA arm on it the truth folder and sample come
from the matching bench case (`score_h2h.py::RUN_CASE`): demo6 → `demo6_union`, demo7 →
`demo7_totalonly`. The native arm's demo7 is `demo7_blockflip` (`NATIVE_CASE`), as in
`native_aa_metrics.tsv`.

### 1.3 Fairness decisions

The user took these decisions on 2026-11-03 (NEXT_STEPS.md, "Decisions"; and "Current plan
(2026-11-04)"):

1. **Candidate edges.** The head-to-head uses `aa_import --imputed none`: RCK-ecDNA sees exactly
   the graph CoRAL and CycleExtractor see.
   - RCK-ecDNA's imputed candidates are written only as `#rck_imputed` comment lines, which the
     other two skip.
   - RCK-ecDNA's native run is reported as a separate `native` arm.
   - `--imputed regenerate` is kept for Task 3's realistic part, never in the head-to-head.
2. **Truth graphs.** The oracle headline uses `graph_reads_paths`; `graph_plain` is secondary.
3. **Modes.** Every table is strict and best-of.
4. **Haplotype-blind comparison.** Cross-method comparisons use the haplotype-blind columns
   (`composite_h0`, blind one-to-one matching). CoRAL and CE do not phase, and RCK-ecDNA's
   graph-route runs are total-only.

**CoRAL flank fix** — (c) new, RCK yubi_dev@bb15c21, `rck/ecdna/coral_flank.py`.

- *The failure.* CoRAL 3.1.1 balances an interval's two end nodes only with its s/t (source and
  sink) edges (`coral/models/constraints.py::get_eulerian_node_constraints` in CoRAL 3.1.1 [CoRAL 3.1.1 code](../../docs/pipeline/references.md#coral_code)).
  A cycle through a breakpoint on an end node is therefore infeasible. Our graphs are cut at
  breakpoints: 69 of 75 graphs end at one. On them, CoRAL explained 0 copy number and returned 0
  cycles (Slurm 32104590; the CoRAL tasks were cancelled).
- *The fix.* `coral_flank` adds, at every run end, a 10 kb flank sequence edge (CN = the node's
  balance residual, no reads) and a concordant join. It writes the result in CoRAL's order, with
  the path constraints renumbered.
- *The effect.* Offline with CoRAL's own model, the truth cycles that are feasible go from 28 to
  66 of 67. Test job 32104711 recovered demo5's truth cycle at CN 16.

Only CoRAL is given the flanked graph. CycleExtractor connects s/t at every node, and RCK-ecDNA
has its own flank joins. The full definition is in
[docs/pipeline/evaluation.md §5.4, CoRAL flank fix](../../docs/pipeline/evaluation.md#54-the-coral-flank-fix-coral_flank).

### 1.4 Scoring (`score_h2h.py`)

> **Provenance:** (c) new, RCK yubi_dev@b4ff12d. **Code:** `score_h2h.py::one`, `::tasks`,
> `::normalise_cycles`, `::merged_flags`, `::summarise`, `::scenario_summary`, `::graph_rows`.

Every run of `graph_route.list` and the native arm go through the same call, `aa_eval.score`
(defined in [evaluation.md](../../docs/pipeline/evaluation.md)), at the scorer's defaults:

- strict rule `auto`: `cn-frac:0.1` for CoRAL and CE (no status field), `status` for RCK-ecDNA;
- species tolerance 10 kb;
- graph tolerance 100 bp.

What each method's answer is:

- **cycles:** CoRAL `coral_amplicon<N>_cycles.txt` (on the flanked graph); CE
  `ce_amplicon<N>_cycles.txt`; RCK-ecDNA the pool's rank 1 re-ranked at λ = 0.2, plus 2 same-fit
  alternatives (`aa_export.export`), as in the bench's λ = 0.2 results;
- **graph:** the input graph every method was given (for CoRAL, the unflanked one). It is used
  for `cn_explained` only; graph-level metrics are computed once per input graph
  (`h2h_graph_metrics.tsv`), not per method;
- **paths:** not passed, so `subwalks_satisfied` is not computed (§5).

**Linear walks (the `Path=` workaround).** CoRAL 3.1.1 and CycleExtractor write a linear walk as
`Path=k;...;Segments=0+,...,0-`. `aa_format.read_cycles` reads `Cycle=` lines only and skips
these. `rck-ecdna` writes its linear walks AA's way, as `Cycle=` lines with segment 0, which are
scored. So that every method's linear walks are treated alike, each answer is first copied to
`<score dir>/input/` (`normalise_cycles`):

- **headline:** `Path=` is rewritten as `Cycle=`, so linear walks are scored, as the scorer
  documents for AA-style paths;
- **`--drop-paths` (sensitivity):** every linear walk is removed: `Path=` lines and `Cycle=` lines
  whose first segment is 0, from every method. Outputs carry the suffix `_cycles_only`.

**Per-truth-species counts** (strict, blind matching), benchmark-specific:

- **missed (strict):** no strict cycle is matched one to one (blind).
- **missed (best-of):** $\mathrm{composite\_h0} \le 0$, aa_eval's best-of rule.
- **merged:** the bench_nested `merged` rule
  ([../bench_nested/README.md](../bench_nested/README.md), §3.2), applied per truth species on
  haplotype-stripped cycles. With $R(t,s)$ the unphased junction recall of truth species $t$ in
  strict cycle $s$ at 10 kb, and $\mathrm{carries}(s) = \{t : R(t,s) \ge 0.5\}$:

  ```math
  \mathrm{merged}(t) = \big[\, \exists\, s \in S^{\mathrm{strict}} : t \in \mathrm{carries}(s) \ \land\ |\mathrm{carries}(s)| \ge 2 \,\big]
  ```

  where $S^{\mathrm{strict}}$ is the method's strict cycles.
- **false:** `n_false` counts unmatched strict cycles. `n_false_linear` counts those of them that
  are paths: CoRAL and CE background walks through the flanks, and CE's linear walks.
- **exact:** AA/CoRAL exact reconstruction at 100 bp, haplotype-blind
  ([evaluation.md §4.4, exact cycle](../../docs/pipeline/evaluation.md#44-species-level-metrics-aa_species_metricstsv)).

Merged, missed-or-merged and false counts are reported for the strict mode only.

## 2. Head-to-head results (2026-11-06, single seed)

**RCK-ecDNA's numbers are single-seed.** In the round trip, its decomposition changed in 29 of 36
cases with an unchanged flow objective. The flow MILP has many equally good integer flows, the
solver's choice among them depends on model size and order, and the decomposition follows that
choice. Seed repeats are pending (NEXT_STEPS step e). Its per-case scores swing more than the
other methods' between truth and inferred graphs. It scores 0.5–0.6 *lower* on the cleaner truth
graph in nest_cross_1to1_s509, nest_cross_2to1_s103 and nest_samehap_3to1_s304, which fits
tie-breaking. Read its case-level numbers with that in mind.

All 225 runs finished (`graph_route.py status`: 225 done, none failed): 75 graphs × 3 methods,
with CoRAL rerun in job 32119533. Every CoRAL run has ≥ 1 cycle. The runs are scored on the 67
truth species of the 36 cases, all circular. Values are the mean composite-H0 over cases,
strict / best-of. RCK-ecDNA* is total-only on the shared graph. Native is its bench run, which has
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

**Main results.**

- **Shared total-only graph.** No method is clearly ahead.
  - On truth graphs, CE leads on strict composite-H0 (0.80, 43 exact); RCK-ecDNA and CE tie on
    best-of (0.82).
  - On the inferred graphs RCK-ecDNA leads (0.69 / 0.76, against CoRAL 0.67 / 0.71 and CE
    0.60 / 0.67).
  - These margins are mostly within RCK-ecDNA's tie-breaking noise.
- **Allele data clearly separates.** The native arm (0.82 / 0.91, 50 exact) beats every
  total-only arm. The gap is largest in the nest_cross scenarios (0.86 against at most 0.63),
  where nested species on different haplotypes cannot be separated from total copy number alone.

By scenario:

- CoRAL and CE solve the fold-back controls (ctrl_dup, 1.00 on truth). There RCK-ecDNA splits the
  truth walk into two cycles at an identical copy number (s708, s710; the same in the native arm).
- CE wins nest_samehap_3to1 and 6to1 on truth graphs (0.90, 0.92), where its first greedy cycle
  satisfies the path constraints. The same greedy step makes one long walk that merges both demo8
  species at the minor CN. On inferred graphs it emits many linear walks (72 false). CE loses the
  most from truth to inferred graph (0.80 → 0.60).
- CoRAL is weakest on nest_cross (0.35 strict on 1to1) and, on truth graphs, on samehap_6to1.

Every arm loses 0.05–0.20 from truth to inferred graph. The inferred graphs keep every breakpoint
(`h2h_graph_metrics.tsv`: `bp_accuracy` 1.0, no false edge, 38 of 39 `error_free_merged`). The
loss therefore comes from segmentation and copy number, not from missing junctions.

**What depends on the CoRAL flank fix and on CoRAL's multiplicity bound.**

- Every CoRAL number depends on `coral_flank`: without it, CoRAL returned 0 cycles.
- No CoRAL *cycle* uses a flank segment. Flanks appear only in CoRAL's linear s/t paths: 22 strict
  paths on truth graphs and 15 on inferred. Of these, 12 and 7 are the one-to-one match of a truth
  species, which gives partial credit and is never exact.
- `aa_format.read_cycles` skips CoRAL/CE `Path=` lines, so the headline reads them as AA-style
  linear walks, as RCK-ecDNA's are scored. Dropping every linear walk (`h2h_*_cycles_only.tsv`)
  gives CoRAL truth strict 0.71 (missed 5 → 15) and CE 0.79; RCK-ecDNA is unchanged. No ranking
  changes.
- nest_samehap_6to1_s804: the CN-40.8 species exceeds CoRAL's read-clustered discordant
  multiplicity bound, and CoRAL's heaviest cycle is CN 6.8. CoRAL scores 0.39 truth / 0.35
  inferred (CE 0.54 / 0.52, RCK-ecDNA 0.72 / 0.71). Leaving s804 out raises CoRAL's overall truth
  strict mean by 0.01 (0.742 → 0.752).

**Not computed.** `subwalks_satisfied` is not computed (§5). CoRAL's own
`Path_constraints_satisfied` is in its cycles files.

Outputs go to `/n/fs/yubi-iw/runs/graph_route/` (project space, not purged; rebuild with the
commands in Appendix A).

## 3. What the AA graph holds

> **Code:** `rck/ecdna/aa_export.py` (`subgraph` subcommand), `rck/ecdna/aa_import.py`,
> `rck/ecdna/aa_format.py::write_graph`.

`aa_export subgraph` writes the graph RCK-ecDNA solves on (`read_subgraph`: phantom segments
contracted, parallel imputed edges merged), with the copy numbers it uses:

- segment CN = `CN_TOTAL`;
- junction CN = reads ÷ junction reads per copy θ. This is the observed adjacency copy number.
  RCK-ecDNA observes junctions this way, not through the rounded `CN` column of `edges.tsv`.

Positions are subgraph coordinates + 1, as for exported answers.

| line | content | seen by CoRAL / CE |
|---|---|---|
| `sequence` | one per segment, in CoRAL's sort order; AverageCoverage = depth | yes |
| `concordant` | reference adjacencies between segments; an unobserved join (imputed grid point) gets the smaller segment CN | yes |
| `discordant` | observed novel adjacencies; fold-backs (both ends on one extremity) as self-loops | yes |
| `path_constraint` | read paths of `paths.tsv` over segments, reference adjacencies and observed novel adjacencies, as CoRAL 3.1.1 tokens; support = `#NOV` | yes |
| `interval` | runs of contiguous segments (what CoRAL infers itself) | ignored |
| `#rck_imputed` | RCK-ecDNA's imputed candidates (no read evidence) | no |
| `#rck_outside` | joins to a region's outer flank (an extremity without a segment); CoRAL and CE add their own s/t ends instead | no |
| `#rck_node`, `#rck_seq`, `#rck_edge`, `#rck_pc`, `#rck_path` | the subgraph's node, segment, adjacency and path ids and orders, so that the import is exact | no |

`aa_import` inverts it:

- `--imputed comment` restores the `#rck_imputed` lines;
- `--imputed none` leaves them out (what CoRAL and CE see);
- `--imputed regenerate` proposes RCK-ecDNA's own closing candidates on the graph
  (`subgraph_construct.py`'s rule adapted to sequence edges; NEXT_STEPS.md);
- `#rck_outside` lines are always restored.

`scale.json` is fitted from the graph:

- junction scale = median of reads ÷ CN;
- a graph without reads gets WEIGHT = CN and scale 1;
- no depth fit from a CoRAL-rewritten graph.

**Why the ids are kept in comments.** RCK's junction CN is reads ÷ θ, and the order of the ILP's
variables depends on the ids. Keeping them in `#rck_*` lines makes the round trip exact.

## 4. Validation done offline (2026-11-03, no Gurobi)

> **Code:** `roundtrip.py`, `check_coral_format.py`; tests `tests/ecdna/test_aa_import.py`,
> `tests/ecdna/test_coral_flank.py`.

- **Round trip, 36/36 cases** (demo5–demo8 + 32 bench; `roundtrip_prepare.tsv`):
  - `--imputed comment` with the original `cn_profile` gives back the normalised graph exactly:
    node, segment and adjacency ids and orders, reads, observations, the junction scale to the
    bit, and the read paths;
  - `--imputed none` is exact apart from the imputed candidates;
  - graph-only `none` differs only in the flank background (2 instead of the `cn_profile`).

  Unit tests: `tests/ecdna/test_aa_import.py`. The Slurm part (RCK-ecDNA on each arm,
  `roundtrip.sbatch`, job 32096382; `roundtrip_compare.tsv` / `.txt`) found:
  - with imputed edges restored, every case is identical;
  - without them, the flow objective is unchanged except in the two cases whose flow used one
    (ctrl_dup_s710, nest_samehap_6to1_s804);
  - the decomposition changes in 29 of 36 cases at an unchanged flow objective. This is a
    property of RCK-ecDNA, not a translator bug (NEXT_STEPS, "Current plan").
- **CoRAL 3.1.1, 36 exported + 78 truth graphs** (`check_coral_format.py`):
  - parse counts equal ours;
  - `infer_amplicon_intervals` gives our intervals and keeps the edge order;
  - CoRAL's writer round trip keeps edges, CN and constraints;
  - all 418 path constraints of the exported graphs were rebuilt independently on CoRAL's own
    adjacencies and written by its `path_to_str`, with our exact tokens;
  - every truth constraint crosses exactly the junctions of its read.
- **CycleExtractor (HiGHS)** on demo5, demo7 and ctrl_dup_s702 (exported) and ctrl_dup_s710
  (truth, with reads and constraints): it parses them and writes cycles that
  `aa_format.read_cycles` reads. On the ctrl_dup_s710 truth graph it returns the truth cycle
  exactly, both fold-backs included.

## 5. Known limitations (recorded, not worked around)

- **Single seed** for every method; see §2. Tie-sensitivity runs are NEXT_STEPS step e.
- **`subwalks_satisfied` is not computed.** `aa_eval.read_paths` needs `chr:pos` ids, and the
  inferred graphs' `paths.tsv` carry the subgraph's `imputed_ext_*` ids.
- **`aa_format.read_cycles` skips `Path=` lines.** `score_h2h.py` works around it (§1.4); the
  scorer itself still skips them.
- **CoRAL ignores breakpoints on interval end nodes** (worked around for CoRAL only). The two outer
  nodes of each run of contiguous sequence edges are CoRAL's end nodes, balanced only by its s/t
  edges (`coral/models/constraints.py::get_eulerian_node_constraints`), so a cycle through a breakpoint there is infeasible.
  Our graphs end at breakpoints (69 of 75), so `graph_route.py` gives CoRAL the graph with 10 kb
  flanks (`python -m rck.ecdna.coral_flank GRAPH OUT [--flank 10000]`); NEXT_STEPS.md has the
  checks.
- **CoRAL's log reports "Total length weighted CN" as 0** in greedy mode, whatever the answer.
  Judge a CoRAL run by its number of cycles.
- **CoRAL re-sorts without renumbering.** `coral cycle` calls `coral/breakpoint/breakpoint_graph.py::infer_amplicon_intervals`, which
  sorts sequence and concordant edges (natural chromosome order) after the path constraints were
  parsed by file index. The exporter and `truth_mixture.py --path-constraints` write CoRAL's
  order. A graph from elsewhere with constraints in another order is silently corrupted by CoRAL
  (ctrl_separate_s603, chr15 + chr7, would have been).
- **CoRAL's writer** puts `lr_nc / size` in AverageCoverage and truncates read counts to integers.
  Never fit a depth scale from a CoRAL-rewritten graph (`aa_import` refuses to).
- **Direction of a path constraint.** CoRAL writes a path in canonical direction. RCK-ecDNA
  matches paths in both directions, and `#rck_pc` restores the read's own direction anyway.
- **What CoRAL and CE cannot see:**
  - allele copy numbers, phase blocks and junction haplotype labels. The import is total-only;
    `_segment_cn.tsv` keeps alleles for scoring;
  - RCK-ecDNA's flank joins. They add their own s/t ends; the import restores the joins even with
    `--imputed none`, as they are RCK-ecDNA's counterpart of those ends.
- **Unobserved concordant joins** (at imputed grid points, 10 in demo7) get the smaller
  neighbouring segment CN in the AA graph; RCK-ecDNA has no observation there.
- **Truth and inferred read paths differ.**
  - The truth graph gets every read of `obs/path.tsv` whose junctions are mixture edges (demo7:
    12). The subgraph's `paths.tsv` was filtered by `rck_inputs.clip_path` (demo7: 3).
  - Synthetic truth reads are $\operatorname{round}(\mathrm{CN} \cdot R)$, so the fitted junction
    scale is close to R but not exactly R (demo7: 6.497 for 6.5).
- **CycleExtractor's fold-back gadget** (§6): its two half-edges each carry the junction's read
  count into CE's multiplicity bound.

## 6. Fold-back junctions

A fold-back is an inverted duplication whose two breakends are the same extremity. It appears as
a discordant self-loop. 14 of the 32 bench graphs have one; the ids are in each `_export.json`
under "foldbacks":

- ctrl_dup s702, s703, s708, s709, s710;
- nest_cross_1to1 s505, s506, s508;
- nest_cross_3sp s402, s403;
- nest_samehap_2to1 s201;
- nest_samehap_6to1 s803, s804, s805.

How each method handles them:

- **CycleExtractor 1.0.0 does not drop them.** This corrects an earlier note in NEXT_STEPS.md.
  `Create_graph` replaces each self-loop by a gadget: `node -> i`, a 1 bp sequence edge
  `i -> i+1`, `i+1 -> node`. The `if start != end` branch only routes loops to it.
  - On ctrl_dup_s702, its cycle turns at the fold-back (`...5+,6+` then `6-,5-...`) and satisfies
    the fold-back path constraint `e6+:2,d4+:1,e6-:2`.
  - On the ctrl_dup_s710 truth graph it recovers the truth cycle with both fold-backs.
- **CoRAL models self-loops directly.** `coral/breakpoint/breakpoint_graph.py::add_discordant_edge` lists a self-loop twice at its node,
  so the cycle model's Eulerian balance counts one traversal twice there. Its CN fit
  (`coral/datatypes.py::DiscordantEdge.is_self_loop`) treats self-loops separately.
- **RCK-ecDNA**: a self-loop's haplotype labelling is same-haplotype only (AA or BB), as in RCK
  [Aganezov & Raphael 2020](../../docs/pipeline/references.md#aganezov2020).

## Appendix A. Commands

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

# flank a graph for CoRAL
$PY -m rck.ecdna.coral_flank GRAPH OUT [--flank 10000]

# check graphs with CoRAL itself (no solver), and CycleExtractor with HiGHS (no licence)
PYTHONPATH=~/RCK $CPY results/graph_route/check_coral_format.py GRAPH... [--paths-from-export]
PYTHONPATH=~/RCK $CPY results/graph_route/check_coral_format.py --ce OUTDIR GRAPH...
PYTHONPATH=~/RCK $CPY results/graph_route/check_coral_format.py --golden OUT.tsv   # reference tokens for test_aa_import.py

# the head-to-head
cd ~/RCK/results/graph_route
python graph_route.py prepare            # offline, no Gurobi: graphs, imports, LIST (done 2026-11-04)
python graph_route.py coral-inputs       # offline: the flanked graphs CoRAL runs on
python graph_route.py list               # one run per line
L=/n/fs/scratch/$USER/graph_route/h2h/logs; mkdir -p $L
sbatch --output=$L/%x_%A_%a.out graph_route.sbatch   # --array=0-224%2: lines 0-149 CoRAL + RCK-ecDNA, 150-224 CE
python graph_route.py status             # offline, afterwards
python score_h2h.py [--jobs 16]          # offline scoring: h2h_*.tsv
python score_h2h.py --drop-paths         # sensitivity: h2h_*_cycles_only.tsv

# the round trip
python roundtrip.py prepare [--arms original rt_comment rt_none_profile rt_none]   # offline
python roundtrip.py list
sbatch roundtrip.sbatch
python roundtrip.py compare              # offline, after the array
```

Run infrastructure — (b) patch, no effect on results:

- `%2` arrays: CoRAL, CycleExtractor (`--solver gurobi`) and RCK-ecDNA all use the academic Gurobi
  WLS licence, which allows 2 sessions;
- Gurobi "too many sessions" errors are retried up to 5 times, 2 min apart
  (`graph_route.py::RETRIES`);
- a run is skipped when its outputs exist;
- environments, CycleExtractor and outputs are in `/n/fs/yubi-iw`;
- Slurm logs stay on scratch;
- `GRAPH_ROUTE_OUT` moves the outputs.

The sbatch scripts request:

| script | memory | time | notes |
|---|---|---|---|
| `graph_route.sbatch` | 4 GB | 3 h | RCK-ecDNA peaked at 366 MB in the round trip |
| `roundtrip.sbatch` | 2 GB | 1 h 15 min | `--array=0-143%2` (36 cases × 4 arms) |

## Appendix B. Files

| file | what |
|---|---|
| `envs/` | lock files of `/n/fs/yubi-iw/conda/envs/coral-3.1.1` (CoRAL and CycleExtractor both run with it) |
| `fixtures/` | graphs written by CoRAL 3.1.1's own writer (`coral_rewrite.py`): the reference layout |
| `coral_rewrite.py` | parse a graph with CoRAL and write it back with CoRAL's writer |
| `check_coral_format.py` | check AA graphs against CoRAL 3.1.1's parser, sorter and writer; `--golden`; `--ce` (CycleExtractor, HiGHS) |
| `roundtrip.py`, `roundtrip.sbatch`, `roundtrip.list` | the round-trip test (subgraph → AA → subgraph → RCK-ecDNA total-only) |
| `roundtrip_prepare.tsv`, `roundtrip_compare.tsv`, `roundtrip_compare.txt` | its offline part (per case, the graph comparison of each arm) and the comparison of the runs |
| `graph_route.py`, `graph_route.sbatch`, `graph_route.list` | the head-to-head runs: CoRAL, RCK-ecDNA (`aa_import --imputed none`) and CycleExtractor on the truth `graph_reads_paths` and inferred graphs, one run per line, `run.json` per run (method, version, solver, time limit, input graph, command, elapsed, attempts) |
| `graph_route_prepare.tsv` | per graph: input file and the import's counts and fitted scale |
| `$G/<case>/runs/<graph>/coral_input/` | the graph CoRAL runs on: the input graph with 10 kb flanks at every interval end |
| `$G/<case>/runs/<graph>/rck_input/amplicon<N>/` | `aa_import --imputed none` of the graph |
| `$G/<case>/runs/<graph>/<method>/amplicon<N>/` | one run, with `run.json` |
| `native_aa_metrics.tsv` | the native arm: RCK-ecDNA's bench runs at λ = 0.2, scored by aa_eval (rows of `results/bench_nested/aa_metrics.tsv`) |
| `score_h2h.py` | scores every run and the native arm with aa_eval (offline) |
| `h2h_table.tsv` | case × graph × method × mode |
| `h2h_scenario.tsv` | the same per scenario |
| `h2h_species.tsv` | per truth species |
| `h2h_graph_metrics.tsv` | graph-level metrics, once per input graph |
| `h2h_*_cycles_only.tsv` | the `--drop-paths` sensitivity tables |
| `$G/<case>/score[_cycles_only]/<graph>/<method>/amplicon<N>/` | per-run scorer output (`aa_eval.write`) |
