# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This file covers the `yubi_dev` branch relative to `rck/version-fix`.

## [Unreleased]

### Added

- **`rck-ecdna`** (`python -m rck.ecdna_run`, package `rck.ecdna`): allele-specific
  reconstruction of ecDNA from one amplicon breakpoint subgraph, built on RCK's diploid
  interval adjacency graph and haplotype-aware balance constraints.
  - `subgraph_io`: reads `subgraph_construct.py` folders (`nodes.tsv`, `edges.tsv`,
    `phase_blocks.tsv`) into a normalized graph:
    - contracts phantom 1 bp segments;
    - merges parallel imputed edges;
    - marks terminals;
    - attaches allowed phasings from junction haplotype labels;
    - splits regions at gaps over 1 Mb.
  - `scale`: converts depth and junction reads to copy numbers. Depth per copy comes from an
    override, HATCHet `best.seg.ucn`, the mode of balanced background blocks, or a
    tumor/normal ratio; junction reads per copy from background extremities. It also
    gives allele-specific segment copy numbers and flank background copy numbers.
  - `ilp_circular`: a two-layer flow MILP adapted from `OptModelMultiClone`:
    - a strictly balanced ecDNA layer;
    - a background layer that ends only at terminals;
    - imputed candidate edges with a usage penalty;
    - optional infinite-sites constraint.
  - `decompose`: minimum-species decomposition of the ecDNA flow:
    - integer weights and edge multiplicities, with exact linearization;
    - lazy connectivity cuts;
    - alternatives from the solution pool or no-good cuts;
    - `low_weight` flags;
    - ranking by residual and number of significant species.
  - `walk`: alternating Hierholzer walks, enumeration of distinct orderings up to rotation
    and reversal, and ranking of orderings by long-read support.
  - `output`: `decompositions.tsv`, `flow_segments.tsv`, `flow_adjacencies.tsv`,
    AmpliconArchitect-style graph and cycles files, `observations.json` and `summary.json`.
  - `evaluate` and `rck-ecdna-eval`: scoring against simulator truth:
    - truth read from `truth/aa`, with per-sample copy numbers from `proportions.tsv`;
    - truth species scoped to the subgraph's regions;
    - haplotype Jaccard, junction recall and precision, order match and copy-number error;
    - errors are reported instead of failing the run.
  - `gurobi_env`: thread defaults from Slurm, a quiet `gurobipy` logger, and model disposal.
- **Linear structures in the ecDNA layer**:
  - penalized linear ends at copy-number boundaries and region flanks
    (`--linear-ends`, `--linear-end-penalty`);
  - linear species and a `topology` column;
  - an amplicon call (`ecDNA` / `linear` / `mixed` / `none`) with `circular_fraction`;
  - linear walks written with AA's `0+`/`0-` source;
  - linear walks are never order-matched to a circular truth.
- **Long reads in the decomposition**: multi-junction paths become read constraints, with a
  capped reward (`--path-weight`) or a required explained share (`--path-fraction`, as in
  CoRAL). After solving, each read is checked in order against the orderings and
  `reads_explained` is reported.
- **Fractional copy numbers** (`refit`): after the integer structure search, species weights
  are refitted against the observations with a linear program (scipy, HiGHS). Species below
  `--refit-min-copies` are dropped. Outputs gain `weight_int`, `data_objective` and `pruned`,
  and the evaluation gains `cn_rel_err_int`.
- **Subgraph-only input mode** (`--inputs auto|subgraph|windows`): reads `scale.json`, the
  `CN_*` edge columns, `cn_profile.tsv` and `paths.tsv` written by amplicon-allele-analysis
  `subgraph_construct.py`, so no windows, extremities or HATCHet input is needed.
  `--baseline-cn` defaults to the subgraph's cell ploidy.
- **Recovery matrix** (`recovery_matrix.tsv`, and `recovery_matrix` / `mean_composite` in
  `evaluation.json`): one row per truth species.
  - Scores against the rank-1 decomposition:
    - segment recall, precision and Jaccard;
    - haplotype Jaccard;
    - junction recall, precision and F1, unphased and phased;
    - phasing accuracy, as given and A/B swap-invariant;
    - order and topology match;
    - copy-number error;
    - an equal-weight `composite` in [0, 1].
  - `best_rank`, `best_inferred` and `best_composite`: the best alternative decomposition.
  - Out-of-scope species are listed but not scored. A missed species scores 0.
  - `evaluation.tsv` gains the new per-pair columns.
- **Figures** (`rck.ecdna.plot`, `rck-ecdna-plot`): at the end of a run, `plots/<name>_ecdna.pdf`
  and one PNG per page.
  - One page per species of the top decomposition. On a cnplot `GenomeAxis` (amplicon-allele-analysis
    `third_party/cnplot`) restricted to the species' segments, it shows:
    - the haplotype-coloured topology (CoRAL-style);
    - the matched truth species;
    - the junction arcs;
    - A/B copy number: observed points, the inferred line and this species' share.
  - Beside the tracks, a ring diagram of the walk.
  - A recovery-matrix heatmap page after an evaluation.
  - Options: `--no-plots`, `--plot-ranks`, `--plot-low-weight` and `--cnplot-path` / `$CNPLOT_PATH`.
  - matplotlib and cnplot are optional. Without them, or on a plotting error, the run logs a
    warning and still succeeds.
- `bg_balance_weight`: a small tie-break toward a heterozygous background.
- **Phase-block flips** (`--block-flip unoriented`, default `off`): allele observations per
  phase-block piece, with the A/B labelling of every block the subgraph left unoriented chosen
  by the flow MILP.
  - `subgraph_io`: reads `segment_blocks.tsv` (written by amplicon-allele-analysis
    `subgraph_construct.py`) and `phase_blocks.tsv`, and maps each piece to its segment by the
    parent edge's key, through phantom contraction. Without `segment_blocks.tsv` the run logs a
    warning and keeps the segment observations.
  - `scale`: one observation per piece. Pieces below `--min-snps-allele` SNPs are total-only.
    Pieces in blocks with orientation 1 or −1 are fixed; pieces in blocks with orientation 0 or
    missing from `phase_blocks.tsv` follow the block's flip.
  - `ilp_circular`:
    - one binary flip per unoriented block, shared by its pieces across segments;
    - piece targets `a + f (b − a)` and `b + f (a − b)`, linear in the flip;
    - the segment's data weight split over its pieces by length, so the objective scale is
      unchanged;
    - a 10⁻⁶ cost per flip, so a block without evidence stays as observed;
    - after solving, `delta_if_reversed` per block, and `unresolved` blocks (reversal within
      `--flip-tol`, default 0.01) kept as observed;
    - optional `--flip-oriented-penalty` lets oriented blocks flip at a cost.
  - The flips are fixed after the flow: the refit, `flow_segments.tsv` (segment observations
    pooled over the pieces after the flips) and the CN-A / CN-B figure points use them.
  - New `block_flips.tsv`; `observations.json` gains `pieces`, and `summary.json` gains
    `block_flip` (mode, counts per status, flipped blocks). With `--block-flip off` no output
    changes.
  - demo7, on the union subgraph rebuilt with `segment_blocks.tsv` (run on ionic; outputs in
    amplicon-allele-analysis `simulations/runs/demo7_nestedcycles_ecdna/rck_out_union_blockflip/`):
    - the flow flips chr4:185 and chr4:186 and keeps chr4:171, chr4:172 and chr4:187 as
      oriented; the 12 blocks of the unamplified 172.4–185.5 Mb gap are `unresolved`;
    - the flow objective is 9.13, with no imputed edge (39.26 with `--block-flip off`, which
      reproduces `rck_out_union/`, and 8.8 with total-only observations);
    - both nested species are recovered as separate circles in the rank-1 decomposition: 16.43
      copies against a true 16.90 and 9.08 against 9.10. Composites are 0.994 and 0.9996 (0.43
      and 0 before), with segment, haplotype, junction and phasing scores of 1.0 and matching
      order. The refit drops one 0.39-copy leftover.
  - Tests:
    - `test_block_flip.py`: inputs without Gurobi; small MILPs and runs;
    - `test_demo7.py`: needs a full licence and `RCK_ECDNA_DEMO7_SUBGRAPH`, and is skipped
      otherwise. It passes on ionic;
    - `tests/ecdna/helpers.py` gains `write_segment_blocks` and `have_full_gurobi`.
  - Documented in `docs/ecDNA.md`: the subgraph folder files, inputs, the model, outputs,
    figures and the tests table.
- **Walk complexity and ambiguity** (`decompose`), for nested species that the copy numbers
  cannot separate:
  - `--step-penalty λ` (default 0, so existing results are unchanged): the ranking's structure
    cost adds λ per step (segment traversal) of each significant species. This prefers separate
    simple walks over walks that re-traverse shared sequence, without pushing towards splits.
    - The penalty is applied in the ranking only. In the MILP it also charged the light
      leftover species and made the solver split true duplications (demo3, demo7). A weight
      threshold inside the MILP was gamed by splitting a merged walk into light copies.
    - With λ > 0 the solution pool keeps every solution within `--ambiguity-window` of the best
      objective instead of within 10%, then drops the ones the prior cannot prefer.
  - Ambiguity report (`--ambiguity-window`, default 2.0; 0 = off): one more MILP fixes the
    flow each edge covers to the top decomposition's and forbids two species on the same
    support. It then enumerates alternatives with exactly the same fit. `summary.json` gains
    `ambiguity` (`ambiguous`, `n_alternatives`, `gap`, `decided_by`, the alternatives) and a
    `steps` count per decomposition; an ambiguous top decomposition logs a warning.
  - Sweep of λ = 0, 0.02, 0.05, 0.1 and 0.2 on ionic, kept in `results/step_sweep/` (summary
    `sweep.tsv`, every run's output, the sweep script and `README.md`):
    - demo3, demo4, demo5, demo6 and demo7 with `--block-flip`: unchanged at every λ;
    - demo7 total-only: the merged walk (composites 0.51 / 0.34) becomes the two true species
      at λ ≥ 0.05 (0.98 / 0.97, no leftover);
    - demo8: the same correct answer, now decided by the step penalty (gap 1.9 at λ = 0.1)
      rather than a 0.004 tie-break.
    - λ = 0.1 is recommended for nested amplicons. The wider pool makes the demo8
      decomposition take 31 s instead of 8 s.
    - `results/step_sweep/sweep_v1_milp_penalty.tsv` keeps the first design's sweep, with the
      penalty in the MILP. It split true duplications on demo3 (composite 0.96 → 0.72) and on
      demo7 with block flips (0.77 for `ecdna2`).
  - Tests: `test_step_penalty.py`. Documented in `docs/ecDNA.md` (section 4).
- **Offline re-ranking over `--step-penalty`** (`--keep-pool`, `rck.ecdna.rerank`,
  `rck-ecdna-rerank`), so one run gives the answer of every λ > 0. λ acts only in the ranking,
  and every λ > 0 widens the pool to the same `--ambiguity-window`. The decomposition MILP and
  its pool therefore do not depend on λ, but the pruning of the widened pool does.
  - `--keep-pool` (needs `--step-penalty` > 0) also writes `pool/`. The run's own outputs are
    unchanged.
    - `pool.tsv`: every pool solution before the pruning, with what the ranking reads.
    - `decompositions.tsv` and `evaluation.tsv`: the same solutions, annotated, refitted and
      evaluated like ranked ones.
    - `pool.json`: the settings, the rank-1 solution per `--pool-lambdas` value (default 0, 0.02,
      0.05, 0.1, 0.2, 0.5), and those solutions' same-fit alternatives.
  - `python -m rck.ecdna.rerank RUN --lambdas ...`: per λ, the rank-1 solution, its refitted
    significant species, steps, residual, composite per truth species, and `ambiguous`,
    `decided_by`, `gap` and `n_alternatives`. It uses the run's own functions.
  - λ = 0 is only emulated, because a λ = 0 run solves a narrower pool whose ties can come out
    differently. On demo3, λ = 0 gives 12.77 copies and every λ > 0 an equivalent 12.70.
  - `decompose`: the pruning (`pool_keeps`, `prune_pool`), the ranking key (`ranking_key`,
    `rank_pool`) and the ambiguity report (`same_fit_structures`, `ambiguity_summary`) are now
    functions of recorded numbers, shared by the run and the re-ranking. `CycleDecomposer.solve`
    returns the unpruned pool, and `decompose()` prunes it as before. Results are unchanged: the
    existing `test_step_penalty.py` cases pass as before.
  - Tests in `test_step_penalty.py`:
    - pool ranking without Gurobi;
    - at 2:1 and 3:1, a pool kept at λ = 0.1 re-ranks to the rank-1 decomposition and ambiguity
      fields of direct runs at λ = 0.02–0.5;
    - the command line refuses `--keep-pool` at λ = 0, and `rerank` reproduces a run's rank 1.
  - The re-ranking was checked locally (restricted licence) against `results/step_sweep/sweep.tsv`
    on demo4 `amplicon_1`, demo5 and demo6. The informative cases, demo3, demo7 total-only and
    demo8, need a full licence and are the first step of `results/bench_nested/`.
- **`results/bench_nested/`**: a minimum benchmark for choosing `--step-penalty` and measuring
  nested ecDNA end to end.
  - Cases: 14 simulations from amplicon-allele-analysis `make_scenarios.py --set nested` (7
    scenarios × 2 checked seeds, with `--block-flip unoriented`) and the 8 step-sweep cases.
  - Each case gets two runs, λ = 0 and a `--keep-pool` run at 0.1. λ = 0.02, 0.05, 0.1, 0.2 and
    0.5 are re-ranked from the pool.
  - `run_case.py` and `bench_nested.sbatch`: a Slurm array with at most 2 concurrent tasks, the
    conda env on `PATH`, and each run skipped when its outputs exist.
  - `collect.py` writes:
    - `rerank_validation.tsv`;
    - `bench_sweep.tsv`: per case and λ, the rank-1 species against the truth, composites, a
      `merged` flag, ambiguity and runtime;
    - `pipeline.tsv`: per bench case at the chosen λ, subgraph scores, the flow, block flips
      against the truth, the RCK scores, and the stage behind each failure;
    - the results section of `README.md`: the smallest λ on the best plateau with no control or
      demo regression, and the λ range that gives the same answers.
  - The default `--step-penalty` stays 0.
  - Run on ionic, 20 of 22 cases. `nest_cross_3sp_s406` and `nest_cross_1to1_s502` were stopped
    at the session's time limit. Their decomposition MILPs take tens of minutes; the three-species
    case was at an 81% gap after 800 s.
    - Re-ranking check: the re-ranked λ > 0 composites and ambiguity fields match
      `results/step_sweep/sweep.tsv` at 32/32 (case, λ), including demo3, demo7 total-only and
      demo8, and the λ = 0 runs at 8/8. λ = 0 emulated from the pool misses demo3 (0.959 against
      0.960), the expected tie; that is why λ = 0 is run.
    - Mean composite over the 20 cases: 0.800 at λ = 0, 0.803 at 0.02, 0.896 at 0.05–0.2 and
      0.898 at 0.5. No control or demo case loses against λ = 0.
    - The rule picks 0.5 (the edge of the grid) because of one case that stays wrong at every λ.
      λ = 0.05 is the smallest value with the large gain. The README records both.
    - Six of 12 bench cases still fail at the chosen λ:
      - `ctrl_separate_s601` and `s602`: the subgraph drops segments of one species;
      - `nest_samehap_3to1_s301` and `nest_cross_3sp_s402`: a species is still merged;
      - `nest_cross_1to1_s501`: the nested species scores 0.56;
      - `ctrl_dup_s702`: one ecDNA phase block wrong after the flow, and an imputed edge.
    - 754 s of RCK wall time in all.
  - `run_case.py`: frees the Gurobi default environment at the end of each task
    (`disposeDefaultEnv`) and retries up to 5 times, 2 min apart. Sessions of finished WLS
    processes stayed counted ("5 active sessions for a baseline of 2"), and 3 tasks failed on it.
  - The runs live on scratch: `results/bench_nested/runs` is a git-ignored symlink.
  - `run_case.py` and `collect.py` read `BENCH_NESTED_ROOT`, as amplicon-allele-analysis
    `bench_nested.py` does, to take the bench cases from a benchmark built elsewhere
    (`configs/seeds.tsv`, `subgraphs/<run>/union` and `runs/<run>` under it). `collect.py` also reads
    `BENCH_NESTED_RESULTS`, which writes its tables and a results-only `README.md` to another
    folder, so a new benchmark does not overwrite this one's. Unset, nothing changes.
  - **New benchmark, bench_nested_v2** (on scratch, not committed). It was rebuilt with
    amplicon-allele-analysis `hap_labels: inherit`. Run on ionic: Slurm 32000320, the 14 bench
    cases. The 8 regression cases are links to the runs above. Tables:
    `/n/fs/scratch/$USER/bench_nested_v2/rck_results/`, with `comparison.md` and `comparison.tsv`
    against the committed tables above, which are unchanged.
    - 9 of the 14 bench cases keep their seed. Their walks and copy numbers are the same, and for
      the 5 cross and control seeds the haplotype labels change. The other 5 are replacement
      seeds of the same scenarios.
    - `ctrl_dup_s702` is fixed: 0.744 → 0.997 at every λ. It has no wrong junction label, no
      imputed edge and no block wrong after the flow.
    - `ctrl_separate_s601` is fixed by `--max-extend`: 0.759 → 0.982.
    - `nest_cross_1to1_s502` finishes in 45 s.
    - `nest_cross_3sp_s402` has clean labels, flips and imputed edges. It is still merged at
      λ ≤ 0.1 (2 species for 3, 0.523) and recovered at λ ≥ 0.2 (0.978). **The λ = 0.05–0.5
      plateau no longer holds.**
    - Mean composite over 22 cases: 0.827 at λ = 0, 0.908 at 0.05–0.1, 0.928 at 0.2 and 0.930 at
      0.5. No control or demo case loses against λ = 0. The rule picks 0.5 again because of
      `nest_samehap_3to1_s301`. λ = 0.2 is the smallest value with every real gain.
    - Failures at λ = 0.1: 6 of 14 bench cases (5 at 0.2):
      - `nest_samehap_3to1_s301`, unchanged;
      - `nest_cross_3sp_s402` (λ ≤ 0.1);
      - `nest_cross_3sp_s403`;
      - `nest_cross_1to1_s502` and `s505`, a third low-weight species;
      - `ctrl_dup_s703`, 2 species for 1.
    - `nest_cross_3sp_s403`, like the old `nest_cross_1to1_s501`, `s502` and
      `nest_cross_3sp_s406`, has a subgraph whose fitted depth per copy is half the true one
      (ploidy 4–6 instead of 2). Its segment copy numbers are about 2× the junction-read ones.
      `pipeline.tsv` labels these "RCK" because its stage rule does not check the scale.
    - 4349 s of RCK wall time, 2860 s of it in `nest_cross_3sp_s403`.
    - Follow-up checks (`results/bench_nested/README.md`):
      - s402 switches to the right answer at λ = 0.11.
      - s502 is resolved by the fitted junction reads per copy (5.11), which is wrong for s505
        and s403.
      - s403 still fails with the right scale.
      - s301 and s703 have a near-correct answer in the pool that loses on the ranking. The exact
        decomposition of a noisy integer flow trades one copy of residual against one species.
- Tests under `tests/ecdna/`, with fixtures from demo3 and demo4: subgraph reading, scale,
  flow, decomposition, walks, refit, subgraph inputs, evaluation, recovery matrix
  (`test_recovery.py`), figures (`test_plot.py`), Gurobi environment, and end-to-end demo3 and
  demo4.
- Documentation:
  - `docs/ecDNA.md`: install, inputs, model, outputs, figures, evaluation and the recovery
    matrix, ionic, tests and subgraph-only mode;
  - a link from `README.md`.
- Princeton ionic examples in `docs/examples/`:
  - `ecdna_ionic.sbatch` and `ecdna_demo3_tests.sbatch`;
  - `ecdna_bench.sbatch`: a benchmark array, with at most 2 concurrent tasks for the academic
    Gurobi licence;
  - `collect_bench.py`: summarizes the benchmark and tallies log errors.
- `setup.py`: the `rck.ecdna` package and the `rck-ecdna`, `rck-ecdna-eval` and
  `rck-ecdna-plot` console scripts.
- Result and log files:
  - `results/plan2/`: demo3, and demo4 normal, tumor1 and tumor2;
  - `bench_summary.txt`, `test_run.log` and `test_demo3.log`.
- Validation runs of `rck-ecdna` on new simulations. Each was run from the subgraph folder
  alone, with `--truth-dir`, on a Gurobi academic WLS license. The outputs are kept in
  amplicon-allele-analysis under `simulations/runs/<run>/`, and that repository's
  `CHANGELOG.md` gives the full scores. No RCK code changed.
  - demo5, one simple cycle (`demo5_cycle1_ecdna/rck_out/`): one circular species recovered
    exactly (16.82 copies against 16.0), composite 0.99.
  - demo6, two unnested cycles on chr17 and chr15, given as one union subgraph
    (`demo6_cycle2_ecdna/rck_out_union/`): exactly two circular species, each matching its
    truth species one to one, composites 0.984 and 0.994.
  - demo7, a 3 Mb cycle with a second cycle nested inside it on chr4
    (`demo7_nestedcycles_ecdna/rck_out_union/` and `rck_out_union_totalonly/`): not
    recovered. With default settings the composites are 0.43 and 0. With total-only
    observations (`--min-snps-allele 100000000`) they are 0.51 and 0.33. With
    `--block-flip unoriented` (above; `rck_out_union_blockflip/`) both species are recovered,
    with composites 0.994 and 0.9996.
  - demo8, demo7's two nested cycles on one haplotype at 3:1, with no haplotype crossings
    (`demo8_nestedcycles_ecdna/rck_out_union/`): both species are recovered at rank 1
    (composites 0.985 and 0.989), but only by the edge-traversal tie-break, 0.004 ahead of a
    merged walk with identical copy numbers. This motivated `--step-penalty` and the ambiguity
    report (above).
- Known limitations found by demo7:
  - **Nested species merge** (when allele data cannot separate them). When two species share
    reference segments, the minimum-species decomposition prefers one walk that traverses
    both. On demo7 with total-only observations it returned a 7.54 Mb circle at 8 copies that
    runs through the outer cycle twice and the nested cycle once. The flow itself was close to
    the truth (16–17 and 8–10 junction copies against 16.9 and 9.1). Long-read paths only
    break the tie when some read separates the species, and demo7's 3 multi-junction reads
    each lie inside one species.
    - With `--block-flip unoriented` the demo7 species no longer merge. The two species share
      9 reference segments but only one segment haplotype, because `ecdna1` is on B and
      `ecdna2` on A where they overlap. Once the allele observations are correctly phased, a
      single walk cannot reproduce the haplotype-specific flow more cheaply than two circles.
    - demo8 tests nested species on one haplotype (3:1). Rank 1 was right, but only by a 0.004
      tie-break over a merged walk with identical copy numbers. `--step-penalty` (above) now
      decides such cases explicitly, and the ambiguity report says when the copy numbers
      cannot. Long reads that span the difference remain the only evidence-based remedy.
  - **Allele observations are trusted as phased** (by default). Segment allele copy numbers
    are taken on the true haplotype scale. When `subgraph_construct.py` cannot orient a phase
    block, here because the only pure junctions there attach to segments shared by both
    species, A and B are swapped against the junction labels. On demo7 this made the flow
    objective 39.3 instead of 8.8 and brought in imputed edges. `--block-flip unoriented`
    (above) lets the flow choose a per-block A/B flip for blocks whose `phase_blocks.tsv`
    orientation is 0 or missing. It needs a subgraph with `segment_blocks.tsv`. On demo7 it
    brings the flow objective to 9.13 with no imputed edge. Blocks the subgraph oriented are
    still trusted, unless `--flip-oriented-penalty` is set. It remains opt-in: `off` is the
    default.

- `rck-ecdna` noise-tolerant decomposition, opt-in; defaults and their outputs are unchanged
  (`docs/ecDNA.md`, "Noise-tolerant residual"). The decomposition covers a noisy integer flow at a
  residual cost of 1 per copy, the same as one species, so a 1–3 copy dip from depth bias makes an
  exact fit with a leftover or a merged species win.
  - `--adj-residual-weight` (the existing `DecompositionParams.adj_residual_weight`, now on the
    command line) and `--seg-residual-scale`, a multiplier on the segment residual weight.
  - `--residual-mode relative` with `--residual-sigma` (0.1): every edge's residual weight is
    divided by σ × max(flow, 1). It applies to junctions too, since a segment dip in a balanced
    flow comes with junction dips.
  - `--rank-residual sum`: ranks decompositions on residual + structure cost. The default ranking
    puts residual first, so lower weights alone only change the pool.
  - `--keep-pool` records the new settings; `rck.ecdna.rerank` uses the run's `--rank-residual`.
  - Tests (`tests/ecdna/test_decompose.py`): the residual weights and ranking key without Gurobi,
    and a walk with a 2-copy dip on one loop, an exact two-species fit by default and one walk
    with cheaper residuals (absolute or relative) and the sum ranking.

- **Method-agnostic scorer** `rck.ecdna.aa_eval` (`python -m rck.ecdna.aa_eval`) for the
  comparison with AmpliconArchitect, CoRAL and the original RCK (`docs/ecDNA.md`, section 6,
  "Method-agnostic scorer"). It reads AA-format cycles and graph files only, against simulator
  truth, with `Haplotypes=` optional. No `rck-ecdna` default or output changes.
  - `rck.ecdna.aa_format`: the AA graph and cycles parser, shared with the truth reader. It also
    reads `0+ … 0-` paths and the optional `Haplotypes=`, `Status=`, `Species=` and
    `Alt_segments=` / `Alt_haplotypes=` keys.
  - Graph level, per truth amplicon, at `--graph-tol` 100 bp, against the truth mixture graph of
    amplicon-allele-analysis `simulations/truth_mixture.py`:
    - AA error-free reconstruction;
    - CoRAL breakpoint-graph accuracy;
    - RCK novel-adjacency FNR and FDR;
    - segment copy-number error, allele-specific with `--segment-cn`;
    - infinite-sites violations.
  - Species level: the `rck-ecdna-eval` scores at `--species-tol` 10 kb, unchanged
    (`score_decomposition` and `score_pair` are called as they are), plus:
    - CoRAL's cycle-interval overlap, cyclic LCS and log2 length error, each also in the form
      CoRAL's code computes it;
    - exact cycle at 100 bp;
    - composite-H0: haplotype-blind, without phasing accuracy.
  - Sample level: species-count error, missed and false species, the copy number explained by
    the 1 and 3 heaviest cycles (CoRAL), and the long-read subwalks satisfied.
  - Two modes in every table. `strict` is the method's answer; `--strict-rule` is `status`,
    `cn-frac:F`, `top-cn:X` or `all`, and the default `auto` is `status`, else `cn-frac:0.1`,
    `rck-ecdna`'s own threshold. `best_of` takes, per truth species, every cycle of the answer
    and of `--alt-cycles`, with each metric maximised on its own, as CoRAL scores.
  - `rck.ecdna.aa_export`: writes an `rck-ecdna` answer for the scorer, a run's rank 1 or a kept
    pool's entry. It writes the low_weight status, the other orderings, the first k same-fit
    alternatives (as `diag/grow_rules.py`'s set metric), the graph and allele copy numbers.
    Positions are the walk's + 1, so the scorer reads exactly the steps `rck-ecdna-eval` reads.
  - Tests: `tests/ecdna/test_aa_eval.py` (24, no Gurobi; 27 since `error_free_merged` and the
    CoRAL-code headline).
- **bench_nested: method-agnostic metrics, opt-in** (`collect.py --aa-metrics`). The flag
  re-scores every case and λ through the AA-format path, appends the main columns (`strict.*`,
  `best_of.*`, `aa_path_matches`) to `bench_sweep.tsv` after the existing ones, and writes
  `aa_metrics.tsv` and a README summary between `aa-metrics` markers. The λ rule and the failure
  rule are unchanged. `diag/aa_regression.py` compares the AA path with today's scores field by
  field.
  - On the 40 cases at λ = 0 and 0.2, every composite part, composite, matched species,
    `n_species`, `n_truth` and `merged` is identical at full precision. With `--aa-metrics`, all
    240 (case, λ) match, and the 19 existing columns of `bench_sweep.tsv`, `pipeline.tsv`,
    `rerank_validation.tsv` and the README results section are unchanged.
  - Results (README, "Method-agnostic metrics (2026-11-03)"), at λ = 0.2:
    - best-of mean composite 0.931 against 0.853 strict;
    - composite-H0 0.835;
    - 55 of 72 truth species exact cycles;
    - breakpoint-graph accuracy 1.000, novel-adjacency FDR 0.007;
    - 1 of 42 amplicons AA error-free: all 532 unmatched inferred segment boundaries are
      imputed extremities of the subgraph.
- **`aa_eval`: `error_free_merged`** (graph level, next to `error_free`, which stays the
  headline and is unchanged). Before AA's rule is applied, `merge_benign_splits` merges, in the
  truth and in the method's graph, neighbouring segments joined by a concordant edge with no
  discordant or source edge at the boundary: AA and CoRAL graphs, like `rck-ecdna`'s, also split
  segments at copy-number changes without a breakpoint. New columns `error_free_merged` and
  `n_{true,inferred,false}_segments_merged`, `n_segments_found_merged`; sample level
  `error_free_merged`. On the bench at λ = 0.2, 38 of 42 truth amplicons pass it (1 of 42
  `error_free`). Tests: a benign copy-number split passes it and fails `error_free`, the same
  split with a discordant edge at the boundary fails both, and `merge_benign_splits` on its own.
- **`aa_eval`: `aa_eval.json` `info.mixture_graphs`**, the truth mixture graph files a run was
  scored against.
- **Graph route (Task 2): one AA breakpoint graph for CoRAL 3.1.1, CycleExtractor, original RCK and
  RCK-ecDNA** (`results/graph_route/`, docs/ecDNA.md section 6 "Graph route"):
  - `aa_format.write_graph`, a graph writer both CoRAL 3.1.1 and CycleExtractor 1.0.0 read (no
    blank line, 7- and 6-field lines, `PathConstraint:` / `AmpliconIntervals:` blocks, `#`
    comments, exact floats); `read_graph` also returns path constraints, intervals, comments and
    sequence coverage / size / reads; `coral_path_string` (CoRAL's `path_to_str`) and `chrom_key`
    (CoRAL's chromosome order).
  - `python -m rck.ecdna.aa_export subgraph`: the normalised subgraph with the copy numbers
    `rck-ecdna` uses (junction CN = reads / θ), read paths as CoRAL path constraints, in CoRAL's
    sort order (`coral cycle` re-sorts edges without renumbering path constraints); imputed
    candidates, region-flank joins and the subgraph's ids in `#rck_*` comment lines.
    `aa_export RUN OUT` is unchanged.
  - `python -m rck.ecdna.aa_import` (`rck/ecdna/aa_import.py`): any AA graph back to a subgraph
    folder for `rck-ecdna --inputs subgraph`, with `scale.json` fitted from the graph and
    `--imputed comment|none|regenerate`; `compare_subgraphs` for round trips. Export then import
    is exact on demo5–demo8 and the 32 bench cases.
  - `results/graph_route/check_coral_format.py` (checks graphs with CoRAL 3.1.1's own parser,
    sorter and writer; reference tokens; CycleExtractor with HiGHS), `roundtrip.py` /
    `roundtrip.sbatch` / `roundtrip.list` (the round-trip test; Slurm part prepared).
  - Tests: `tests/ecdna/test_aa_import.py` (20, no Gurobi), reference tokens written by CoRAL in
    `tests/ecdna/data/graph_route/`.
  - `subgraph_io.build_graph` reads an optional `AID` on `extremity` rows (aa_import uses it to keep
    reference-join ids); subgraph_construct.py never writes one, so its subgraphs read as before.
  - `results/graph_route/envs/`: CycleExtractor now runs with the coral-3.1.1 environment's python
    (same numpy, networkx, gurobipy, highspy); the separate `cycle-extractor` environment and its
    lock files are removed.

### Changed

- **Environments, tools and run outputs now live in project space, `/n/fs/yubi-iw`** (2026-11-04;
  home is a 5 GB quota and scratch purges after two weeks). Scripts and docs point there:
  - conda environments: `$HOME/.conda/envs/<env>` -> `/n/fs/yubi-iw/conda/envs/<env>`
    (`gurobi-env`, `coral-3.1.1`) in `bench_nested.sbatch`, `diag/rck_diag.sbatch`,
    `diag/summarize.sh`, `results/step_sweep/sweep.sbatch`, `results/graph_route/roundtrip.sbatch`,
    the graph-route README, NEXT_STEPS and `envs/README.md`, and docs/ecDNA.md;
  - CycleExtractor: `~/tools/CycleExtractor` -> `/n/fs/yubi-iw/tools/CycleExtractor`
    (`check_coral_format.py`'s `CE`);
  - run outputs: `/n/fs/scratch/$USER/{bench_nested,bench_nested_v2,bench_nested_aa,graph_route}`
    -> `/n/fs/yubi-iw/runs/...`: the defaults of `diag/aa_regression.py` (`DEFAULT_OUT`, used by
    `collect.py --aa-metrics`) and `roundtrip.py` (`GRAPH_ROUTE_OUT`), `roundtrip.list`, the
    `diag/rck_noise` LIST files, and the paths in the bench_nested and graph_route READMEs,
    NEXT_STEPS, `comparison.md` and `diag/rck_noise/SUMMARY.md`. `aa_regression.mixture_dir`
    builds a missing truth mixture in place for simulations under `/n/fs/yubi-iw/runs/` too, as
    it did on scratch.
  - Slurm logs and temporary files stay on scratch (`roundtrip.sbatch` and NEXT_STEPS log paths,
    `rck_diag.sbatch`'s submit folder), as does the not-yet-built `bench_nested_sigma0` root.
    Past run records (`diag/rck_noise/runs/*/summary.json`) keep the paths they were run with.
  - Text that said outputs stay on scratch and are purged after two weeks (bench_nested README,
    graph_route README and NEXT_STEPS, `diag/rck_noise/SUMMARY.md`, docs/ecDNA.md) now
    describes project space; `envs/README.md` says the coral-3.1.1 environment is in project space.
  - Outside the repository: the 55 symlinks inside the copied runs
    (`bench_nested_v2/{grow/rck_runs,grow/runs,rck_runs,rebuild/runs}`) still pointed at scratch
    (rsync copies links as they are) and now point at their `/n/fs/yubi-iw/runs/` counterparts
    (undo script `/n/fs/yubi-iw/runs/.undo_symlink_repoint_2026-11-04.sh`). The old
    `~/.conda/envs/{coral-3.1.1,gurobi-env}` and `~/.conda/pkgs` are removed (home: 590 MB ->
    2.5 GB free).
  - Checked with the new environments, without Slurm or solves: CycleExtractor's README sample
    (HiGHS) gives the same cycle, `2+,4+,6+` up to rotation (the start segment varies between runs
    in the old environment too), `Copy_count=85.18981770205714`; `coral cycle --help`; `import
    gurobipy` (13.0.3) and an empty model read the WLS licence; `tests.ecdna.test_aa_import`,
    `test_aa_eval`, `test_subgraph_io` (60 tests) pass; `collect.py` and `collect.py --aa-metrics`
    (to a separate `BENCH_NESTED_RESULTS`) reproduce `bench_sweep.tsv`, `aa_metrics.tsv`,
    `pipeline.tsv`, `rerank_validation.tsv` and the README results byte for byte, and
    `roundtrip.py compare` reproduces `roundtrip_compare.tsv` / `.txt`, reading no file on
    scratch; the sbatch scripts resolve `python` to the new gurobi-env.

- **`aa_eval`: CoRAL's code versions are the headline cycle metrics.** `cycle_overlap_coral_code`
  and `cyclic_lcs_coral_code` (what produced CoRAL's published numbers; the code's LCS
  normaliser penalises merged walks) now come first in `aa_species_metrics.tsv`. The
  supplementary-methods versions are renamed `cycle_overlap` → `cycle_overlap_supp` and
  `cyclic_lcs` → `cyclic_lcs_supp` (values unchanged). Sample level: `mean_cycle_overlap` and
  `mean_cyclic_lcs` become `mean_cycle_overlap_supp` and `mean_cyclic_lcs_supp`, and
  `mean_cycle_overlap_coral_code` and `mean_cyclic_lcs_coral_code` are added. `collect.py
  --aa-metrics` reports the code versions first and the supplement's as "secondary"; at
  λ = 0.2 the code's cyclic LCS is 0.871 strict and 0.958 best-of (supplement 0.901 and 0.986).
  Re-scoring all 240 (case, λ) left every other column of `aa_metrics.tsv` and
  `bench_sweep.tsv` identical.
- **bench_nested: demo5–demo8 truth mixtures from the repository.** `diag/aa_regression.py`
  uses `truth/aa/mixture_amplicons.tsv` wherever it exists, so demo5–demo8 read
  amplicon-allele-analysis's tracked copies instead of
  `/n/fs/scratch/$USER/bench_nested_aa/mixtures/` (demo4's is still built there). The copies are
  byte-identical, and the re-scored demo cases changed only in their paths.
- `evaluate.read_aa_cycles` parses with `rck.ecdna.aa_format`. Its output is unchanged; an
  unlisted segment or an AA path (`0+ … 0-`) still raises `ValueError`, with a reworded message.

- `rck_run.py`: `MIPGapAbs` is set only from the new `--run-g-mip-gap-abs` option. It was
  previously set to the relative `--run-g-mip-gap` value, which forced near-exact
  optimality.

- `results/bench_nested/`: **bench_nested_v2 replaces the first benchmark** (2026-10-29).
  - `bench_sweep.tsv`, `pipeline.tsv` and `rerank_validation.tsv` are v2's (Slurm 32000320), and
    `comparison.md` and `comparison.tsv` (v1 against v2) are now tracked here. The configs and
    subgraphs are in amplicon-allele-analysis. `runs` points at
    `/n/fs/scratch/$USER/bench_nested_v2/rck_runs`. The first benchmark's numbers stay in
    `README.md`, and its tables in git history (9b74b87).
  - Scratch purges files older than two weeks. The remaining run outputs (`runs`, `rck_runs`,
    `diag`) are not kept anywhere else; the v2 configs, subgraphs and tables are also in
    `~/backups/bench_nested_v2_2026-10-29.tar.gz` on ionic.
  - `collect.py`:
    - a `scale` failure stage, between `subgraph` and `RCK`: depth per copy more than 30% from the
      cases' median, or segment and junction scales against the truth more than 1.5× apart.
      v2 `nest_cross_3sp_s403` is now `scale`. `pipeline.tsv` gains the scale columns;
    - the λ rule: a larger λ is preferred only if at least 2 cases gain more than 0.005 from it,
      so one case at the grid edge can no longer pick λ. v2 now gives λ = 0.2 (was 0.5) and the
      first benchmark 0.05 (was 0.5). `--lambda` forces the reported λ.
  - `diag/scale_ratio.py`: `scales(run)` can be imported (`collect.py` uses it); the output is the same.
  - `diag/rck_diag.sbatch`: `DIAG_POOL_LAMBDA` sets the pool run's λ and `DIAG_LAMBDA0=0` skips the
    λ = 0 run. `diag/summarize.sh` skips the missing λ = 0 run.
  - **Step 2 diagnostic of the noise-tolerant decomposition** (Slurm 32025976; `diag/rck_noise/`:
    `SUMMARY.md`, `noise_summary.tsv`, the LIST files, and each run's small outputs in `runs/`).
    - The four decomposition failures (s301, s703, s505, s403 at ploidy 2), each at residual
      weight 1, 0.5 and 0.25 and relative σ 0.1, at λ 0.2. Also s703 and s301 on the path-fixed
      subgraphs at `--path-weight` 0.5 and 2 and `--path-fraction 0.9`.
    - **No setting passes**: s301 and s703 keep their rank 1 (0.693, 0.712) at every setting, and
      s505 stays at 0.765 (the best in its pool).
    - The 0.97 answers of s301 and s703 lose on leftover species counted as significant, not on
      residual; both they and rank 1 fit the integer flow exactly. Refitted to the observations,
      every s703 pool solution fits identically (`data_objective` 5.097) and s301's differ by
      0.002, so neither a refit nor decomposing against the observations would separate them. The
      reads are already explained by the wrong answers.
    - s403 at relative σ 0.1 stays at 0.545 (pool best 0.729, not picked).
    - **Still running when committed (2026-10-29 15:55):** array 32025976 tasks 15
      (`nest_cross_3sp_s403_w0.25`) and 18 (`nest_cross_1to1_s505_w0.5`), with 19–21 (s403 at w0.5
      and w1, s505 at w1) pending; due by about 17:00. Their rows in `SUMMARY.md` and
      `noise_summary.tsv` say pending. Check them with `sacct -j 32025976`, then rerun
      `noise_summary.py --archive` and complete `SUMMARY.md` (`NEXT_STEPS.md`, "Check first").
    - **Completed 2026-10-31**: all 22 tasks COMPLETED. The 5 runs pending at commit are archived
      in `diag/rck_noise/runs/`, and `SUMMARY.md` is complete. s505 stays at 0.765 (the best in its
      pool) and s403 at 0.545 (pool best 0.729, not picked) at every residual setting: no setting
      changes rank 1 or adds a better pool answer. `nest_cross_3sp_s403_w0.25` (task 15) hit
      Gurobi's `--time-limit 1800` in its pool search (optimum proven, pool possibly incomplete),
      so its row is described, not filled. `SUMMARY.md` also corrects the s301 species counts the
      ranking sees (2 against 3 significant, from the integer weights).
  - `diag/noise_summary.py`: per case and setting, the rank-1 answer, the best mean composite in
    the pool and whether rank 1 picked it, the inexact pool fits and the runtime; `--archive` copies
    each run's small outputs into the repository.
  - **Assessment of the next steps** (2026-10-31; `NEXT_STEPS.md`, "Assessment (2026-10-31)", with
    a revised plan; the 2026-10-29 plan is kept). Offline only, with no Gurobi runs and no RCK
    default changed.
    - `diag/leftover_recount.py` (new): re-flags every saved pool solution's species under a
      leftover rule and re-ranks with `rank_pool`. Rules: the global `--min-species-weight-frac`;
      an overlap-relative rule (σ × the heaviest overlapping species; overlap on phased or
      unphased segments, any or > 50% of the lighter species' length); a junction-novelty rule.
      The default rule reproduces every case's rank 1 (mean 0.9283). Output:
      `diag/rck_noise/leftover_recount.tsv`, 22 bench pools and 5 diagnostic pools.
      - **The step 2b leftover prior does not help.** At σ 0.15 it is case for case the global
        0.15 threshold (mean 0.928 → 0.910; s703 fixed; s502 and s505 broken), because every
        leftover overlaps the heaviest species. At σ 0.2 the mean is 0.894.
      - s505 and s502 break through the 1-copy floor of max(1, 0.1 × heaviest).
      - The junction-novelty variant fixes s301 and s703 but breaks s502, s505 and s402 (0.916).
    - `diag/obs_vs_truth.py` (new): segment-depth and junction-read errors against the truth copy
      numbers, after the flow's background (`diag/rck_noise/obs_vs_truth.tsv`). Junction reads are
      more accurate than depth in s703 and s502, not in general.
    - **Findings.**
      - RCK's same-fit alternatives already contain s703's and s301's right answers (first and
        second alternative); `pool_keeps` drops s301's from the ranked list at λ 0.2.
      - The paired normal's depth per copy picks the right ploidy class in 14 of 14 bench
        subgraphs (s403's class 6 is at 0.48× the normal).
      - The composite scores s403's 2× scale error (0.592) above the right scale (0.545).
      - s703's answers can only be told apart by reads over ≈ 63 kb; the simulated reads average
        15 kb.
    - **Revised plan.**
      1. Ambiguity reporting plus a set metric.
      2. A σ = 0 (no depth bias) control run.
      3. A re-baseline with the path fix and a paired-normal ploidy rule.
      4. More seeds and a light-nested scenario before tuning any prior.
      - Step 2b is dropped and θ deferred.
      - The run LISTs and sbatch commands are written out; nothing was submitted.
    - **Agreed order and step c** (2026-10-31, later; `NEXT_STEPS.md`, end of the Assessment):
      - Order: re-baseline the 22 cases, then more seeds; the σ = 0 control is optional; ambiguity
        reporting comes later.
      - **Correction:** the paired-normal ploidy rule is a hidden diploid prior. Across ploidy
        classes, depth per copy × cell ploidy is the tumour's mean depth, so a matched normal
        cannot choose the class. The bench subgraphs use an explicit `--rck-ploidy 2` instead.
      - The 14 bench subgraphs were rebuilt on scratch: graphs unchanged, more paths, and s403 at
        ploidy 2.
      - RCK on all 22 cases was submitted as Slurm 32041913 (outputs in
        `/n/fs/scratch/$USER/bench_nested_v2/rebuild/rck_runs`); this folder's tables are
        unchanged until it is collected.
  - **Step c done: the re-baseline is the new baseline** (2026-11-01; `README.md`, "Re-baseline
    (2026-10-31)"; `NEXT_STEPS.md`, "Step c, result").
    - Slurm 32041913 completed all 22 tasks (2762 s of RCK wall time).
    - Rank 1 at λ = 0.2 is unchanged in 21 of 22 cases; the path fix moves no composite.
    - `nest_cross_3sp_s403` goes from 0.592 to 0.545. Its scale is now right, and `pipeline.tsv`
      stage `scale` becomes `RCK`. Mean at λ = 0.2: 0.9283 → 0.9262. λ = 0.2 is still recommended.
    - s403's rank 1 is the 2026-10-27 ploidy-2 diagnostic's exactly: the same observations and
      flow, and the same rank-1 walks except the path columns. Its pool best is again 0.729 and
      is not picked, so its failure is in the decomposition, not the scale.
    - `bench_sweep.tsv`, `pipeline.tsv` and `rerank_validation.tsv` are the re-baseline's, copied
      from `rebuild_2026-10-31/` (`rerank_validation.tsv` is unchanged). `collect.py` reproduces
      them byte for byte from the adopted subgraphs. v2's tables are in git history (7778991),
      and its numbers in `README.md` ("v2 results before the re-baseline").
    - `rebuild_2026-10-31/comparison_vs_v2.txt` is now the row-by-row v2 comparison. 7778991
      committed a copy of `collect.log` under that name.
  - **Step d prepared** (`NEXT_STEPS.md`, "Step d, prepared"):
    - more seeds for nest_samehap_3to1, nest_cross_1to1 and ctrl_dup (6 accepted each), and the
      new nest_samehap_6to1 scenario (amplicon-allele-analysis CHANGELOG);
    - root `/n/fs/scratch/$USER/bench_nested_v2/grow/` links the existing simulations and RCK runs,
      so only the new seeds are simulated and run;
    - the commands, with the updated σ = 0 control (step b, at ploidy 2), are written out.
  - **Step d1 done, d2 submitted** (2026-11-01; `NEXT_STEPS.md`, "Checking d1, result" and
    "Step d2, submitted"):
    - d1 (Slurm 32049944) completed: 6 accepted seeds each for the four grown scenarios, 18 new
      subgraphs, no subgraph-stage failure, all at `--rck-ploidy 2`.
    - d2, RCK on the 18 new cases only, submitted as Slurm array 32050408 (indices
      14–17, 22–25, 30–39, `%2`). The 22 existing cases are links to `rebuild/rck_runs` and count
      as done, so none reruns. Outputs in `/n/fs/scratch/$USER/bench_nested_v2/grow/rck_runs`;
      this folder's tables are unchanged until it is collected and reviewed.
  - **Step d done; the grown bench adopted** (2026-11-01; `NEXT_STEPS.md`, "Step d2, result"):
    - Slurm 32050408: 18 of 18 COMPLETED. Over 40 cases (32 bench, 8 regression), the mean
      composite at λ = 0.2 is 0.8531 (0.9262 over the 22-case re-baseline). λ = 0.2 stays the
      recommendation. Failures: 17 of 32 bench cases.
    - The 22 existing cases are unchanged: `bench_sweep.tsv` and `rerank_validation.tsv` are
      identical row for row. In `pipeline.tsv` only `depth_vs_median` moves (≤ 0.005), because
      the bench median depth per copy now spans 32 cases.
    - `bench_sweep.tsv`, `pipeline.tsv` and `rerank_validation.tsv` are the grown bench's (the
      run record is in `grow_2026-11-01/`). `README.md`: 8 scenarios, 6 seeds for the 4 grown
      ones, 40 cases; the 22-case results are kept. The local `runs` link points at
      `bench_nested_v2/grow/rck_runs`.
  - **Failure analysis of the grown bench and a revised plan** (`NEXT_STEPS.md`, "Grown bench:
    failure analysis and revised plan (2026-11-01)"). It uses the saved pools only; no new runs.
    - New offline scripts `diag/grow_failures.py` (per-case anatomy of rank 1 against the pool's
      best) and `diag/grow_rules.py` (re-ranks all pools under candidate rules with
      `rank_pool`). Outputs are in `grow_2026-11-01/`.
    - nest_samehap_6to1: the nested species is split into 2–3 sub-threshold low_weight walks,
      so rank 1 has one significant species; the data cannot tell this from the right answer.
    - ctrl_dup s708 and s710: the truth traverses a segment 4 times, above `--max-multiplicity`
      3, so it cannot be expressed as one species.
    - The output re-flags low_weight on refitted weights, while the ranking uses integer
      weights. They disagree in 8 of 40 rank-1 answers, both ways.
    - Plan: ambiguity reporting first (set metric 0.923 at k = 2, no case down), then a
      `--max-multiplicity 4` check, then the extended σ = 0 control.
  - `.gitignore`: `/runs` instead of `runs`, so that only the runs symlink is ignored and not
    `diag/rck_noise/runs/`.

### Fixed

- `evaluate`: truth haplotypes of simulator runs in which the walk traverses one junction with
  several phasings. `truth/adjacency_clone.tsv` keeps one `PHASING` per junction, its most
  common, and it was applied to every traversal. amplicon-allele-analysis runs made with per-step
  haplotype labels (every run before its `hap_labels: inherit`) can carry, for example, AA on one
  copy of a junction and BB on its duplicate.
  - `truth_phasings` reads every phasing with copies from `truth/rck/rck.acnt.tsv`, and
    `species_junctions` keeps them all (`PHASING` alone without that file).
  - Each traversal takes one of its junction's phasings: the assignment with the fewest
    haplotype conflicts along the walk, then the one that uses the most of the junction's
    phasings. A crossing between reference neighbours may also be traversed on one haplotype,
    as a plain reference join the truth does not list.
  - Of the 22 bench_nested and step-sweep truth walks, only `ctrl_dup_s702` changes: 0.995 Mb of
    its 5.37 Mb (the reverse copy of chr16:57.78–58.77 Mb) was on B, and is on A. Re-scored
    against it, the committed λ = 0 result has phasing accuracy 1.0 (0.815 before) and haplotype
    Jaccard 1.0 (0.687), composite 0.781 (0.744). Its remaining loss is a linear species
    without the true order.
  - Tests in `test_evaluate.py`; documented in `docs/ecDNA.md` (section 6).
- `rck_run.py`: `--scnb` copy-number boundaries are read from the `--scnb` file. They were
  previously read from the `--scnt` path, so `--scnb` was silently ignored.
- ecDNA depth-scale mode estimate: wrong results when the smoothing kernel was longer than the
  histogram.
- A circle traversed several times is reported as one circle at the summed weight (gcd
  normalization).
- Decompositions are ranked by significant species, so one molecule plus small leftovers
  ranks above the molecule split in two.
- Gurobi models are disposed after solving, and test temp folders work on a full `/tmp` and
  on NFS.
- `evaluation.json` is strict JSON: missing scores in `top_decomposition` are `null`, where they
  were previously `NaN`.
- `--truth-dir` evaluation works for any simulated run and subgraph: several species,
  deleted segments and per-sample copy numbers.

[Unreleased]: https://github.com/RunpengLuo/RCK/compare/rck/version-fix...yubi_dev
