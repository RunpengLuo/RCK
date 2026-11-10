# Allele-specific ecDNA reconstruction: pipeline overview

This pipeline reconstructs extrachromosomal DNA (ecDNA) species from bulk tumour sequencing data,
resolved by haplotype. Each species is reported as an ordered walk over phased segments, with a
copy number. The pipeline has three stages, plus an evaluation and benchmarking layer that compares
it with AmpliconArchitect (AA), CoRAL, CycleExtractor and original RCK.

```
 Stage 1  simulator (amplicon-allele-analysis simulations/)
          clone tree + somatic CNAs/SVs + phased germline SNPs + nested ecDNA species
          ├─► observations: bin depth, allele counts, novel-adjacency / extremity / path read counts
          └─► truth: RCK karyotypes, AA graphs and cycles, truth mixture graph per amplicon
                 │
 Stage 2  subgraph construction (simulations/subgraph_construct.py, rck_inputs.py)
          amplicon discovery → breakpoint subgraph → allele annotation → γ, θ, CN columns, paths
                 │
 Stage 3  rck-ecdna (RCK rck/ecdna/)
          two-layer flow MILP → minimum-species decomposition → refit → walks and orderings
                 │
 Evaluation  rck-ecdna-eval, rck.ecdna.aa_eval, aa_export / aa_import, coral_flank
          benchmarks: results/bench_nested/, results/graph_route/, results/step_sweep/
```

## Repositories and commits

| Part | Repository, branch | Commit documented |
|---|---|---|
| Stage 1, Stage 2 | [RunpengLuo/amplicon-allele-analysis](https://github.com/RunpengLuo/amplicon-allele-analysis/tree/yubi_dev) `yubi_dev` (diffed against `main`) | eea71aa |
| Stage 3, evaluation, benchmarks | [RunpengLuo/RCK](https://github.com/RunpengLuo/RCK/tree/yubi_dev) `yubi_dev` (diffed against `rck/version-fix`) | b4ff12d |
| Original RCK, ported to current libraries | `rck/version-fix` = upstream + 0b41976 | 0b41976 |
| Original RCK | [aganezov/RCK](https://github.com/aganezov/RCK) `master` (the fork's `master` is identical) | 4c6c288 |

Both repositories' `CHANGELOG.md` files record each change and the reasoning behind it, in date
order. These docs describe the methods as they now stand.

## Documents

| Document | What it covers |
|---|---|
| [Stage 1: simulator](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/README.md) | Generative model S0–S11, truth files, truth mixture graph, benchmark scenarios and seed checks |
| [Stage 2: subgraph construction](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/subgraph_construct.md) | Discovery, graph construction, allele annotation, region growth, γ/θ fitting, paths, checks against the truth |
| [Stage 3: `rck-ecdna`](../ecDNA.md) | Normalised graph, observations, flow MILP, decomposition, ranking and ambiguity, refit, walks |
| [Evaluation and metrics](evaluation.md) | Every metric (`rck-ecdna-eval`, `aa_eval`), translators to and from AA format, the CoRAL flank fix |
| [bench_nested](../../results/bench_nested/README.md) | Nested-ecDNA benchmark: cases, λ rule, failure stages, results |
| [graph_route](../../results/graph_route/README.md) | Head-to-head of CoRAL, CycleExtractor and RCK-ecDNA on the same AA graphs |
| [step_sweep](../../results/step_sweep/README.md) | The `--step-penalty` sweep on demo3–demo8 |
| [Glossary](glossary.md) | Code names and the published terms they correspond to |
| [References](references.md) | Verified bibliography and the list of unverified items |

Design history, kept unchanged:
[simulator_design.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/simulator_design.md)
(what AA and CoRAL simulate, and the plan for the extension) and
[simulation.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/simulation.md)
(early notes on the CoRAL and AA generative processes). The decision logs are
`results/bench_nested/NEXT_STEPS.md` and `results/graph_route/NEXT_STEPS.md`.
The text of every doc before this rewrite is kept beside it as `*_old.md`.

## How every step is classified

Each algorithmic step falls in exactly one category. Every stage doc gives it in a
**Provenance** line under the step's heading.

- **(a) existing**: taken from a named published method (RCK, AA, CoRAL, ecSimulator, MASCoTE,
  HATCHet, ASCAT, CycleExtractor), with the paper and/or the function it comes from.
- **(b) patch**: needed to run on current library versions or on Princeton ionic (Gurobi WLS
  licence, conda, Slurm, storage paths). A patch changes no result.
- **(c) new or modified**: written or changed for this pipeline. For each such step, the stage doc
  gives the failure that motivated it, its measured effect, and any published work behind the
  choice. Two `rck_run.py` changes are labelled "(c) fix to original RCK": they change results
  when the affected option is used.

Simulator steps that were already on amplicon-allele-analysis `main` are classified one by one,
with the commit (`main@…` or `yubi_dev@…`). Subgraph construction ports Yubi Mamiya's unpublished
Colab notebooks (`colo320_graphconstruction_run1_chrdiscovery`, `…_multichr`), which count as (c).

## Provenance table

Commits are on `yubi_dev` unless marked `main@` or `version-fix@`. "S2 §4.3" means section 4.3 of
the Stage 2 doc.

### Stage 1: simulator (amplicon-allele-analysis `simulations/`)

| Step | Cat. | Source | Code | Commit | Doc |
|---|---|---|---|---|---|
| S0 genomic grid within the whitelisted regions | (c) | — | `src/grid.py::build_grid` | main@2da6045 | S1 §S0 |
| S0 heterozygous germline SNPs | (c) modified | MASCoTE `buildHaplotypes` | `src/grid.py::_place_snps` | main@2da6045 | S1 §S0 |
| S0 phase blocks and switch errors | (c) | — | `src/grid.py::_phase` | main@2da6045 | S1 §S0 |
| S1 clone tree from a Prüfer sequence | (a) | MASCoTE `RandomTree.buildRandom` | `src/clones.py::prufer_tree` | main@2da6045 | S1 §S1 |
| S1 branch lengths and classes | (c) | — | `src/clones.py::build_clones` | main@2da6045 | S1 §S1 |
| S2 somatic event classes (WGD, chromosome, arm, focal) plus balanced SVs | (c) modified | MASCoTE; ecSimulator | `src/clones.py::_apply` | main@2da6045 | S1 §S2 |
| S2 Poisson event-rate matrix | (c) | — | `src/clones.py::build_clones` | main@2da6045 | S1 §S2 |
| S2 infinite-sites rejection | (c) modified | RCK | `src/clones.py::_blocked`, `::_take` | main@2da6045 | S1 §S2 |
| S3 ecDNA interval sizes and placement, breakpoints, rearrangements | (a) | ecSimulator `compute_interval_sizes`, `select_interval_regions`, `assign_bps`, `conduct_EC_SV` | `src/structures.py::interval_sizes`, `::select_intervals`, `::assign_breakpoints`, `::rearrange` | main@2da6045 | S1 §S3 |
| S3 nested species (`nest_species`) | (c) | nested ecDNA as in Decoil | `src/structures.py::select_nested_intervals` | df53d6b | S1 §S3 |
| S3 `per_step` haplotype labels | (c) | — | `src/structures.py::label_haplotypes` | main@2da6045 | S1 §S3 |
| S3 inherited haplotype labels (`hap_labels: inherit`) | (c) | — | `src/structures.py::inherit_haplotypes` | 618697f | S1 §S3 |
| S3 nested species on one haplotype | (c) | — | `src/structures.py::build_structures` | b566528 | S1 §S3 |
| S4 common segment refinement | (c) | RCK `docs/Segments.md` | `src/segments.py::build_segments` | main@2da6045 | S1 §S4 |
| S5 clone proportions | (a) | HATCHet / MASCoTE | `src/mixture.py::_draw_proportions` | main@2da6045 | S1 §S5 |
| S5 ecDNA abundance ψ, φ; per-cell negative binomial | (c) | motivated by Lange et al. 2022 | `src/mixture.py::_abundance` | main@2da6045 | S1 §S5 |
| S6 genome-length correction from clone to read proportions | (a) | MASCoTE `MixBAMs` | `src/mixture.py::build_mixture` | main@2da6045 | S1 §S6 |
| S7 bin read counts | (c) modified | CoRAL `BreakpointGraph.compute_cn_lr` | `src/emit_counts.py::_bin_reads` | main@2da6045 | S1 §S7 |
| S8 allele counts | (c) modified | HATCHet | `src/emit_counts.py::_snp_counts` | main@2da6045 | S1 §S8 |
| S9 Poisson novel-adjacency counts | (a) | CoRAL `compute_cn_lr` | `src/emit_counts.py::_adjacency_tables` | main@2da6045 | S1 §S9 |
| S9 θ_bp, allele split, long-read paths | (c) | — | `src/emit_counts.py::build_counts`, `::_path_table` | main@2da6045 | S1 §S9 |
| S11 figures | (c) | cnplot | `src/plots.py` | main@2da6045 | S1 §S11 |
| Self-checks | (c) | — | `src/checks.py::run_checks` | main@2da6045 | S1 §3 |
| RCK-format truth | (a) | RCK formats | `src/writers.py::_write_rck` | main@2da6045 | S1 §4.1 |
| RCK 1-based adjacency coordinates | (c) fix | RCK convention | `src/writers.py::_write_rck` | 8cd63f4 | S1 §4.1 |
| AA-format truth per species | (a) | ecSimulator `ecSim_IO.write_bpg_file`, `write_cycles_file` | `src/writers.py::_write_aa` | main@2da6045 | S1 §4.2 |
| ecSimulator path inside the repository | (b) | — | `src/writers.py` `ECSIM_SRC` | 5bdb8e5 | S1 §4.2 |
| Truth mixture graph per amplicon | (c) | ecSimulator amplicon convention | `truth_mixture.py::build`, `::merge_silent` | 40c6c73 | S1 §4.3 |
| Mixture reads, coverage, CoRAL path constraints | (c) | CoRAL `chr_sort_key`, `path_to_str` | `truth_mixture.py::coral_order`, `::write_amplicon` | 5e2da37 | S1 §4.3 |
| Benchmark and nested scenarios | (c) | — | `make_scenarios.py` | 785159e, d1a1dce, f1b09b6 | S1 §5 |
| Seed acceptance checks | (c) | — | `bench_nested.py::check_seed` | d1a1dce, 618697f | S1 §5.3 |

### Stage 2: subgraph construction (amplicon-allele-analysis `simulations/`)

| Step | Cat. | Source | Code | Commit | Doc |
|---|---|---|---|---|---|
| Discovery: junction clustering, snapping, path links, ranking | (c) | notebook `…_chrdiscovery` | `subgraph_construct.py::discover_amplicons`, `::load_paths` | aacc399 | S2 §1 |
| Depth signal, bin aggregation, normal correction | (c) | — | `::load_windows`, `::aggregate_windows` | aacc399 | S2 §2.1 |
| Amplification threshold (median + k·MAD) | (c) | robust scale as in Rousseeuw & Croux 1993 | `::build_full_graph` | aacc399 | S2 §2.2 |
| Depth segmentation (binary segmentation, penalised L2) | (a) technique in (c) code | ruptures; Scott & Knott 1974; Yao 1988 | `::build_full_graph` | aacc399 | S2 §2.3 |
| Imputed extremities and candidate novel adjacencies | (c) | notebook `…_multichr` | `::build_full_graph` | aacc399 | S2 §2.4 |
| Full breakpoint graph | (c) | RCK terms | `::build_full_graph` | aacc399 | S2 §2.5 |
| Node statistics | (c) | — | `::node_statistics` | aacc399 | S2 §2.6 |
| Segment BAF (Jeffreys interval) and haplotype call | (c) | Brown et al. 2001; ASCAT/HATCHet SNP pooling | `::baf_stats`, `::hap_call` | aacc399 | S2 §3.2 |
| Junction read label | (c) | — | `::annotate_alleles` | aacc399 | S2 §3.3 |
| Phase-block orientation by pure-junction votes | (c) | [citation needed] | `::annotate_alleles`, `::orient_sum` | aacc399 | S2 §3.4 |
| Full `phase_blocks.tsv`, `segment_blocks.tsv` | (c) | — | `::phase_block_table`, `::segment_block_table` | 8720eaf | S2 §3.5 |
| Seeds, `--observed-seeds-only` | (c) | — | `::trim_graph` | aacc399 | S2 §4.1 |
| Region growth along amplified runs (`--max-extend`) | (c) | — | `::amplified_runs`, `::grow_interval` | 27f6232 | S2 §4.3 |
| Region merge and trimming | (c) | — | `::trim_graph` | aacc399 | S2 §4.4 |
| Purity and ploidy grid fit, giving γ | (c) modified | ASCAT; HATCHet; PELT | `rck_inputs.py::_fit`, `::estimate_scales` | bf79d23 | S2 §5.1 |
| θ from background extremities | (c) | — | `rck_inputs.py::junction_scale` | bf79d23 | S2 §5.2 |
| Observed copy-number columns | (c) | — | `rck_inputs.py::span_cn`, `::annotate_edges` | bf79d23 | S2 §5.3 |
| `cn_profile.tsv`, flank copy numbers | (c) | — | `rck_inputs.py::cn_profile` | bf79d23 | S2 §5.4 |
| Long-read paths, `clip_path` | (c) | — | `rck_inputs.py::subgraph_paths`, `::clip_path` | bf79d23, e3fcba5 | S2 §5.5 |
| `verify_subgraph.py` checks, every truth phasing, summed-haplotype segment calls | (c) | — | `verify_subgraph.py::score_folder`, `::read_label` | 517e4ab, 618697f, e3fcba5 | S2 §6.2 |
| `check_segment_blocks.py` | (c) | — | `check_segment_blocks.py::check` | 8720eaf | S2 §6.3 |

### Stage 3: `rck-ecdna` (RCK `rck/ecdna/`)

| Step | Cat. | Source | Code | Commit | Doc |
|---|---|---|---|---|---|
| DIAG, haplotype labellings, haplotype-aware balance, ≤ 1 labelling per novel adjacency | (a) | RCK `OptModelMultiClone` | `ilp_circular.py::OptModelCircular._balance` | 7d8dff9 | S3 §4.1 |
| Graph normalisation (phantoms, parallel edges, terminals, regions) | (c) | — | `subgraph_io.py::build_graph` | c257e54 | S3 §2 |
| Labellings from `junction_hap` | (c) | — | `subgraph_io.py::_junction_phasings` | c257e54 | S3 §2.5 |
| γ and θ estimation (mode, ratio, HATCHet input) | (c) | reads HATCHet output | `scale.py::estimate_depth_scale`, `::estimate_junction_scale` | c257e54, 41f5804 | S3 §3.1–3.2 |
| Allele split of segment copy number | (a) | ASCAT / HATCHet | `scale.py::build_observations` | c257e54 | S3 §3.3 |
| Subgraph-only input mode | (c) | — | `scale.py::load_subgraph_scale` | a987ffa | S3 §3 |
| Two-layer flow MILP (ecDNA + background) | (c) modified | RCK | `ilp_circular.py::OptModelCircular` | 7d8dff9 | S3 §4 |
| Terminal slack for the background layer | (a) | RCK telomeres | `::_balance` | 7d8dff9 | S3 §4.3 |
| Imputed edges, objective | (c) | — | `::_objective` | 7d8dff9 | S3 §4.5 |
| Linear ends | (c) | AA / CoRAL linear amplicons | `::linear_end_nodes` | 10416aa | S3 §4.4 |
| Phase-block flips | (c) modified | RCK per-fragment allele flip | `::_piece_terms`, `::block_report` | a422c39 | S3 §4.6 |
| Infinite sites (`--enforce-is`) | (a) | RCK | `::_infinite_sites` | 7d8dff9 | S3 §4.3 |
| Decomposition MILP | (c) modified | CoRAL `minimize_cycles` | `decompose.py::CycleDecomposer.build` | 3c82929 | S3 §5.1 |
| Weight × multiplicity linearisation | (a) | Glover 1975 | same | 3c82929 | S3 §5.1 |
| Lazy connectivity cuts | (a) | Dantzig et al. 1954 | `::_callback` | 3c82929 | S3 §5.2 |
| Noise-tolerant residual | (c) | — | `::_omega` | 2fdf73d | S3 §5.3 |
| Long-read constraints | (c) modified | CoRAL subpath constraints | `::read_constraints` | c5e0b9c | S3 §5.4 |
| Solution pool and no-good cuts | (a) | Gurobi | `::solve` | 3c82929 | S3 §5.5 |
| gcd normalisation, `low_weight`, ranking by significant species | (c) | — | `::_flag_low_weight`, `::ranking_key` | 41f5804, f84faff | S3 §5.6–5.7 |
| Step penalty λ (ranking only), pool pruning | (c) | — | `::ranking_key`, `::prune_pool` | 90ac98a | S3 §5.7 |
| Ambiguity report | (c) | — | `::ambiguity_summary` | 90ac98a | S3 §5.8 |
| Offline re-ranking | (c) | — | `rerank.py`, `::rank_pool` | 949f796 | S3 §5.9 |
| Amplicon call | (c) | — | `::classify` | 10416aa | S3 §5.10 |
| Refit LP | (c) | HiGHS via SciPy | `refit.py::refit_weights` | c5e0b9c | S3 §6 |
| Alternating Euler walk | (a) | Hierholzer & Wiener 1873 | `walk.py::alternating_circuit` | 3c82929 | S3 §7 |
| Ordering enumeration and long-read ranking | (c) | — | `walk.py::enumerate_circuits`, `::path_support` | 3c82929 | S3 §7 |
| AA-format outputs (`Haplotypes=` is new) | (a) | AA | `output.py` | 3c82929 | S3 App. D |
| Figures | (c) | CoRAL-style, cnplot | `plot.py` | 4e53747 | S3 App. E |
| `--scnb` path and `MIPGapAbs` in `rck_run.py` | (c) fix to original RCK | — | `rck/rck_run.py` | b26bb96 | S3 §8 |

### Evaluation and benchmarks

| Step | Cat. | Source | Code | Commit | Doc |
|---|---|---|---|---|---|
| Truth walks, haplotypes (several phasings), copy number, scope | (c) | AA format | `evaluate.py::truth_from_simulation`, `::_label_haplotypes` | c144584, 184ec95, 3ba7627 | E §2 |
| Jaccard, junction and phasing scores, order, CN error, matching | (c) | Jaccard 1912; CoRAL breakpoint matching | `evaluate.py::score_pair`, `::junction_scores` | c144584, 4e53747 | E §3.1–3.6 |
| Composite and recovery matrix | (c) | — | `evaluate.py::composite_score`, `::recovery_matrix` | 4e53747 | E §3.7–3.8 |
| `bp_accuracy` | (a) | CoRAL | `aa_eval.py::graph_metrics` | e4c2409 | E §4.3 |
| Novel-adjacency FNR / FDR | (a) | RCK | `aa_eval.py::graph_metrics` | e4c2409 | E §4.3 |
| `error_free` | (a) | AA | `aa_eval.py::graph_metrics` | e4c2409 | E §4.3 |
| `error_free_merged` | (c) modified | AA | `aa_eval.py::merge_benign_splits` | 41686bc | E §4.3 |
| Cycle overlap, cyclic LCS (code and supplement versions), length error | (a) | CoRAL code and supplement | `aa_eval.py::cycle_overlap_coral_code`, `::cyclic_lcs`, `::length_error` | e4c2409, 41686bc | E §4.4 |
| Exact cycle, composite-H0 | (c) | — | `aa_eval.py::_exact`, `::composite_h0` | e4c2409 | E §4.4 |
| CN explained by top-k cycles, subwalks satisfied | (a) modified | CoRAL | `aa_eval.py::cn_explained`, `::subwalks_satisfied` | e4c2409 | E §4.5 |
| Strict and best-of modes | (c); best-of after CoRAL | CoRAL `score_reconstruction` | `aa_eval.py::strict_select`, `::best_of_rows` | e4c2409 | E §4.6 |
| `aa_export`, `aa_import`, round trip | (c) | CoRAL / CycleExtractor formats | `aa_export.py`, `aa_import.py` | e4c2409, 2f921a5 | E §5.1–5.3 |
| CoRAL flank fix | (c) | works around CoRAL `get_eulerian_node_constraints` | `coral_flank.py::add_flanks` | bb15c21 | E §5.4 |
| Set metric | (c) | — | `results/bench_nested/diag/grow_rules.py` | de5a805 | E §6.2 |
| bench_nested harness, λ rule, failure stages, `merged` | (c) | — | `results/bench_nested/collect.py::choose`, `::failure_stage` | 949f796, 2fdf73d | bench_nested §2–3 |
| Graph route arms and scoring | (c); CoRAL and CycleExtractor are (a) | CoRAL; CycleExtractor | `results/graph_route/graph_route.py`, `score_h2h.py` | 1c0f3b4, b4ff12d | graph_route §1 |

## Patches (category b)

None of these changes a result.

| Patch | Where | Commit |
|---|---|---|
| gurobipy ≥ 9: `addConstr(lhs, sense, rhs)` → `addLConstr`; module `gurobi` → `gurobipy` | `rck/core/ilp_gurobi.py` | version-fix@0b41976 |
| networkx ≥ 2.4: local `connected_component_subgraphs` helper | `rck/core/graph.py`, `rck/utils/adj/*.py` | version-fix@0b41976 |
| PyVCF → pyvcf3; raw strings for two escape sequences; networkx floor | `setup.py`, `rck/rck_run.py` | version-fix@0b41976 |
| Gurobi thread defaults from Slurm, quiet logger, model disposal | `rck/ecdna/gurobi_env.py` | 41f5804, 8f0f810 |
| Test temporary folders that work on a full `/tmp` and on NFS | `tests/ecdna/helpers.py` | 8f0f810 |
| `evaluation.json` writes `null` instead of `NaN` (strict JSON) | `rck/ecdna/evaluate.py` | 4e53747 |
| WLS licence: `disposeDefaultEnv`, retries, at most 2 concurrent array tasks | `results/bench_nested/run_case.py`, `*.sbatch`, `graph_route.py` | 949f796, 716f4c5, 1c0f3b4 |
| Environments, tools and run outputs moved to `/n/fs/yubi-iw` | sbatch files, scripts, READMEs | RCK 89a5456, AMP eea71aa |
| ecSimulator looked up in the repository's `third_party/` | AMP `simulations/src/writers.py` | AMP 5bdb8e5 |
| CycleExtractor run with the coral-3.1.1 environment's Python | `results/graph_route/envs/` | 2f921a5 |

## Why the new steps exist

The main failures that drove (c) changes, and what each change did. The stage docs give the full
numbers.

| Failure (case) | Cause | Change | Effect |
|---|---|---|---|
| demo7: nested species not recovered (composites 0.43 / 0) | Phase blocks the subgraph could not orient gave segment allele observations with A and B swapped | `segment_blocks.tsv` (S2 §3.5) + `--block-flip unoriented` (S3 §4.6) | Flow objective 39.3 → 9.13; composites 0.994 / 0.9996 |
| demo8: right answer only by a 0.004 tie-break over a merged walk | Copy numbers cannot separate nested species on one haplotype | `--step-penalty` λ in the ranking + ambiguity report (S3 §5.7–5.8) | Same answer decided by a gap of 1.9 at λ = 0.1; demo7 total-only 0.51/0.34 → 0.98/0.97 at λ ≥ 0.05 |
| bench `ctrl_dup_s702`, `nest_cross_3sp_s402`: impossible truth labellings (AA + BB on one junction, crossing fold-backs) | Simulator drew a haplotype for every step of the walk | `hap_labels: inherit` (S1 §S3) | s702 0.744 → 0.997 |
| bench `ctrl_separate_s601/s602`: amplified segments dropped | Candidate regions were ±25 kb around seeds | `--max-extend` region growth (S2 §4.3) | s601 segments 13/16 → 16/16, composite 0.759 → 0.982 |
| bench_nested_v2: `paths.tsv` empty or thin | Paths dropped unless every boundary was a node | `clip_path` (S2 §5.5) | 0–10 → 6–23 paths per case |
| `nest_cross_3sp_s403`: segment copy numbers 2× the truth | Ploidy fit chose class 6 | Explicit `--rck-ploidy 2` for the bench (S2 §5.1) | Scale right; RCK still fails (0.545), so the failure is in the decomposition |
| CoRAL explained 0 copies on our graphs | CoRAL balances an interval's end nodes only through s/t edges | `coral_flank` (E §5.4) | Truth cycles feasible in CoRAL's model 28 → 66 of 67 |
| AA `error_free` passed 1 of 42 amplicons | Imputed segment boundaries count as false segments | `error_free_merged` (E §4.3) | 38 of 42 pass |
