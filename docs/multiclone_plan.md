# Plan: multi-clone RCK-ecDNA (branch `yubi_CNP`)

*Status: plan agreed 2026-11-13. Revised 2026-11-13 after step 3 (risk 1): the karyotype is no
longer split into a chromosome and an ecDNA layer. RCK's `OptModelMultiClone` runs as is, and the
decomposition separates circles from linear chromosomes, heaviest walk first, as in
AmpliconArchitect and CoRAL (decision 2, sections 4.1 and 4.3, risk 1). Steps 1–4 are done; status
in [`multiclone_NEXT_STEPS.md`](multiclone_NEXT_STEPS.md).*

## 0. Goal

Make the RCK-ecDNA pipeline as close as possible to the original RCK (Aganezov & Raphael 2020), with only the modifications ecDNA needs.

- The subgraph construction stage writes **one input folder, the annotated DIAG**. It uses RCK's own input formats and holds everything RCK needs.
- RCK-ecDNA runs RCK's `OptModelMultiClone`, changed as little as possible (with its additions off,
  only an infeasible edge case differs). From it, it infers:
  - per-clone integer karyotypes, exactly as RCK does;
  - a decomposition of every clone's karyotype into shared walks, heaviest first: closed walks
    (cycles) are ecDNA species, and walks between telomere positions (paths) are chromosomal;
  - each species' copies per clone and its proportions.

On the simulator, the tumor clones' copy-number profiles and proportions come from `truth/`. They are treated **as if HATCHet had inferred them and handed them to us**. This is the only deliberate exception to the "read only `obs/`" rule, and it must be documented as such.

### Repositories and branches

| repo | new branch | created from |
|---|---|---|
| https://github.com/RunpengLuo/amplicon-allele-analysis | `yubi_CNP` | `yubi_dev` @ `5447a99` |
| https://github.com/RunpengLuo/RCK | `yubi_CNP` | `yubi_dev` @ `1860fbe` |

Work only on `yubi_CNP` in both repositories. `yubi_dev`, and the existing bulk route (`OptModelCircular`), must keep working unchanged.

---

## 1. What the original RCK does

The source is `rck/core/ilp_gurobi.py::OptModelMultiClone` and `rck/rck_run.py` in the RCK repo.

- **Data term.** The only data term is per-clone allele-specific segment copy numbers (`--scnt`). The objective is Σ_clone Σ_s L_s(δ_A + δ_B), with δ measured against the input (`define_constraints_on_deltas`, `define_objective`).
  - RCK never uses bulk depth, read counts or clone proportions; uᵢ does not appear in the model.
  - Unmixing the bulk signal into clones is HATCHet's job, upstream of RCK.
- **Inferred karyotypes.** RCK infers integer haplotype-specific segment and adjacency copy numbers for every clone.
  - Balance holds at every node, and copy-number excess is allowed only at declared telomere positions (`define_constraints_on_nodes`).
  - A novel adjacency has one labelling shared by all clones (`define_constraints_for_novel_labeling`).
  - At least 1 − fp of the novel adjacencies must be used (`define_constraints_nov_adjacency_overall_presence`, `nov-adj-fp`).
- **Machinery the ecDNA pipeline reinvented.**
  - The allele flip per **fragment** (`FRAGMENT_ALLELE`, shared across clones) is what the phase-block flip (`--block-flip`) does now.
  - **Molecule** adjacency groups (`--adjacency-groups`, `AdjacencyGroupType.MOLECULE`) are what long-read paths (`paths.tsv`) are.

---

## 2. Decisions

1. **Clone copy-number source: `truth/rck/rck.scnt.tsv`**, the file with the ecDNA copies folded in. Do *not* use `truth/clones.tsv`.
   - `rck.scnt.tsv` is what HATCHet would report and what RCK consumes; see `simulations/src/writers.py::_clone_cn_with_ecdna` in amplicon-allele-analysis.
   - `clones.tsv` holds chromosomal copies only. It would be right only for a bulk-background model, which this plan does not use.
   - **Drop the normal clone.** RCK models tumor clones only.
   - **Don't leak the true phase.** `rck.scnt.tsv` is in the true haplotype frame, but HATCHet's output is allele-specific, in an arbitrary frame. Rewrite each clone's (A, B) in the **observed phase-block frame** of `obs/`, and define RCK's fragments as those phase blocks. RCK's own fragment flip then has to find the phase, as in the paper.
   - **Clone proportions:** the `U` column (cell fraction) of `truth/proportions.tsv`, `KIND = clone`, tumor sample, tumor clones only. Do not use `V`, which is each clone's share of reads after length correction.
2. **Per-clone integer karyotypes, as RCK infers them; circles vs chromosomes from the decomposition.**
   - *Revised 2026-11-13.* The original decision split each clone's karyotype into a chromosome
     layer and an ecDNA layer (c = x^chr + x^ec, balance per layer, telomeres in the chromosome layer
     only). Step 3 showed that the split is not identified (risk 1, section 5). A circulation balances
     in either layer, and RCK's objective sees only the totals. On demo5, all 20 ecDNA copies came
     out in the chromosome layer, and at the same objective the ecDNA layer could hold 0 to 20.
   - Original RCK has no mechanism that tells a circle from chromosomal copies. It infers copy
     numbers, and its only rule about ends, telomere positions, does not apply to a cycle.
   - AmpliconArchitect and CoRAL do not split the copy numbers. They decompose the whole amplicon
     into walks. A walk that closes is a cycle; a linear walk must end at designated source/sink
     points. They extract the heaviest walk first.
   - So the karyotype model is RCK's own, and the decomposition (section 4.3) does the separation:
     - cycles are ecDNA;
     - linear walks may end only at RCK's telomere positions (the outer flank ends; CoRAL's s–t
       paths) and are chromosomal;
     - walks are extracted heaviest first (AA's rule, in CoRAL's greedy form).
   - The new route has none of the bulk model's parts: the single-sample model, `--bg-max`, the
     background flank prior (β), `bg_balance_weight` (ε_bal), γ and the depth-bias option.
3. **Chromosomal SVs are inferred per clone, as RCK does.**
   - Novel adjacencies get y_{k,j,φ} per clone, with RCK's one-labelling-across-clones rule. They are not read from the truth.
   - Restore RCK's 1 − fp usage constraint as the default, and keep `--novel-penalty` as an option.
4. **Per-clone ecDNA abundance is inferred.** It is identifiable because the input is per clone.
   - All clones' karyotypes are decomposed together.
   - Walk topologies are shared, and each walk σ has integer copies ψ_{σ,k} in each clone, with
     Σ_σ ψ_{σ,k} m_σ(e) = c_k(e) for every clone k and edge e (RCK's segment and adjacency copy
     numbers). The ecDNA species are the cycles among them.
   - Proportion per topology is reported two ways:
     - per clone: ψ_{σ,k} / Σ_σ ψ_{σ,k};
     - bulk: Σ_k u_k ψ_{σ,k}, normalised over species. This is the only place the clone proportions are used.
   - The truth to score against is in the `ecdna_psi` rows of `truth/proportions.tsv` (`CLONE`, `PSI`).
5. **No noise first; add perturbation later.**
   - Truth copy numbers are exact integers. Real HATCHet output is rounded, capped (`cnmax`) and averaged over clones in high-copy amplicons.
   - A later option should perturb the per-clone table, by rounding error, a copy-number cap and swapping clone labels, to test robustness.

---

## 3. Subgraph construction (amplicon-allele-analysis `yubi_CNP`)

### 3.1 New module `simulations/clone_profiles.py` (the HATCHet stand-in)

- `load_truth_profiles(run_dir, sample)` reads `truth/rck/rck.scnt.tsv`, where the `extra` column holds `cn={'normal': {'A':..,'B':..}, 'clone1': {...}, ...}`, plus `truth/proportions.tsv` (`U`). It returns a long table with columns `chrom, start, end, clone, u, cn_A, cn_B`, tumor clones only, in the true frame.
- `load_hatchet_profiles(best_seg_ucn, sample)` returns the same table from real HATCHet output. RCK's `rck/ecdna/scale.py::_parse_hatchet` already parses this format and can be reused or copied.
- `to_observed_frame(table, phase_blocks)` re-expresses (A, B) per phase block in the observed `obs/` frame, before any junction-vote orientation. That is the frame RCK's fragment flip has to resolve.
- The docstring must state that the truth reader is a deliberate HATCHet stand-in and the only use of `truth/` in the pipeline.

### 3.2 Changes to the subgraph

- **Refine segments at copy-number changes.** Wherever any clone's copy number changes inside a subgraph segment, split the segment. Each split adds an extremity pair joined only by a reference adjacency. This is required because RCK assumes copy numbers are constant within a segment.
- **Add flank segments.** For each region, add one background segment on each side, about 1 Mb, carrying its clone copy numbers. Their outer ends are the telomere positions, so a linear (chromosomal) walk must pass through a flank with normal copy number. This replaces `region_flank_cn` and the background prior, in a way RCK itself would handle.
- The existing phases 1–2, seeds, trimming and allele annotation stay unchanged.

### 3.3 Output: the annotated DIAG folder, in RCK's formats

Written next to the existing subgraph outputs, for example `<subgraph>/rck_input/`:

| file | RCK flag | contents |
|---|---|---|
| `rck.scnt.tsv` | `--scnt` | refined subgraph segments plus flanks; `extra` = `cn={clone1:{A,B},...}` in the observed frame; tumor clones only |
| `rck.adj.tsv` | `--adjacencies` | reference and novel adjacencies in RCK's adjacency format; `extra` carries `junction_hap`, `#NOV`, `#REF`, `imputed` |
| `rck.fragments.tsv` | `--fragments` | phase blocks as fragments; junction-vote orientation kept in `extra` for diagnostics only |
| `rck.telomeres.tsv` | `--telomere-positions` | subgraph terminals, i.e. the outer ends of the flank segments; the only places a linear walk may end |
| `rck.adjgroups.tsv` | `--adjacency-groups` | molecule groups, one per long-read path from `paths.tsv` with at least 2 junctions |
| `rck.scnb.tsv` | `--scnb` | copy-number bounds per clone, segment and haplotype, from the input, as RCK's `--pre-scnb` would make them |
| `clone_proportions.tsv` | (new) | `clone, u`, tumor clones only; used only for bulk reporting |

- `edges.tsv`, `nodes.tsv`, `CN_TOTAL`/`CN_A`/`CN_B`, `scale.json` (ASCAT γ, θ) and the other files stay. They keep the old route working and serve as diagnostics. The new RCK route does not read them, except θ, if the optional junction-read term is on.
- New options on `subgraph_construct.py`:
  - `--clone-profiles {none,truth,hatchet}`;
  - `--truth-dir` or `--hatchet-seg` for the source;
  - `--clone-sample`;
  - `--flank-bp` (default 1,000,000).
- With `none`, the outputs must be byte-identical to `yubi_dev`.

### 3.4 Checks, benchmark wiring, docs and tests

- **First milestone:** the DIAG folder must load and solve with the **unmodified** original `rck` (`rck/rck_run.py`) on a small case such as demo5. That checks the format before any ecDNA changes.
- `verify_subgraph.py`: report the number of refined segments and flank segments, and check that the per-clone copy numbers in the DIAG match `truth/rck/rck.scnt.tsv` after undoing the frame change.
- `bench_nested.py`: add a variant that builds `subgraphs/bench_nested_cnp/` with `--clone-profiles truth`. The existing subgraphs stay untouched.
- `simulations/subgraph_construct.md`: add a section on clone copy-number annotation and the RCK-format output. Edit the "Only `obs/` is read" statement to name the one deliberate exception.
- Tests in `simulations/tests/`: parsing `rck.scnt.tsv`, the frame conversion, segment refinement, flank segments, the formats of the written files, and a check that the `none` mode is unchanged.

---

## 4. RCK-ecDNA (RCK `yubi_CNP`)

### 4.1 Model: RCK's `OptModelMultiClone`, with a thin subclass for the additions

*Revised 2026-11-13; the two-layer version is described in risk 1 (section 5).*

`rck/core/ilp_gurobi_ecdna.py::OptModelMultiCloneEcDNA(OptModelMultiClone)`:
- **With its options off, it builds exactly RCK's model:** the same variables, constraints and
  objective. A test checks this against `OptModelMultiClone`. On demo5 the route's
  `rck.scnt.tsv` and `rck.acnt.tsv` are byte-identical to `rck_run.py --no-pre`'s.
- **One exception:** RCK's 1 − fp rule is skipped when there are no observed novel adjacencies.
  RCK would write 0 ≥ 1 − fp, which is infeasible (demo5 `amplicon_2` has only imputed ones).
- **No layers.** Which copies form circles is decided by the decomposition (4.3).

### 4.2 ecDNA additions, each off by default and documented as a deviation

- **Junction-read term** (θ, `--mc-junction-weight`). It needs θ from `scale.json` and is applied
  per clone with u_k.
- **Imputed adjacencies** with `--mc-imputed-penalty`. With them off, a circularising junction the
  SV caller missed can't be recovered. The benchmark will show what that costs. They are left out
  of the 1 − fp rule.
- **Linear ends** (`--mc-linear-ends`): extra extremities where copies may end, like RCK's
  telomere positions but at a cost per copy. In the decomposition, these are extra places a linear
  walk may end.
- *Dropped 2026-11-13:* the per-telomere tie-break. It could not fix risk 1, because a cycle uses
  no telomere.

### 4.3 Decomposition, heaviest first (`rck/ecdna/`)

*Revised 2026-11-13.* AmpliconArchitect's rule: extract the heaviest walk, subtract it, repeat.
CoRAL implements it in its greedy mode (`coral/cycle_decomposition.py::maximize_weights_greedy`,
CoRAL 3.1.1): a one-walk MILP that maximises the walk's length-weighted copy number (plus a
long-read term), with the walk's copies capped by the remaining copy number on every edge. It
stops when less than `p_total_weight` (10%) of the length-weighted copy number is left or the
walk weight drops below a resolution. The multi-clone version:

- **Input:** every clone's RCK karyotype, i.e. c_{k,s,h} and N_{k,a,φ}, on the DIAG. The
  remaining copy number R_k(e) starts at those values.
- **One iteration** finds one walk σ on the haplotype-labelled DIAG:
  - either closed (a **cycle**), or with its two ends at telomere positions (a **path**; also at
    `--mc-linear-ends` extremities when that option is on);
  - alternating segment and adjacency edges, with integer multiplicities m_σ(e) up to
    `--max-multiplicity`, as in `decompose.py`;
  - integer copies ψ_{σ,k} ≥ 0 per clone, with ψ_{σ,k} m_σ(e) ≤ R_k(e) for every clone k and edge e.
- **Heaviest** (`--mc-heaviest`; decided in step 4, 2026-11-13):
  - `copies` (default): AA's rule. First maximise the copies, Σ_k ψ_{σ,k}. Then maximise CoRAL's
    length-weighted copy number Σ_k ψ_{σ,k} Σ_s m_σ(s) L_s plus CoRAL's long-read term (RCK
    molecule groups). Remaining ties go to fewer segment traversals.
  - `length`: CoRAL's weight alone.
  - *Why `copies`:*
    - The DIAG's 1 Mb flanks let a long chromosome path that loops through the amplicon outweigh a
      short ecDNA under CoRAL's weight. CoRAL's own graphs have no such flanks.
    - On demo7, `copies` scores composite 0.747 (2/4 ψ exact), against 0.256 (0/4) for `length`.
    - With the junction labelling held to `junction_hap` (risk 5), `copies` is exact (1.0, 4/4)
      and `length` is unchanged (0.256).
  - The walk is then subtracted from every clone's R_k, and the search repeats.
- **Stop:** when R_k = 0 everywhere, or at `--max-species` walks. RCK's karyotype is an exact
  integer balanced flow, so it can be used up exactly. CoRAL's 10% rule is kept as an option for
  noisy input.
- **Calls:**
  - cycles are ecDNA species, with `--min-ecdna-copies` as in the bulk route;
  - paths are chromosomal.
- **What it gives on demo5:** the 20-copy cycle has about 20 × 1.9 Mb of weight, against at most
  1 × ~6 Mb for any path through the flanks. So it is extracted first, and the remaining single
  copy per haplotype becomes telomere-to-telomere paths.
- **Reuse:** `decompose.py`'s walk constraints (alternation, multiplicity, lazy connectivity cuts,
  linear ends through one virtual end node, long-read constraints) with one walk slot and the
  objective above. `refit.py` drops out, since the weights are integers per clone.
- The bulk route's minimum-species decomposition stays as it is.

### 4.4 Command line, outputs, evaluation, docs

- `rck-ecdna --model {bulk,multiclone}` (default `bulk`), with `--rck-input <dir>` pointing at the DIAG folder. `OptModelCircular` and the bulk route are not changed.
- Outputs:
  - RCK's native per-clone `scnt` and `acnt`, exactly as `rck_run.py` writes them;
  - `decompositions.tsv`, one row per walk in extraction order, with `kind` (`cycle`/`path`),
    `psi_<clone>`, `proportion_<clone>` and `proportion_bulk` (proportions over cycles only);
  - `summary.json` with `model: multiclone`.
- `rck/ecdna/evaluate.py` and `docs/pipeline/evaluation.md`: score ψ_{σ,k} and the proportions against the `ecdna_psi` rows of `truth/proportions.tsv`.
- `docs/ecDNA.md`: a new section describing the multi-clone model, with a table of every deviation from `OptModelMultiClone`.
- Tests in `tests/`: tiny synthetic DIAGs of 1–2 clones and a single-cycle ecDNA. Check that the default model is RCK's, that the flips match RCK's behaviour, and that the decomposition recovers ψ (the cycle with e copies per clone, and one telomere-to-telomere path per haplotype).

### 4.5 Benchmark

- `results/bench_nested/run_case.py`: add `--model multiclone` and a `--subgraphs` root for the `_cnp` subgraphs, with results under `results/bench_nested_cnp/`. Same 40 cases.
- demo3 and demo5–8 have `truth/rck/rck.scnt.tsv`.
- demo4 has 3 tumor clones and 2 tumor samples, and RCK takes one per-clone copy-number set. Check that the samples agree before adding demo4.

---

## 5. Risks

1. **Splitting copies between a chromosome and an ecDNA layer is not identified.** *Confirmed in step 3; resolved by the 2026-11-13 revision.*
   - The flanks did pin linear telomere-to-telomere copies to the flank copy number.
   - But a circulation balances in either layer, so ecDNA copies could sit in the chromosome layer as a closed cycle at the same objective. On demo5 that is what the solver returned.
   - A per-telomere cost could not help.
   - Resolution: no layers in the model; circles vs chromosomes come from the heaviest-first decomposition (4.3), whose linear walks must end at telomeres.
2. **MILP size.** RCK's model grows with the number of tumor clones; the heaviest-first decomposition solves one small one-walk MILP per extracted walk, with per-clone copies. That should be fine for the bench (2 tumor clones). Check solve times on small cases first.
3. **Optimism.** Exact truth copy numbers make this route look better than it would on real HATCHet output (decision 5).
4. **Heaviest-first is greedy.** It need not give the fewest walks.
   - AA and CoRAL accept this; CoRAL post-processes greedy solutions with its full model.
   - Nested species on one haplotype are the likely failure: a merged walk could outweigh the inner species. The bulk route handles this with a minimum-species decomposition and the step penalty.
   - Measure it on small cases in step 5 (demo7, a nested bench case) before the benchmark.
   - If needed, a CoRAL-style post-processing pass is the documented fallback.
   - *Step 4:* on demo7 the nested species are recovered exactly once the labelling is right (risk 5). The error seen comes from the labelling, not from the greedy.
5. **RCK's junction labelling is free.** *Found in step 4.*
   - RCK does not read `junction_hap`. A novel adjacency's labelling (AA/AB/BA/BB) together with the fragment flips can change at the same objective.
   - On demo7, RCK labelled the inner species' junctions sim00005–8 AB/BA. The truth and `junction_hap` both say AA.
   - The karyotype then has no single-haplotype inner cycle, and the decomposition merges it into the outer species.
   - Holding labellings to `junction_hap` (A → AA, B → BB; `mixed` free) recovers demo7 exactly. It would be a new deviation from RCK (hard, or as a soft cost), and needs Yubi's decision.

## 6. Order of work

1. amplicon-allele-analysis: `clone_profiles.py`, segment refinement, flanks, DIAG writer, tests.
2. Validate the DIAG with the **unmodified** `rck` on demo5 (small).
3. RCK: `OptModelMultiCloneEcDNA` (RCK's model plus off-by-default additions), the CLI route, outputs, tests on tiny synthetic graphs.
4. RCK: joint heaviest-first decomposition of the clones' karyotypes (cycles and telomere-to-telomere paths), evaluation.
5. Small end-to-end runs (demo5, demo7, one or two bench cases).
6. Full benchmark: **only with the owner's explicit approval** (section 7).

## 7. Compute rule

Testing on small cases (unit tests, tiny synthetic graphs, single demo cases with short time limits) is fine. **Do not run large Gurobi jobs without Yubi's explicit approval.** That includes the 40-case bench, sweeps, many cases in parallel, long `--time-limit` values and SLURM submissions.
