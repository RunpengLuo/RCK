# Plan: multi-clone RCK-ecDNA (branch `yubi_CNP`)

*Status: plan agreed 2026-11-13; not yet implemented.*

## 0. Goal

Make the RCK-ecDNA pipeline as close as possible to the original RCK (Aganezov & Raphael 2020), with only the modifications ecDNA needs.

- The subgraph construction stage writes **one input folder, the annotated DIAG**. It uses RCK's own input formats and holds everything RCK needs.
- RCK-ecDNA runs a version of RCK's `OptModelMultiClone` that is changed as little as possible. From it, it infers:
  - per-clone integer karyotypes, split into a chromosomal layer and an ecDNA layer;
  - ecDNA species topologies;
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
2. **Per-clone integer karyotypes, split into two layers.**
   - For clone k: c_{k,s,h} = x^chr_{k,s,h} + x^ec_{k,s,h}.
   - Balance holds in each layer separately.
   - Telomere excess is allowed only in the chromosome layer, so the ecDNA layer is a circulation.
   - This removes from the new route: the bulk single-sample model, `--bg-max`, the background flank prior (β), `bg_balance_weight` (ε_bal), γ and the depth-bias option.
3. **Chromosomal SVs are inferred per clone, as RCK does.**
   - Novel adjacencies get y_{k,j,φ} in either layer, with RCK's one-labelling-across-clones rule. They are not read from the truth.
   - Restore RCK's 1 − fp usage constraint as the default, and keep `--novel-penalty` as an option.
4. **Per-clone ecDNA abundance is inferred.** It is identifiable because the input is per clone.
   - All clones' ecDNA layers are decomposed together.
   - Species topologies are shared, and each species σ has integer copies ψ_{σ,k} in each clone, with Σ_σ ψ_{σ,k} m_σ(e) = x^ec_k(e) for every clone k and edge e.
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
- **Add flank segments.** For each region, add one background segment on each side, about 1 Mb, carrying its clone copy numbers. These pin the chromosome layer: a chromosome-layer copy must pass through a flank with normal copy number. This replaces `region_flank_cn` and the background prior, in a way RCK itself would handle.
- The existing phases 1–2, seeds, trimming and allele annotation stay unchanged.

### 3.3 Output: the annotated DIAG folder, in RCK's formats

Written next to the existing subgraph outputs, for example `<subgraph>/rck_input/`:

| file | RCK flag | contents |
|---|---|---|
| `rck.scnt.tsv` | `--scnt` | refined subgraph segments plus flanks; `extra` = `cn={clone1:{A,B},...}` in the observed frame; tumor clones only |
| `rck.adj.tsv` | `--adjacencies` | reference and novel adjacencies in RCK's adjacency format; `extra` carries `junction_hap`, `#NOV`, `#REF`, `imputed` |
| `rck.fragments.tsv` | `--fragments` | phase blocks as fragments; junction-vote orientation kept in `extra` for diagnostics only |
| `rck.telomeres.tsv` | `--telomere-positions` | subgraph terminals, i.e. the outer ends of the flank segments (chromosome layer only) |
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

### 4.1 Model: a minimal subclass of `OptModelMultiClone`

New class in `rck/core/`, for example `OptModelMultiCloneEcDNA(OptModelMultiClone)`. It overrides only what the layers require:

- **Variables.** Every segment copy-number variable c_{k,s,h} becomes x^chr_{k,s,h} + x^ec_{k,s,h}. Every adjacency copy-number variable y_{k,j,φ} becomes y^chr + y^ec.
- **Balance (`define_constraints_on_nodes`).** Balance is written separately for each layer:
  - chromosome layer: RCK's rule as is, with excess allowed at telomere positions;
  - ecDNA layer: strict equality everywhere.
- **Labelling and presence.** RCK's rules apply to the sum of the two layers. One labelling per novel adjacency is shared by all clones and both layers, and the 1 − fp constraint is restored.
- **Everything else is RCK's code unchanged:** the objective on the total c = x^chr + x^ec, the deltas, the fragment flips, the scnb bounds, the reference labelling and the molecule groups.

### 4.2 ecDNA additions, each off by default and documented as a deviation

- **Junction-read term** (θ, `--alpha-novel`/`--alpha-ref`). It needs θ from `scale.json` and is applied per clone with u_k.
- **Imputed adjacencies** with `--imputed-penalty`. With them off, a circularising junction the SV caller missed can't be recovered. The benchmark will show what that costs.
- **Linear ends** (`--linear-ends`).
- **A tie-break per telomere copy**, only if risk 1 in section 5 shows up.

### 4.3 Decomposition (`rck/ecdna/decompose.py`)

- The one genuinely new piece: decompose all clones' ecDNA layers together. Topologies are shared; ψ_{σ,k} are integers.
- Reuse the existing cycle search, connectivity cuts, long-read constraints, step penalty and ranking where possible.
- `refit.py` mostly drops out for this route, since the weights are integers per clone.

### 4.4 Command line, outputs, evaluation, docs

- `rck-ecdna --model {bulk,multiclone}` (default `bulk`), with `--rck-input <dir>` pointing at the DIAG folder. `OptModelCircular` and the bulk route are not changed.
- Outputs:
  - RCK's native per-clone `scnt` and `acnt`, with a `layer` field (`chr`/`ec`);
  - `decompositions.tsv` with `psi_<clone>`, `proportion_<clone>` and `proportion_bulk`;
  - `summary.json` with `model: multiclone`.
- `rck/ecdna/evaluate.py` and `docs/pipeline/evaluation.md`: score ψ_{σ,k} and the proportions against the `ecdna_psi` rows of `truth/proportions.tsv`.
- `docs/ecDNA.md`: a new section describing the multi-clone model, with a table of every deviation from `OptModelMultiClone`.
- Tests in `tests/`: build the model on a tiny synthetic DIAG of 1–2 clones and a single-cycle ecDNA. Check that the ecDNA layer is a circulation, that the flips match RCK's behaviour, and that the decomposition recovers ψ.

### 4.5 Benchmark

- `results/bench_nested/run_case.py`: add `--model multiclone` and a `--subgraphs` root for the `_cnp` subgraphs, with results under `results/bench_nested_cnp/`. Same 40 cases.
- demo3 and demo5–8 have `truth/rck/rck.scnt.tsv`.
- demo4 has 3 tumor clones and 2 tumor samples, and RCK takes one per-clone copy-number set. Check that the samples agree before adding demo4.

---

## 5. Risks

1. **Splitting copies between the two layers can be ambiguous.** Within a clone, extra copies could run linearly through the chromosome layer from one telomere to another instead of closing in the ecDNA layer. Original RCK has no reason to prefer either, because it never separates the two.
   - The flank segments (section 3.2) are the faithful fix.
   - If they're not enough, the minimal tie-break is a small cost per telomere copy (section 4.2), documented as a deviation.
2. **MILP size.** The model grows by about (number of tumor clones) × 2 layers. That should be fine for the bench (2 tumor clones). Check solve times on small cases first.
3. **Optimism.** Exact truth copy numbers make this route look better than it would on real HATCHet output (decision 5).

## 6. Order of work

1. amplicon-allele-analysis: `clone_profiles.py`, segment refinement, flanks, DIAG writer, tests.
2. Validate the DIAG with the **unmodified** `rck` on demo5 (small).
3. RCK: `OptModelMultiCloneEcDNA`, the CLI route, outputs, tests on tiny synthetic graphs.
4. RCK: joint decomposition, evaluation.
5. Small end-to-end runs (demo5, demo7, one or two bench cases).
6. Full benchmark: **only with the owner's explicit approval** (section 7).

## 7. Compute rule

Testing on small cases (unit tests, tiny synthetic graphs, single demo cases with short time limits) is fine. **Do not run large Gurobi jobs without Yubi's explicit approval.** That includes the 40-case bench, sweeps, many cases in parallel, long `--time-limit` values and SLURM submissions.
