# Multi-clone RCK-ecDNA: status and next steps (branch `yubi_CNP_ionic`)

*Written 2026-11-13; updated the same day after step 3 and after the plan revision (no layers;
heaviest-first decomposition). The plan is [`multiclone_plan.md`](multiclone_plan.md); this file says where
it stands and what comes next. Changes and their reasons are in the `[yubi_CNP_ionic]` sections of
the amplicon-allele-analysis `CHANGELOG.md` (steps 1–2) and of RCK's `CHANGELOG.md` (step 3).*

## Where things are

| repo | branch | ionic checkout | from |
|---|---|---|---|
| RCK | `yubi_CNP_ionic` | `/n/fs/ugrad/ug28/ym3292/RCK-cnp` (git worktree of `~/RCK`) | `origin/yubi_CNP` @ `a58f167` |
| amplicon-allele-analysis | `yubi_CNP_ionic` | `/n/fs/ugrad/ug28/ym3292/amplicon-allele-analysis-cnp` (git worktree of `~/amplicon-allele-analysis`) | `yubi_dev` @ `5447a99` |

- The GitHub `yubi_CNP` branches are not used: this host has no GitHub credentials.
  - RCK's could be fetched, since the repo is public.
  - amplicon-allele-analysis's could not be fetched.
  - Neither repo can be pushed from this host.
- Commits on `yubi_CNP_ionic`, none pushed: amplicon-allele-analysis up to `2811efe` (steps 1–2)
  and `36d7ff8` (docs); RCK `60ca24d` (`--scnb` fix), `18dd969` (lazy PyVCF), `8e6fb74` (model,
  first version with layers), `d3b0448` (route, outputs, tests), `a496a12` (docs), `ea99498`
  (layers removed: the model is RCK's own) and the docs commit after it.
- `~/RCK` and `~/amplicon-allele-analysis` stay on `yubi_dev`. In `~/RCK`, `yubi_dev` is 3 commits
  ahead of GitHub (23b5ecf, 97db161, cfbc8ca), none of them pushed. They are not on
  `yubi_CNP_ionic`.

| plan step (section 6) | status |
|---|---|
| 1. amplicon-allele-analysis: clone profiles, refinement, flanks, DIAG writer, tests | **done** (not pushed) |
| 2. the DIAG with the unmodified `rck` on demo5 | **done**: optimal, objective 0, all post-checks pass |
| 3. RCK: `OptModelMultiCloneEcDNA`, CLI route, outputs, tests | **done**, simplified after risk 1: RCK's own model plus off-by-default additions (section 2) |
| 4. RCK: heaviest-first decomposition of the clones' karyotypes (cycles = ecDNA, telomere-to-telomere paths = chromosomes), evaluation against `ecdna_psi` | **next** (section 3) |
| 5. small end-to-end runs (demo5, demo7, one or two bench cases) | open |
| 6. full benchmark | only with Yubi's explicit approval |

## 0. Decisions taken in step 1 (accepted by Yubi, 2026-11-13)

These differ from, or go beyond, the plan:
1. **Imputed adjacencies are in `rck.adj.imputed.tsv`, not in `rck.adj.tsv`.** Plain RCK would
   treat them as SV calls and its 1 − fp usage rule would force them in. The multi-clone route
   should read the file only with `--imputed-penalty` (plan 4.2).
2. **The observed phase frame is decided per block from the obs SNP counts** (log-likelihood of
   the expected BAF against its mirror), not from `truth/phasing.tsv`.
3. **Copy-number bounds are [0, max(10, 2 × max input CN)]**, because RCK's default upper bound of
   10 would cap ecDNA.
4. **Segments are also cut at phase-block ends**, since RCK needs fragments aligned with segments.

## 1. RCK issues found in step 2 (fixed in step 3)

- **`rck_run.py:224` cannot read `--scnb`.** `scnb[clone_id].update(scnb_from_file[clone_id])`
  passes a `SegmentCopyNumberBoundaries`, but `.update()` iterates its argument as a dict. The
  DIAG works around it by putting the bounds in `rck.scnt.tsv` (`cnb=`).
  - **Fixed in `60ca24d`** with `scnb[clone_id].update(scnb_from_file[clone_id]._records)`. The DIAG
    still writes the bounds in both places, which is harmless.
- **`rck_run.py` imports PyVCF (`vcf`)** through `rck/utils/adj/convert.py`, and the gurobi env
  does not have it. Step 2 used an empty `vcf` module on `PYTHONPATH`
  (`<scratchpad>/stubs/vcf/__init__.py`).
  - **Fixed in `18dd969`:** the import is lazy in `convert.py`. No stub is needed now.
- **Step 2 command,** for reference (D = `<subgraph>/rck_input`):

  ```bash
  PYTHONPATH=$HOME/RCK-cnp python $HOME/RCK-cnp/rck/rck_run.py --scnt $D/rck.scnt.tsv \
      --adjacencies $D/rck.adj.tsv --fragments $D/rck.fragments.tsv \
      --telomere-positions $D/rck.telomeres.tsv --adjacency-groups $D/rck.adjgroups.tsv \
      --no-pre --workdir <out> --run-g-time-limit 120 --run-g-threads 2 --post-check-all
  ```

  `--no-pre` is required: RCK's preprocessing would merge same-copy-number segments and replace
  `--fragments`.

- **Found in step 3, not fixed in plain RCK:** with no novel adjacencies, RCK's 1 − fp rule is
  0 ≥ 1 − fp, which is infeasible (demo5 `amplicon_2` has 0 observed and 18 imputed). The subclass
  skips the rule in that case.

## 2. Step 3: done (plan 4.1, 4.2, 4.4)

What was built is described in `docs/ecDNA.md` §11, with the deviation table in §11.1 and the
reason for having no layers in §11.2, and in RCK's `CHANGELOG.md`.
- **`rck/core/ilp_gurobi_ecdna.py::OptModelMultiCloneEcDNA`:** with its options off, it is RCK's
  `OptModelMultiClone`: the same variables, constraints and objective (tested).
  - The only exception: RCK's 1 − fp rule is skipped when there are no observed novel adjacencies,
    where RCK is infeasible.
  - The off-by-default additions are `--mc-imputed-penalty`, `--mc-novel-penalty`,
    `--mc-junction-weight` and `--mc-linear-ends`.
- **`rck/ecdna/multiclone.py` and `rck-ecdna --model multiclone --rck-input <dir>`:** the DIAG is
  read as `rck_run --no-pre` reads it. The outputs are `rck.scnt.tsv` / `rck.acnt.tsv` as `rck_run`
  writes them, plus `summary.json`.
- **`tests/ecdna/test_multiclone.py`:** 11 tests, a few seconds in total.
- **demo5** (`--time-limit 60`, 2 threads, 0.1 s): `rck.scnt.tsv` and `rck.acnt.tsv` are
  byte-identical to the step-2 plain-RCK run. With `--mc-imputed-penalty 2 --mc-junction-weight 0.5`
  (first version), the imputed edges were not used and the junction residual was 0.385, matching a
  hand calculation.
- **Bulk route** on demo5: identical to before, except log timestamps and `runtime_s`.
- **Risk 1, and why there are no layers.**
  - The first version split each clone into a chromosome layer and an ecDNA layer. The split was not
    identified: a circulation balances in either layer, and on demo5 all 20 copies came out in the
    chromosome layer.
  - Original RCK has no way to separate them. Following AA and CoRAL, the separation moved to the
    decomposition (Yubi's decision, 2026-11-13).
  - The demo5 DIAG for checks is rebuilt with `subgraph_construct.py --sim-dir
    runs/demo5_cycle1_ecdna --amplicons 0 --observed-seeds-only --strong-sv-support 20
    --jump-sv-support 20 --rck-junction-reads-per-copy 6.5 --no-plot --clone-profiles truth` (no
    Gurobi, about 2 min). The default arguments give a different subgraph.

## 3. Step 4: heaviest-first decomposition and evaluation (plan 4.3, 4.4)

Yubi's decision (2026-11-13): use AmpliconArchitect's heaviest-cycle-first rule.
- **Source of the rule.** CoRAL 3.1.1 implements it as its greedy mode,
  `coral/cycle_decomposition.py::maximize_weights_greedy` with `models/concrete.py::get_greedy_objective`
  (installed in `/n/fs/yubi-iw/conda/envs/coral-3.1.1`). AA's own source is not on ionic.
- **What CoRAL's greedy mode does.** It solves a one-walk MILP that maximises the walk's
  length-weighted copy number, Σ_seq x_e · w · len_e, plus a bonus for unsatisfied long-read
  subpaths. The walk's copies are capped by the remaining copy number on every edge. It subtracts
  the walk and repeats, until less than 10% (`p_total_weight`) of the length-weighted copy number
  is left or the weight drops below `resolution`. CoRAL can then post-process the greedy solution
  with its full model.

**To build** (`rck/ecdna/`, for example a new `greedy.py` driven from `multiclone.py`):
1. **Input:** every clone's RCK solution, c_{k,s,h} and N_{k,a,φ} on the DIAG, as the remaining
   copy number R_k(e).
2. **One iteration:** find one walk σ on the haplotype-labelled DIAG, with integer copies
   ψ_{σ,k} ≥ 0 per clone and ψ_{σ,k} m_σ(e) ≤ R_k(e).
   - The walk is either closed (a cycle) or has both ends at RCK telomere positions (a path; also
     at `--mc-linear-ends` extremities when that option is on).
   - Its multiplicities m_σ(e) are integers up to `--max-multiplicity`.
   - Objective: maximise Σ_k ψ_{σ,k} Σ_s m_σ(s) L_s, plus CoRAL's long-read term. Break ties
     towards fewer segment traversals, using the step penalty as a tie-break.
3. **Subtract and repeat.** Stop when R_k = 0 for every clone, or at `--max-species` walks; RCK's
   karyotype is exact, so it can be used up exactly. Keep CoRAL's 10% stop as an option.
4. **Reuse `decompose.py`:** alternation, multiplicities, lazy connectivity cuts, linear ends
   through one virtual end node, and long-read constraints. Use one walk slot and the objective
   above. Check whether the ψ × m products need the exact linearisation already in `decompose.py`.
   `refit.py` is not used.
5. **Outputs:** `decompositions.tsv`, one row per walk in extraction order:
   - `kind` (`cycle`/`path`) and the walk itself;
   - `psi_<clone>`;
   - `proportion_<clone>` = ψ_{σ,k} / Σ over cycles of ψ_{·,k};
   - `proportion_bulk` = Σ_k u_k ψ_{σ,k}, normalised over cycles. This is the only use of
     `clone_proportions.tsv` outside the optional junction term.

   Cycles with at least `--min-ecdna-copies` are ecDNA calls.
6. **`rck/ecdna/evaluate.py` and `docs/pipeline/evaluation.md`:** score the cycles' ψ and proportions
   against the `ecdna_psi` rows of `truth/proportions.tsv` (`CLONE`, `PSI`).
7. **Tests** (tiny synthetic DIAGs, `tests/ecdna/test_multiclone.py` has a writer):
   - the single-cycle DIAG gives a cycle with ψ = e per clone, then one telomere-to-telomere path
     per haplotype;
   - two clones and two species with different ψ per clone recover ψ exactly;
   - a tie between "a path looping twice + a cycle of e − 1" and "a single-pass path + a cycle of
     e" is broken towards the latter.
8. **demo5** (single run, short time limit): the expected first walk is the 20-copy cycle on the
   amplified haplotype; it weighs about 20 × 1.9 Mb, against at most about 1 × 6 Mb for any path.
9. **Docs:** `docs/ecDNA.md` §11 (the decomposition, and its deviations from CoRAL's greedy model)
   and the CHANGELOG.

**Known risk (plan risk 4):** greedy is not minimal. A merged walk over nested species on one
haplotype could outweigh the inner species. Check it in step 5 on demo7 and a nested bench case;
CoRAL-style post-processing is the fallback.

**Stop and report after step 4.**

## 4. Step 5: small end-to-end runs

- Build the DIAG for demo7 and one or two bench cases (no Gurobi), then run
  `rck-ecdna --model multiclone` with a short time limit:
  - demo5 (done in step 2, plain RCK), demo7, one nested bench case (for example s103, which
    the bulk route solves) and one the bulk route fails (for example s402).
  - The demo5 subgraph command is in amplicon-allele-analysis `subgraph_construct.md`
    section 5.6.
- demo3 and demo5–8 have `truth/rck/rck.scnt.tsv`. For demo4 (3 tumor clones, 2 tumor samples),
  first check that the samples' clone copy numbers agree (plan 4.5).

## 5. Step 6: full benchmark (only with explicit approval)

- `amplicon-allele-analysis/simulations/bench_nested.py --cnp rebuild` builds
  `subgraphs/bench_nested_cnp/` (no Gurobi, about 35 s per case).
  - It needs the bench simulation runs, which are untracked and live in
    `~/amplicon-allele-analysis/simulations/runs/bench_nested`. Run it with `BENCH_NESTED_ROOT`
    or from that checkout.
- `results/bench_nested/run_case.py`: add `--model multiclone` and a `--subgraphs` root for the
  `_cnp` subgraphs; results go under `results/bench_nested_cnp/`. Same 40 cases. This is a Slurm
  array and needs Yubi's approval. Gurobi WLS has 2 sessions, so use `%1` or `%2`.

## 6. Later

- **Perturbation option** (plan decision 5): round, cap (`cnmax`) and swap clone labels in the
  per-clone table to test robustness, since exact truth copy numbers are optimistic.
- **`--clone-profiles hatchet`** on real data (COLO320) once a `best.seg.ucn` exists.
