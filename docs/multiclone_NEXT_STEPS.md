# Multi-clone RCK-ecDNA: status and next steps (branch `yubi_CNP_ionic`)

*Written 2026-11-13. The plan is [`multiclone_plan.md`](multiclone_plan.md); this file says where
it stands and what comes next. Changes and their reasons are in the amplicon-allele-analysis
`CHANGELOG.md`, section `[yubi_CNP_ionic]`.*

## Where things are

| repo | branch | ionic checkout | from |
|---|---|---|---|
| RCK | `yubi_CNP_ionic` | `/n/fs/ugrad/ug28/ym3292/RCK-cnp` (git worktree of `~/RCK`) | `origin/yubi_CNP` @ `a58f167` |
| amplicon-allele-analysis | `yubi_CNP_ionic` | `/n/fs/ugrad/ug28/ym3292/amplicon-allele-analysis-cnp` (git worktree of `~/amplicon-allele-analysis`) | `yubi_dev` @ `5447a99` |

- The GitHub `yubi_CNP` branches are not used: this host has no GitHub credentials.
  - RCK's could be fetched, since the repo is public.
  - amplicon-allele-analysis's could not be fetched.
  - Neither repo can be pushed from this host.
- `~/RCK` and `~/amplicon-allele-analysis` stay on `yubi_dev`. In `~/RCK`, `yubi_dev` is 3 commits
  ahead of GitHub (23b5ecf, 97db161, cfbc8ca), none of them pushed. They are not on
  `yubi_CNP_ionic`.

| plan step (section 6) | status |
|---|---|
| 1. amplicon-allele-analysis: clone profiles, refinement, flanks, DIAG writer, tests | **done** (not pushed) |
| 2. the DIAG with the unmodified `rck` on demo5 | **done**: optimal, objective 0, all post-checks pass |
| 3. RCK: `OptModelMultiCloneEcDNA`, CLI route, outputs, tests | next |
| 4. RCK: joint decomposition across clones, evaluation against `ecdna_psi` | open |
| 5. small end-to-end runs (demo5, demo7, one or two bench cases) | open |
| 6. full benchmark | only with Yubi's explicit approval |

## 0. Decisions to confirm before step 3

These were taken in step 1 and differ from, or go beyond, the plan:
1. **Imputed adjacencies are in `rck.adj.imputed.tsv`, not in `rck.adj.tsv`.** Plain RCK would
   treat them as SV calls and its 1 − fp usage rule would force them in. The multi-clone route
   should read the file only with `--imputed-penalty` (plan 4.2).
2. **The observed phase frame is decided per block from the obs SNP counts** (log-likelihood of
   the expected BAF against its mirror), not from `truth/phasing.tsv`.
3. **Copy-number bounds are [0, max(10, 2 × max input CN)]**, because RCK's default upper bound of
   10 would cap ecDNA.
4. **Segments are also cut at phase-block ends**, since RCK needs fragments aligned with segments.

## 1. Known RCK issues found in step 2 (not fixed)

- **`rck_run.py:224` cannot read `--scnb`.** `scnb[clone_id].update(scnb_from_file[clone_id])`
  passes a `SegmentCopyNumberBoundaries`, but `.update()` iterates its argument as a dict. The
  DIAG works around it by putting the bounds in `rck.scnt.tsv` (`cnb=`).
  - Fix candidate: `scnb[clone_id].update(scnb_from_file[clone_id]._records)`.
  - This is a bug fix, not a model change; do it in step 3 as its own commit.
- **`rck_run.py` imports PyVCF (`vcf`)** through `rck/utils/adj/convert.py`, and the gurobi env
  does not have it. Step 2 used an empty `vcf` module on `PYTHONPATH`
  (`<scratchpad>/stubs/vcf/__init__.py`).
  - Fix candidates: install `PyVCF3` in the env (ask first), or make the import lazy in
    `convert.py`.
- **Step 2 command,** for reference (D = `<subgraph>/rck_input`):

  ```bash
  PYTHONPATH=$HOME/RCK-cnp:<stubs> python $HOME/RCK-cnp/rck/rck_run.py --scnt $D/rck.scnt.tsv \
      --adjacencies $D/rck.adj.tsv --fragments $D/rck.fragments.tsv \
      --telomere-positions $D/rck.telomeres.tsv --adjacency-groups $D/rck.adjgroups.tsv \
      --no-pre --workdir <out> --run-g-time-limit 120 --run-g-threads 2 --post-check-all
  ```

  `--no-pre` is required: RCK's preprocessing would merge same-copy-number segments and replace
  `--fragments`.

## 2. Step 3: `OptModelMultiCloneEcDNA` (plan 4.1, 4.2, 4.4)

Code in `rck/core/` (for example `rck/core/ilp_gurobi_ecdna.py`), subclassing
`rck/core/ilp_gurobi.py::OptModelMultiClone`.
1. **Read `OptModelMultiClone` in full first.** Look at `define_variables` (l. 282),
   `define_constraints_on_nodes` (l. 704), `get_lin_expr_for_node_haplotype` (l. 719),
   `define_constraints_nov_adjacency_overall_presence` (l. 549),
   `define_constraints_for_novel_labeling` (l. 557), `define_constraints_on_deltas` (l. 751),
   `define_objective` (l. 782) and the getters.
   - Write down which variables are indexed by (clone, segment, haplotype) and (clone,
     adjacency, phasing). These get a layer index.
2. **Variables.** Split each segment and adjacency copy-number variable into a chromosome-layer
   and an ecDNA-layer variable. Keep the original variable as their sum, with one equality
   constraint per pair, so every untouched RCK method still sees the total.
3. **Balance.** Override `define_constraints_on_nodes`:
   - the chromosome layer gets RCK's rule (excess allowed at telomere positions);
   - the ecDNA layer gets strict equality at every node (a circulation).
4. **Unchanged RCK code:** labelling (one labelling per novel adjacency, shared across clones
   and layers), the 1 − fp presence rule, the deltas and objective on the total, the fragment
   flips, the bounds and the molecule groups. List any method that turns out to need a change in
   the deviation table (plan 4.4).
5. **Off-by-default additions** (plan 4.2), each a flag and a row in the deviation table:
   - the junction-read term (θ from `scale.json`, applied per clone with u_k);
   - `--imputed-penalty` (reads `rck.adj.imputed.tsv`);
   - `--linear-ends`.
   - The telomere tie-break is added only if risk 1 shows up (point 8).
6. **CLI.** `rck-ecdna --model {bulk,multiclone}` (default `bulk`) with `--rck-input <dir>`.
   - The bulk route (`OptModelCircular`, current defaults) must stay byte-identical; check it
     with an existing small test.
   - Reuse the I/O of `rck_run.py` (`read_segments_from_file`, `extract_scnt_from_segments`,
     `extract_scnb_from_segments`, ...) rather than new readers, so the route reads exactly what
     plain RCK reads.
7. **Outputs.**
   - RCK's native per-clone `rck.scnt.tsv` and `rck.acnt.tsv`, each with a `layer` field
     (`chr`/`ec`);
   - `summary.json` with `model: multiclone`.
   - `decompositions.tsv` comes in step 4.
8. **Tests** (`tests/`, tiny synthetic DIAGs written with amplicon-allele-analysis
   `rck_diag.py` conventions; Gurobi, seconds each):
   - 1 clone, one single-cycle ecDNA plus flanks: the ecDNA layer is a circulation, the
     chromosome layer is one copy per haplotype through the flanks, and the objective is 0;
   - 2 clones with different ecDNA copies: per-clone layers;
   - a flipped phase block: the fragment flip restores it, as plain RCK did on demo5;
   - **risk 1:** whether extra copies leak into the chromosome layer as linear telomere-to-
     telomere paths. If they do, add the tie-break (a small cost per telomere copy),
     documented as a deviation, and stop and report first.
9. **Then demo5** (single case, short `--time-limit`): compare with the plain-RCK result of
   step 2 (total copy numbers should agree) and look at the layer split.
10. **Docs:** a new section in `docs/ecDNA.md` on the multi-clone model, with a table of every
    deviation from `OptModelMultiClone`; a CHANGELOG entry.

**Stop and report after step 3.**

## 3. Step 4: joint decomposition and evaluation (plan 4.3, 4.4)

- **`rck/ecdna/decompose.py`:** decompose all clones' ecDNA layers together. Topologies are
  shared, with integer copies ψ_{σ,k} per clone and Σ_σ ψ_{σ,k} m_σ(e) = x^ec_k(e) for every clone
  k and edge e.
  - Reuse `CycleDecomposer`'s cycle search, connectivity cuts, read constraints, step penalty
    and ranking where possible. The new part is per-clone weights on one set of topologies.
  - `refit.py` should drop out, since the weights are integers per clone.
- **Outputs:** `decompositions.tsv` with `psi_<clone>`, `proportion_<clone>` =
  ψ_{σ,k} / Σ_σ ψ_{σ,k}, and `proportion_bulk` = Σ_k u_k ψ_{σ,k}, normalised over species. This is
  the only place `clone_proportions.tsv` is used.
- **`rck/ecdna/evaluate.py` and `docs/pipeline/evaluation.md`:** score ψ and the proportions
  against the `ecdna_psi` rows of `truth/proportions.tsv` (`CLONE`, `PSI`).
- **Tests:** two clones and two species with different ψ per clone; recover ψ exactly.

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
