# Multi-clone RCK-ecDNA: status and next steps (branch `yubi_CNP_ionic`)

*Written 2026-11-13; updated the same day after step 3, the plan revision (no layers;
heaviest-first decomposition), step 4, and the junction-read labelling groups. The plan is [`multiclone_plan.md`](multiclone_plan.md); this file says where
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
  (layers removed: the model is RCK's own), `690d089` (docs), `fc22a34` (step 4: heaviest-first
  decomposition, ψ evaluation), `8faaace` (docs), `140e938` (within-block phasing score) and the
  docs commit after it; amplicon-allele-analysis `8c89156` (`--diag-labeling-groups`) and
  `c6f24cc` (docs), `54226ed` (per-end labels) and `380b3fc` (docs).
- `~/RCK` and `~/amplicon-allele-analysis` stay on `yubi_dev`. In `~/RCK`, `yubi_dev` is 3 commits
  ahead of GitHub (23b5ecf, 97db161, cfbc8ca), none of them pushed. They are not on
  `yubi_CNP_ionic`.

| plan step (section 6) | status |
|---|---|
| 1. amplicon-allele-analysis: clone profiles, refinement, flanks, DIAG writer, tests | **done** (not pushed) |
| 2. the DIAG with the unmodified `rck` on demo5 | **done**: optimal, objective 0, all post-checks pass |
| 3. RCK: `OptModelMultiCloneEcDNA`, CLI route, outputs, tests | **done**, simplified after risk 1: RCK's own model plus off-by-default additions (section 2) |
| 4. RCK: heaviest-first decomposition of the clones' karyotypes (cycles = ecDNA, telomere-to-telomere paths = chromosomes), evaluation against `ecdna_psi` | **done** (section 3); `junction_hap` passed as RCK labelling groups within phase blocks (section 3a) |
| 5. small end-to-end runs (demo5, demo7, one or two bench cases) | demo5 and demo7 done in step 4; bench cases open |
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

## 3. Step 4: done (plan 4.3, 4.4)

Described in `docs/ecDNA.md` §11.3 (with a table of deviations from CoRAL's greedy model) and in
`docs/pipeline/evaluation.md` §3.9.
- **`rck/ecdna/greedy.py`:**
  - Builds the DIAG as an `ECGraph` and every clone's RCK karyotype as phased edge copy numbers,
    with telomere excess as linear ends.
  - Each round solves one one-walk MILP, reusing `decompose.py`'s alternation, multiplicities,
    lazy connectivity cuts and virtual end node, with integer ψ per clone capped by what is left.
  - The walk is subtracted and the search repeats until nothing is left.
  - Cycles are ecDNA; paths, which end only at telomeres, are chromosomes.
- **Heaviest (`--mc-heaviest`):**
  - `copies` (default): AA's rule. The most copies first, then CoRAL's length-weighted copy number
    (plus the long-read bonus over RCK molecule groups), then fewer traversals.
  - `length`: CoRAL 3.1.1's greedy weight alone.
- **Outputs:** `decompositions.tsv` with `kind`, `order`, `psi_<clone>`, `proportion_<clone>` and
  `proportion_bulk`; paths have status `chromosomal`, which `evaluate.read_inferred` skips.
  `evaluate.psi_scores` writes `evaluation_psi.tsv`.
- **Tests:** `tests/ecdna/test_greedy.py`, 4 tests. The full suite is 267 tests; only the 7 old
  failures remain, which were already there at `32f2244`.
- **Runs** (60 s limit, 2 threads, about 1 s each; DIAGs rebuilt in the scratchpad):

| case | rule | walks | mean composite | ψ exact |
|---|---|---|---|---|
| demo5 | copies | cycle ψ 20, then 4 paths | 1.000 | 1/1 |
| demo5 | length | the same | 1.000 | 1/1 |
| demo7 (truth ψ 26/13, 14/7) | copies | cycle 14/7 (inner merged with part of the outer), cycle 12/6, 2 paths | 0.747 | 2/4 |
| demo7 | length | cycle 13/6, a 175-edge looping path, a path, cycle 0/1 | 0.256 | 0/4 |
| demo7, labelling held to `junction_hap` (scratch probe, not in the repo) | copies | cycle 26/13, cycle 14/7, 2 paths | 1.000 | 4/4 |
| demo7, the same | length | as above | 0.256 | 0/4 |

- **Conclusion on the rule.** CoRAL's greedy weight does not give better results here; AA's copy
  rule does. CoRAL's weight lets a long low-copy chromosome path through the DIAG's 1 Mb flanks
  outweigh a short ecDNA, and that path loops through the amplicon.
- **demo7 DIAG** (no Gurobi, 23 s): `subgraph_construct.py --sim-dir runs/demo7_nestedcycles_ecdna
  --union --observed-seeds-only --strong-sv-support 20 --jump-sv-support 20
  --rck-junction-reads-per-copy 6.5 --no-plot --clone-profiles truth`. It has 49 segments,
  10 novel and 19 imputed adjacencies, 3 molecule groups and 2 clones.

## 3a. Junction read labels as RCK labelling groups (done 2026-11-13, Yubi's decision)

- **The problem.** RCK fixes a junction's labelling only where allele-specific copy numbers force
  it. Between phase blocks that only novel adjacencies connect, which is typical for ecDNA, flips
  and labellings are free at the same objective. On demo7, the inner species' junctions came out
  AB/BA (truth AA), and the species merged.
- **The fix, input-side.** amplicon-allele-analysis `subgraph_construct.py --diag-labeling-groups
  block` (`rck_diag.py::labeling_groups`; its `subgraph_construct.md` §5.6):
  - one RCK labelling group (`agt=L`) per phase block and `junction_hap` label;
  - never across blocks;
  - groups that RCK's refinement would merge across blocks or labels are dropped.

  RCK's model is unchanged (`docs/ecDNA.md` §11.4).
- **Evaluation:** `*_blockwise` scores, with one A/B swap per phase block shared by all species
  (`docs/pipeline/evaluation.md` §3.10; `rck-ecdna-eval --phase-blocks`). The multi-clone route
  passes the DIAG's fragments automatically.
- **Runs** (DIAGs rebuilt in the scratchpad, 60 s limit, 2 threads, under 1 s each):

| case | groups | cycles ψ | ψ exact | composite | within blocks |
|---|---|---|---|---|---|
| demo5 | none / 1 | 20 | 1/1 | 1.000 | 1.000 |
| demo7 | none | 14/7 (merged), 12/6 | 2/4 | 0.747 | 0.747 |
| demo7 | 1 (chr4:171, A, 8 ends) | 26/13, 14/7 | 4/4 | 0.693 | 1.000 |

- **Still open:** phase *between* blocks is undetermined by the data (hence the within-block
  score). Junction ends alone in their block stay free.
- **Per-end labels (done 2026-11-13, replacing one label per junction).** The first version gave
  a junction's single `junction_hap` to both ends, which is wrong for a junction that crosses
  haplotypes.
  - amplicon-allele-analysis `54226ed`/`380b3fc`: the simulator writes
    `obs/adjacency_allele_ends.tsv` (phased reads per end, in the observed frame of the end's
    block), and `--diag-labeling-groups block` builds the groups from per-end labels.
  - **Correction.** The earlier measurement (16 of 88 A/B labels on crossing junctions, 3
    opposite) used subgraphs built from bench runs of 2026-10-25, which predate the simulator's
    2026-10-27 junction-read change (`618697f`). The bench runs in `/n/fs/yubi-iw/runs` are stale
    in this respect: 6 of the 10 checked cases no longer regenerate their `obs/` (the pre-change
    code gives the same as the new code, so this session's change is not the cause). On
    regenerated runs, the old labels were all right, and crossing junctions were `mixed`, i.e.
    left out.
  - **Per-end accuracy:** 234 of 234 labelled ends right on demo5, demo7 and 10 regenerated bench
    cases (58 ends of crossing junctions). Scratch scripts `end_accuracy.py`, `bench_case.sh`.
  - **demo7 with per-end groups** (5 groups): ψ 4/4, composite 1.000 absolute and within blocks.
    demo5 unchanged, 1.000.
  - **Before the benchmark:** regenerate the bench simulations with the current simulator. The
    runs on disk lack `adjacency_allele_ends.tsv`, and 6 of the 10 checked are stale.
- **Next: step 5.** One or two nested bench cases (for example s103, and s402, which the bulk route
  fails), with `--clone-profiles truth --diag-labeling-groups block`, as single runs with short
  time limits. Check that `junction_hap` purity holds on the bench, and compare with and without
  the groups.

**Stop and report after step 4.** (Done 2026-11-13.)

## 4. Step 5: small end-to-end runs

- demo5 and demo7 were run in step 4 (section 3). Still open: build the DIAG for one or two bench
  cases (no Gurobi), then run `rck-ecdna --model multiclone` with a short time limit, for example
  s103 (which the bulk route solves) and s402 (which it fails).
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
