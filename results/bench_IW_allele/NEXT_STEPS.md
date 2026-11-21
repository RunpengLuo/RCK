# bench_IW_allele, item 1 (allelic siblings): following up after the method arrays

*Written 2026-11-20 when the arrays were submitted. The design and the results so far are in
[`PLAN.md`](PLAN.md) ("Item 1 in detail", "Resume here"). This file is the how-to: check the arrays, fix
failures, score, and what is still to be written.*

## RUNNING: RCK MC vs CoRAL topology figures (Slurm array 32468172, submitted 2026-11-21)

**What it is.** Yubi asked (2026-11-21) for topology plots, Jaccard (species) matrices and coverage figures
of RCK MC and CoRAL side by side, on cases where RCK MC recovers nested ecDNA with both alleles amplified.
`make_compare_plots.py` draws them, one case per array task (`compare_plots.sbatch`, 2 CPUs, 4 GB, 30 min, no
Gurobi). Output: `topology_plots/compare/<k>_<case>/`, one folder per case, plus `topology_plots/compare/README.md`
(the case list and why each was chosen).

| task | folder | set | why chosen (scores from the benchmark tables) |
|---|---|---|---|
| 0 | `1_sib_opp_shared_s18004` | iw_allelic | opposite-allele sibling, equal clone weights; RCK MC both exact, CoRAL one cycle over both (sibling missed) |
| 1 | `2_sib_opp_mirror_s17004` | iw_allelic | opposite-allele sibling, mirrored clones; RCK MC both exact; CoRAL 1 of 2 exact + a false species |
| 2 | `3_sib_mixed_mirror_s19012` | iw_allelic | sibling opposite on some intervals, same on others; RCK MC both exact; CoRAL misses the sibling |
| 3 | `4_sib_opp_shared_s18006` | iw_allelic | opposite-allele sibling, reuse stratum; RCK MC both exact; CoRAL none exact |
| 4 | `5_hard_S4_s16005` | iw_intpsi | 4 nested species, both alleles on 56% of ecDNA; RCK MC 3 of 4 exact; CoRAL one cycle, 3 missed |
| 5 | `6_hard_S4_s16002` | iw_intpsi | 4 nested, both alleles on 39%; RCK MC 3 of 4 exact; CoRAL 1 missed, none exact |
| 6 | `7_hard_S4_s16006` | iw_intpsi | (extra) 4 nested, both alleles on 39%; RCK MC 2 of 4 exact; CoRAL 2 missed |
| 7 | `8_hard_S2_s15003` | iw_intpsi | (extra) 2 nested, both alleles on 62%; RCK MC 1 of 2 exact (J_seg 0.82); CoRAL J_seg 0.21 |

Files in each folder (details in `make_compare_plots.py`'s docstring):
- `topology_compare.png`: **the main figure.** One genome axis. Rows: the truth species, then RCK MC's species,
  then CoRAL's. Haplotype A blue, B orange; CoRAL has no alleles, so grey. RCK MC's alleles are shown after the
  scorer's one swap per phase block. "→ ecdnaN" is the truth species each answer species is paired with.
- `truth_topology.png`, `coral_truth_<species>_cycles.png`: the truth alone.
- `rck_mc_rank1_ec<k>.png`, `rck_mc_recovery_matrix.png`: rck-ecdna-plot pages of the RCK MC run.
- `coral_cycles.png`, `coral_graph.png`: CoRAL's own figures of CoRAL's answer and of the graph it decomposed.
  `coral_rck_mc_cycles.png`: RCK MC's species in CoRAL's style.
- `rck_mc_species_matrix_block.*`, `rck_mc_species_coverage_block.*`: RCK MC within phase blocks (alleles must
  match). `rck_mc_species_*_blind.*`, `coral_species_*_blind.*`: both methods haplotype-blind, for a like-for-like
  view.
- `README.txt`: the case, why it was chosen, RCK MC and CoRAL scores, and for the sibling cases the per-base
  sibling-vs-outer check.
- `species_matrix.log`: what the matrix script printed (pairings, coverage shares).

### When it has finished (from `~/RCK-cnp/results/bench_IW_allele`)

1. **Did every task finish?** Each task log ends with `<case>: done -> <folder>`, and lists any step that failed:
   ```bash
   grep -h "done ->\|FAILED" /n/fs/yubi-iw/runs/bench_IW_allelic/logs/slurm/iwa-compare_32468172_*.out
   sacct -j 32468172 -X --format=JobID,State,Elapsed,MaxRSS     # if needed; wait >= 60 s between sacct calls
   ```
   Expect 8 `done` lines and no `FAILED`. A failed step leaves the other figures in place; the traceback is in
   that task's log.
2. **Are all the files there?** Each folder should have `topology_compare.png`, `truth_topology.png`,
   `coral_cycles.png`, `coral_graph.png`, `coral_rck_mc_cycles.png`, at least one `rck_mc_rank1_ec*.png`, and 6
   species-matrix / coverage PNGs (`rck_mc_*_block`, `rck_mc_*_blind`, `coral_*_blind`):
   ```bash
   for d in topology_plots/compare/*_*/; do echo "$d $(ls $d*.png | wc -l) png"; done
   ```
3. **Redo one case** (task index K from the table): `sbatch --array=K --output=/n/fs/yubi-iw/runs/bench_IW_allelic/logs/slurm/%x_%A_%a.out compare_plots.sbatch`,
   or locally `PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH python make_compare_plots.py --index K`
   (1–2 min). It overwrites that case's folder only.
4. **Confirm the results**, case by case (open `topology_compare.png`, then the matrices):
   - RCK MC's rows should repeat the truth's rows: the same blocks in the same order, the same colours. For the
     sibling cases, the outer and the sibling should be on opposite colours (s19012: opposite on some intervals,
     the same on others), matching `README.txt`'s 0 bp disagreeing.
   - CoRAL's rows should show the failure the table says: one long grey cycle running through both species
     (s18004, s16005), a missing sibling (s19012), or extra / partial cycles.
   - `rck_mc_species_matrix_block.png`: an outlined diagonal of 1.00 / 1.00 for the exact species.
     `coral_species_matrix_blind.png`: "missed" columns or low J_seg. Expected misses: the hard_S4 cases have one
     RCK MC species that is not exact (3 of 4 exact; 2 of 4 on s16006), and s15003 has 1 of 2.
   - Coverage figures: RCK MC's truth bars in one colour each; CoRAL's shared by one cycle across species.
5. **Commit the figures** (about 8 folders of PNG / PDF):
   `git add results/bench_IW_allele/topology_plots/compare && git commit -m "bench_IW_allele: RCK MC vs CoRAL topology figures (array 32468172)"`.
   Push from a machine with GitHub credentials (ionic has none).

**State 2026-11-21: §2–§5 done; item 1 finished.** All 48 runs ok. CoRAL and CE on sib_opp_mirror_s17006 were
rerun (jobs 32468001 / 32468002) on an **input patch**: a copy of rck_bulk's graph with its 4 zero-read
discordant edges set to 1 read, because CoRAL and CE cannot take in 0 reads (`allelic_bench.py`,
`patched_graph`; PLAN.md "Item 1 arrays: status and first read"). §5 compared rck_mc, CoRAL and CE only
(Yubi, 2026-11-21). Results in PLAN.md "Item 1 results"; tables in `item1/`; figures
`figures/item1_species_*`. To redo §5 from the scored runs (local, seconds, no Gurobi):

```bash
python allele_bench.py --set iw_allelic run --jobs 16 && python allele_bench.py --set iw_allelic collect
python rck_style.py --set iw_allelic run --jobs 16 && python rck_style.py --set iw_allelic collect
python item1_results.py && python make_item1_figures.py
```

Use the within-block pairing (`item1_results.py`), not the runs' own `evaluation_psi.tsv` (absolute frame; it
mispairs s17008 and s19006). The text below is the original how-to for §1–§5.

Everything below runs from this folder with the gurobi-env python, and nothing in it uses Gurobi unless it
says so:

```bash
cd ~/RCK-cnp/results/bench_IW_allele
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
```

## 1. What is running

Four 12-task arrays, one task per case, one task at a time (`%1`: the Gurobi WLS licence has 2 sessions).
Submitted 2026-11-20 ~22:00 EST with `./submit_allelic.sh --submit`:

| job | method | role | time limit per task | starts |
|---|---|---|---|---|
| 32462471 | `rck_mc` | the method under test (multi-clone route, true clone profiles) | 15 min | at once |
| 32462470 | `rck_bulk` | **graph builder only**: CoRAL / CE decompose its graph; its own answers are not scored or reported | 90 min | at once |
| 32462472 | CoRAL 3.1.1 | baseline, decomposition-only on rck_bulk's graph | 90 min | case *i* after rck_bulk task *i* succeeds |
| 32462473 | CycleExtractor | baseline, as CoRAL | 60 min | case *i* after rck_bulk task *i* succeeds |

Expected: `rck_mc` a few minutes in all (1–5 s per solve), `rck_bulk` minutes to an hour, CoRAL / CE a few
hours in all (bench_IW_intpsi: ~1–2 h per baseline array). Tasks 1 and 8 of `rck_mc` (the pilot) skip as done.

Where things are:
- task logs: `/n/fs/yubi-iw/runs/bench_IW_allelic/logs/slurm/iwa-<method>_<job>_<task>.out`
  (`submitted.txt` there lists the job ids)
- runs: `/n/fs/yubi-iw/runs/bench_IW_allelic/methods/<case>/<method>/` (CoRAL / CE in `.../<method>/amplicon1/`),
  each with `run.json` (command, attempts, elapsed, status)
- task index → case: `python allelic_bench.py list --method rck_mc` (the same order for every method)

## 2. Is it done?

```bash
python allelic_bench.py status
```

This reads the run folders only (no Slurm calls). Done looks like `done 12` for each of the four methods;
anything else is listed with its index, case and state. If you also want Slurm's view, wait at least 60 s
between `squeue` / `sacct` calls (the controller rate-limits):

```bash
sacct -j 32462470,32462471,32462472,32462473 -X --format=JobID,JobName,State,Elapsed,MaxRSS
```

## 3. If a task failed

Read the task log and `<run>/run.json` first; don't resubmit blind. The usual causes:

| symptom | cause | what to do |
|---|---|---|
| log ends with repeated Gurobi licence errors, `run.json` attempts 15 | 4 tasks asked for the 2 licence sessions for longer than the 15 × 2 min retries (CoRAL / CE join while rck_bulk / rck_mc run) | resubmit (below) once fewer arrays are running |
| CoRAL / CE tasks stay `PENDING (DependencyNeverSatisfied)` | that case's `rck_bulk` task failed | fix and rerun the rck_bulk task, then `scancel` the stuck tasks and resubmit them |
| `TIMEOUT` in `sacct` | the run needed longer than the array's time limit | rerun that task with a longer `--time` |
| a Python error in the log | a bug or a bad input | look at the traceback before anything else |

Resubmitting is a new Gurobi submission (2 sessions: check that no other array of yours holds them).
`submit_allelic.sh` skips runs that already finished, so resubmitting everything only redoes what is missing:

```bash
./submit_allelic.sh            # prints the four sbatch commands
./submit_allelic.sh --submit   # submits them
```

To redo one task only, take its line from the printed commands and change `--array` to the task index, e.g.
`--array=5` (CoRAL / CE without the `--dependency` if that case's rck_bulk run is already done).

## 4. Score and collect (local, seconds, no Gurobi)

```bash
python allelic_bench.py score --jobs 8     # every finished run -> <run>/score*/ (cached; delete to rescore)
python allelic_bench.py collect            # -> bench_iw_allelic_table.tsv, _pooled.tsv, _proportions.tsv here
python sibling_relation.py                 # -> sibling_relation.tsv here (rck_mc)
```

These are bench_IW_intpsi's scores, unchanged (`results/bench_IW/bench_iw.py` docstring). In
`bench_iw_allelic_table.tsv`, one row per case × method variant:
- **Read the `rck_mc` rows and the CoRAL / CE rows.** CoRAL / CE appear twice, `walks_kept` and
  `walks_dropped` (their linear walks counted as species or not). Skip the `rck_bulk` rows: it is only the
  graph builder.
- Useful columns: `n_missed` / `n_false` (species), `phased_exact_cycles`, `phased_jseg_best_mean` /
  `phased_jadj_best_mean` (segment / adjacency Jaccard with haplotypes), `species_detection_strict_phased_f1`,
  `cn_recall` / `cn_precision` (ecDNA copies), `proportion_tv`. `rck_mc` is scored within phase blocks
  (`phase_frame` = `within_block`); `n_missed_absolute` / `n_false_absolute` are the absolute-frame counts.
- Each `rck_mc` run also has `evaluation_psi.tsv` (ψ per clone, truth vs inferred) and `recovery_matrix.tsv`
  (per truth species; `*_blockwise` columns are within phase blocks).

**`sibling_relation.py` is the check that matters most for this benchmark.** The within-block scores choose
the A/B swap per phase block for each truth species separately, so they could forgive a sibling put on the
outer species' allele. The script compares the sibling-vs-outer relation (opposite / same haplotype) on every
shared base, truth vs answer, which needs no swap. Good looks like `disagree_bp` 0 with `matched` True. It
already gives 0 bp disagreeing on the 3 finished cases (s17003, s17004, s19002).

## 5. Still to be written (item 1 scoring in PLAN.md)

Not done yet; worth asking Claude to do once the runs are scored ("resume item 1 from NEXT_STEPS.md §5"):
1. **Tables by scenario and stratum**: mean over the 4 cases of each scenario, and clean vs reuse (n = 2 each,
   descriptive only).
2. **`sib_opp_mirror` vs `sib_opp_shared`**: the allele model with and without the clone-wise allele signal in
   the input (the control's purpose).
3. **ψ per clone** for `rck_mc`: from each run's `evaluation_psi.tsv`, absolute error per species and clone.
4. **Proportions matched with haplotypes** (`evaluate._match(..., with_hap=True)`); `bench_iw_allelic_proportions.tsv`
   matches haplotype-blind.
5. **Part A's metrics on the new set** (`rck_style.py`: labelled-adjacency FNR / FDR, Eq. 11, `is_*`), and the
   post-hoc phasing comparison (`allele_bench.py`). Both are wired to bench_IW_intpsi today and need the
   `iw_allelic` set added, as `allelic_bench.py` does for bench_iw.
6. **Species matrix and figures** (`make_species_matrix_figure.py`, a figure per scenario), and the results
   written into PLAN.md.

## Reference

| file | what |
|---|---|
| `allelic_bench.py` | bench_iw pointed at set `iw_allelic`: `list`, `run`, `status`, `score`, `collect` |
| `allelic_method.sbatch`, `submit_allelic.sh` | the arrays (one method per array, `%1`) |
| `sibling_relation.py` | the per-base sibling-vs-outer relation check |
| `item1_results.py` | §5 tables (`item1/`): by scenario / stratum, mirror vs shared, ψ per clone, haplotype-matched proportions, post-hoc and Part A summaries (rck_mc, CoRAL, CE) |
| `make_item1_figures.py` | species-coverage figure per scenario (`figures/item1_species_*`) |
| `allele_bench.py --set iw_allelic`, `rck_style.py --set iw_allelic` | post-hoc phasing and Part A metrics on the item 1 set (`*_allelic.tsv`) |
| `make_topology_plots.py --truth RUN_DIR ...` | truth figures of simulator runs |
| amp `simulations/configs/bench_IW_allelic/` | configs, `seeds.tsv` (accepted cases, `ext_stratum`), `dryrun.tsv` |
| amp `simulations/subgraphs/bench_IW_allelic{,_cnp}/` | the bulk and multi-clone subgraphs of the 12 cases |
| `/n/fs/yubi-iw/runs/bench_IW_allelic/runs/` | simulator runs (truth) |
