# bench_IW_allele, item 1 (allelic siblings): following up after the method arrays

*Written 2026-11-20 when the arrays were submitted. The design and the results so far are in
[`PLAN.md`](PLAN.md) ("Item 1 in detail", "Resume here"). This file is the how-to: check the arrays, fix
failures, score, and what is still to be written.*

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
| `make_topology_plots.py --truth RUN_DIR ...` | truth figures of simulator runs |
| amp `simulations/configs/bench_IW_allelic/` | configs, `seeds.tsv` (accepted cases, `ext_stratum`), `dryrun.tsv` |
| amp `simulations/subgraphs/bench_IW_allelic{,_cnp}/` | the bulk and multi-clone subgraphs of the 12 cases |
| `/n/fs/yubi-iw/runs/bench_IW_allelic/runs/` | simulator runs (truth) |
