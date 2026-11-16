# bench_IW: status and next steps (branch `yubi_CNP_ionic`)

*First written 2026-11-14 at the end of the session that built bench_IW; "Resume here" rewritten
the same evening, after the methods array was checked and Yubi chose whole-number psi. Changes
and their reasons are in the `[yubi_CNP_ionic]` sections of both repos' `CHANGELOG.md` (entries
dated 2026-11-14, evening).*

## Resume here (next session)

> **Latest (2026-11-15): bench_IW_intpsi is finished, scored, collected and written up** (steps 5–7).
> All 96 runs done; results in `RESULTS_intpsi.md` and the three `bench_iw_intpsi_*.tsv` here.
> Open items for Yubi are under "After step 7" at the end of the next-steps list. Manuscript figures
> and Table 1: `make_figures.py` and `RESULTS.md`; what reviewers will likely ask for is in
> "Additions reviewers will likely ask for". Everything is
> committed on `yubi_CNP_ionic` (amplicon-allele-analysis unchanged at `0dc2d69`). Nothing pushed.

### Check the arrays (done 2026-11-15)

**Result of the check (2026-11-15, 13:00).**
- All 96 tasks COMPLETED with exit code 0; `status`: 96 done. Wall time 21:55–00:04 (2.1 h).
- No time-outs, out-of-memory kills or cancellations; peak MaxRSS 796 MB (coral task 9).
- No CE crashes. CE used at most 4 Gurobi threads (117 solves at 4, 46 small ones at 1).
- All 48 CoRAL / CE `run.json` have `graph_source: bulk+paths`, `path_constraints_dropped: 0` and
  `baseline: decomposition-only`.
- `rck_mc`: all 364 solves optimal, longest 4.2 s.
- Stopped at their own time limit, still scored on the incumbent: `rck_bulk` decomposition on
  hard_S4_s16003 (gap 79%; 85% on bench_IW) and CoRAL on simple_S4_s13002 (gap 9%).
- **Different from expected:**
  - Licence contention was heavy. 42 of 96 runs needed more than one attempt (up to 9; none
    used all 15). Gurobi reported "5 active sessions for a baseline of 2", so something
    outside these arrays also held sessions.
  - `rck_mc` simple_S4_s13004 took 731 s of its 15 min limit, almost all waiting for the licence.
  - Summed `elapsed_s` (licence waits included): `rck_bulk` 2.04 h, `coral` 1.98 h, `rck_mc`
    1.57 h, `ce` 1.23 h (6.8 h, against about 4 h expected).
- Nothing resubmitted.

The original instructions, kept for reference:


| array | method | Slurm job | tasks | time limit | depends on |
|---|---|---|---|---|---|
| A | `rck_bulk` | **32280369** | 0–23 %1 | 90 min | – |
| B | `rck_mc` | **32280370** | 0–23 %1 | 15 min | – |
| C | `coral` | **32280371** | 0–23 %1 | 90 min | `aftercorr:32280369` |
| D | `ce` | **32280372** | 0–23 %1 | 60 min | `aftercorr:32280369` |

- Task i = case i of `python bench_iw.py --set iw_intpsi list --method <m>`. Task 0
  (simple_S1_s11001) and task 23 (hard_S4_s16006) were run locally in the pilot and exit at once
  ("done already").
- Submitted 21:55 with `./submit_intpsi.sh --submit` (`/n/fs/yubi-iw/runs/bench_IW_intpsi/logs/slurm/submitted.txt`).
  Expected about 2–2.5 h of wall time.
- Logs: `/n/fs/yubi-iw/runs/bench_IW_intpsi/logs/slurm/iwi-<method>_<jobid>_<task>.out`. Runs:
  `/n/fs/yubi-iw/runs/bench_IW_intpsi/methods/<case>/<method>[/amplicon1]/run.json`.

```bash
cd ~/RCK-cnp/results/bench_IW
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
squeue -u $USER                                            # anything left? (>= 60 s between calls)
sacct -j 32280369,32280370,32280371,32280372 -X --format=JobID%20,JobName,State,ExitCode,Elapsed,MaxRSS \
  | awk 'NR>2{print $2, $3}' | sort | uniq -c              # states per array
sacct -j 32280369,32280370,32280371,32280372 --format=JobID%20,State,ExitCode,Elapsed,MaxRSS \
  | grep -E 'FAILED|TIMEOUT|OUT_OF_ME|CANCELLED|PENDING'   # anything not COMPLETED
python bench_iw.py --set iw_intpsi status                 # done / missing / failed per (case, method)
grep -l 'licence' /n/fs/yubi-iw/runs/bench_IW_intpsi/logs/slurm/*.out | wc -l   # tasks that waited for the licence
```

What to look for, and what to do:
- **All 96 `done` in `status`:** go to step 6 (score and collect).
- **A failed `rck_bulk` task:** its case's CoRAL / CE tasks stay PENDING with reason
  `DependencyNeverSatisfied`; cancel them (`scancel 32280371_<i> 32280372_<i>`).
  - Diagnose from its `.out` and `<case>/rck_bulk/run.json` first.
  - Then resubmit only the affected tasks with the same flags as `submit_intpsi.sh`, e.g.
    `sbatch --array=<i>%1 ... --export=ALL,SET=iw_intpsi,METHOD=rck_bulk bench_iw_method.sbatch`,
    and CoRAL / CE for that case with `--dependency=aftercorr:<new bulk job>`. Finished runs are
    skipped. Ask Yubi before resubmitting more than a handful.
- **Licence failures** (`attempt k failed (licence)` x 15 in the `.out`): resubmit those tasks;
  nothing to fix.
- **CE crashes** (`NameError` inside `CE.py`, as twice on bench_IW): a method failure. Record it, do
  not patch CE, do not resubmit.
- **TIMEOUT:** note which; do not raise limits without Yubi.
- **`MaxRSS`:** all should be far below 2G (peak 754 MB on bench_IW).
- Check that CE ran with 4 threads (`grep -h 'Thread count was' <case>/ce/amplicon1/ce.stdout`) and
  that every CoRAL / CE `run.json` has `graph_source: bulk+paths` and
  `input_graph_counts.path_constraints_dropped: 0`.

**State (2026-11-14, evening).**
- **bench_IW's methods array (32267060) is finished and checked.** It is not the benchmark's
  result: bench_IW has a mismatch between the multi-clone input and the scoring truth (below).
  - 92 of 96 tasks completed; no time-outs, no out-of-memory kills (peak 754 MB of 4 GB).
  - Failed:
    - CoRAL on simple_S1_s11004 and hard_S2_s15004 (tasks 14 and 74): the Gurobi licence was
      busy on all 5 attempts;
    - CycleExtractor on simple_S4_s13002 and hard_S2_s15005 (tasks 39 and 79): crashes inside
      CycleExtractor 1.0.0 itself (`NameError: Final_closed_walk`, `CE.py:3757`;
      `NameError: segments`, `CE.py:3954`). A resubmission fails the same way.
  - Stopped at the time limit: `rck_bulk` decomposition on hard_S4_s16003 (gap 85%); CoRAL on
    hard_S4_s16002 and simple_S4_s13002 (gap 9%).
  - `rck_mc`: every solve optimal, longest 3.95 s (hard_S4_s16003). This is the timing Yubi
    asked for before fixing the multi-clone time limit: **60 s stays**.
  - CycleExtractor starts Gurobi with 32 threads on a 4-CPU allocation.
  - Compute per method (summed run.json elapsed, licence waits included): `rck_bulk` 1.09 h,
    `coral` 1.87 h, `ce` 0.97 h, `rck_mc` a few minutes.
  - Not scored or collected.
- **Found: the multi-clone input and the scoring truth disagreed.**
  - psi per clone is real-valued, and the reads and `bench_metrics`' truth use it
    (PHI_s = sum_k u_k psi_k).
  - The clone karyotypes in `truth/rck/`, which are `--clone-profiles truth`, fold in
    round(psi).
  - So `rck_mc` could at best report sum_k u_k round(psi_k), up to 4.7% from PHI_s.
- **Decision (Yubi): whole-number psi by construction (option A).** The rejected alternatives
  and the reasons are in both CHANGELOGs.
- **Done for A (amplicon-allele-analysis, `yubi_CNP_ionic`):**
  - simulator flag `l3_mixture.integer_psi` (commit `2363469`, 7 tests);
  - set `iw_intpsi`;
  - **bench_IW_intpsi**: 24 runs re-simulated from bench_IW's own seeds, both subgraphs built
    (commits `0a164dc` data, `0dc2d69` docs).
  - All 48 subgraphs pass their checks; `check-shared`: 24 runs, 0 differ.
  - `truth/rck/` copy numbers, which are the multi-clone input, are byte-identical to bench_IW's:
    round(psi) did not change. What moved to match them is the reads and the scoring truth.
  - One 1-bp DIAG sliver differs in hard_S2_s15005 (chr5:152,315,000; CHANGELOG). Not changed.
  - The same 24 seeds were accepted. Structures, clones, `segment_max_traversals` and species
    sizes are identical to bench_IW; only psi, PHI and the reads changed.
  - psi is whole (3–85) and PHI_s = sum_k u_k psi_k in every run.
  - Details: that repo's CHANGELOG and `simulations/README.md` 5.5.
- Nothing pushed. Yubi pushes.

**Decisions that apply from here** (Yubi's; they override the earlier text below where they differ):
1. Clone-level ecDNA copies are whole numbers by construction (bench_IW_intpsi). The scoring
   truth stays PHI_s = sum_k u_k psi_k, which is now also what the multi-clone input carries.
2. CoRAL and CycleExtractor run on **RCK bulk's exported graph** (`<case>/rck_bulk/aa/*_graph.txt`),
   not on the union subgraph exported directly as array 32267060 did.
   - Label them *decomposition-only baselines*.
   - Each case's graph-route tasks depend on that case's bulk task.
3. The results table has a **`strict_rule` column** (`Status` or `cn-frac:0.1`).
   `bench_metrics.json` already records it; check that `collect` carries it.
4. **Scores:**
   - Headline: ecDNA-only, haplotype-blind, bp-weighted copy-number recall.
   - For `rck_bulk` and `rck_mc` also **total copy-number recall** (background plus ecDNA).
   - Proportions over strict cycles, with the `low_weight`-inclusive version too.
5. **Per run, record** `segment_max_traversals`, species sizes and the strict rule.
6. **Multi-clone time limit:** 60 s, from the timing above.
7. **Multi-clone species counts: within phase blocks** (2026-11-14 night, after the pilot).
   - `rck_mc`'s strict pairing, `n_missed` / `n_false`, M1, composites and phased species matrix
     are scored within phase blocks (`aa_eval --phase-blocks <DIAG>/rck.fragments.tsv`;
     `docs/pipeline/evaluation.md` §3.10): one A/B swap per block, shared by all species.
   - *Why:* closest to the original RCK, whose ILP leaves one allele flip free per fragment
     (`rck/core/ilp_gurobi.py`, `fragment_flipping`). The multi-clone DIAG's fragments are the
     phase blocks, so labels are identified only within them.
   - The absolute counts stay as `n_missed_absolute` / `n_false_absolute`; blind columns
     unchanged. The swap is chosen to agree best with the truth (state this when reporting).
   - `rck_bulk` keeps absolute scoring.
8. **CoRAL / CE linear walks: report both, side by side** (2026-11-14 night).
   - Each answer is scored with its linear walks kept (`method_variant` `coral:walks_kept`,
     `score/`) and dropped (`coral:walks_dropped`, `score_nolinear/`).
   - *Known limit of both:* on the bulk route's total-copy graph the background can also come back
     as a closed cycle (pilot: 2 copies over 45.6 Mb in hard_S4_s16006), which neither variant
     removes.
   - AmpliconClassifier (option d) was not adopted; it is not installed here.

**Rules.**
- Small local runs are fine. New Slurm arrays need Yubi's go-ahead.
- The licence allows 2 Gurobi sessions: never run local Gurobi jobs while arrays hold both.
- Wait at least 60 s between `squeue` / `sacct` calls.
- Do not modify bench_IW's configs, seeds, runs or array outputs (kept for the record).
- Do not touch `~/RCK` or `~/amplicon-allele-analysis` (`yubi_dev`).

**Next steps, in order:**
1. **Done (2026-11-14 night, `9620a8a`): driver points at bench_IW_intpsi.** `bench_iw.py --set
   iw_intpsi` switches the case list, simulator runs, subgraphs, method outputs
   (`/n/fs/yubi-iw/runs/bench_IW_intpsi/methods`) and the results TSVs (`bench_iw_intpsi_*.tsv`).
   Default `iw` unchanged.
2. **Done (2026-11-14 night): CoRAL / CE on the bulk route's graph, as decomposition-only
   baselines.** Details and checks: CHANGELOG, "bench_IW_intpsi driver; graph route on the bulk
   route's graph".
   - Graph: `rck.ecdna.aa_export.baseline_graph`, the bulk run's `aa/*_graph.txt` (+1 shift,
     fitted copy numbers, every edge it used) plus the direct subgraph export's long-read path
     constraints, renumbered by coordinate.
   - *Why the path constraints* (Yubi chose to keep them, to stay close to the original tools):
     CoRAL uses them by default as MILP constraints, CE reads them, and the bulk route reads the
     same paths.
   - Written once per case to `<case>/graph_bulk/`. A coral / ce run exits 2 without running
     until the case's `rck_bulk` is done.
   - `run.json`: `baseline: decomposition-only`, `graph_source`, the graph's counts.
   - CE's Gurobi threads are capped by a `gurobi.env` in its run folder; **to confirm in the
     pilot.**
   - Licence errors: up to 15 retries, 2 min apart.
   - `bench_iw_method.sbatch` + `submit_intpsi.sh` (prints the four commands; `--submit` submits).
   - Offline check on bench_IW's 24 bulk runs: all 389 path constraints carried, none dropped;
     3,391 tokens point at the same coordinates.
3. **Done (2026-11-14 night): total copy-number recall for the RCK routes, per-run fields.**
   Details: CHANGELOG, "total copy-number recall (RCK routes) and per-run fields".
   - `total_cn_*` over the union subgraph's segments, and `total_cn_ecdna_*` over the truth ecDNA
     bases.
   - Truth frame: cell-averaged tumour sample, normal cells included. `rck_bulk` from
     `flow_segments.tsv`; `rck_mc` = sum_k u_k (A_k + B_k) + normal share.
   - **When reporting:** `total_cn_*` over the whole subgraph is mostly background and purity (one
     whole-number bulk karyotype cannot match 1.2-type cell-averaged values). Lead with
     `total_cn_ecdna_*`.
   - The multi-clone total mostly reflects its given profiles, except where the model moves them
     (hard_S4_s16003 on bench_IW).
   - Table: `segment_max_traversals`, `footprint_bp`, `species_size_bp` (+ min / max),
     `strict_rule`, `baseline`.
4. **Done (2026-11-14, 21:30): local pilot on bench_IW_intpsi.** Cases simple_S1_s11001 and
   hard_S4_s16006 (the cheapest hard S=4 case in array 32267060), the four methods one at a time
   (at most one Gurobi session, 4 threads, `nice`).
   - The 8 runs are in `/n/fs/yubi-iw/runs/bench_IW_intpsi/methods/` and count as done; the
     arrays skip them.
   - Logs: `/n/fs/yubi-iw/runs/bench_IW_intpsi/logs/pilot_<case idx>_<method>.log`.

   | case | method | time | cn_recall | cn_precision | proportion_tv | missed / false | blind J_seg / exact | total_cn_ecdna_recall |
   |---|---|---|---|---|---|---|---|---|
   | simple_S1 | rck_bulk | 3 s | 1.00 | 0.91 | 0 | 0 / 0 | 1.00 / 1 | 1.00 |
   | simple_S1 | rck_mc | 2 s | 1.00 | **1.00** | 0 | 0 / 0 | 1.00 / 1 | 1.00 |
   | simple_S1 | coral | 53 s | 1.00 | 0.19 | 0.17 | 0 / 1 | 1.00 / 1 | – |
   | simple_S1 | ce | 1 s | 1.00 | 0.19 | 0.29 | 0 / 2 | 1.00 / 1 | – |
   | hard_S4_s16006 | rck_bulk | 135 s | 0.85 | 1.00 | 0.25 | 1 / 0 | 0.62 / 1 | 0.90 |
   | hard_S4_s16006 | rck_mc | 3 s | 1.00 | 1.00 | 0 | 2 / 2 (phased) | 1.00 / 2 | 1.00 |
   | hard_S4_s16006 | coral | 13 s | 0.77 | 0.23 | 0.55 | 2 / 0 | 0.27 / 0 | – |
   | hard_S4_s16006 | ce | 3 s | 0.92 | 0.26 | 0.56 | 0 / 2 | 0.67 / 0 | – |

   - **Infrastructure checks passed:**
     - CE now uses at most 4 Gurobi threads (`gurobi.env`; it used 32 before);
     - CoRAL / CE read `graph_bulk/` with every path constraint carried (4/4, 26/26);
     - `run.json` and the table carry `baseline`, `strict_rule`, traversals and sizes;
     - CE did not crash;
     - scoring 8 runs took 2 s.
   - **Whole-number psi works as intended:** on simple_S1, `rck_mc` ecDNA precision is now 1.00
     (0.96 on bench_IW, where its input carried round(psi) and the truth did not). Its psi
     evaluation is exact (8 / 4).
   - **Found: phased matching understates `rck_mc`.** On hard_S4_s16006 its answer is correct
     haplotype-blind: each truth species has one cycle with J_seg 1.00 and the same weight, and
     copy-number recall, precision and proportion_tv are exact.
     - But `aa_eval`'s strict matching (and so `n_missed` / `n_false`, M1) and `rck-ecdna`'s own
       psi evaluation (`evaluation_psi.tsv`) match by **phased** Jaccard.
     - The multi-clone route's haplotype labels are meaningful only within phase blocks, so
       they pair species wrongly: 2 missed / 2 false, phased J_seg 0.42.
     - Not changed: a metric decision for Yubi (blind matching, or within-block phasing, for
       `rck_mc`'s species counts).
   - Times are far below the array limits (bulk 135 s on hard S=4).

5. **Submitted 2026-11-14 21:55 (Yubi's go-ahead): jobs 32280369 (rck_bulk), 32280370 (rck_mc),
   32280371 (coral), 32280372 (ce)**; see "Check the arrays" above. `./submit_intpsi.sh` prints
   the commands, `./submit_intpsi.sh --submit` submits. Four 24-task arrays, one task per case, 4 CPUs and 2 G
   each (peak 754 MB on bench_IW):

   | array | method | depends on | time limit | wall time (array 32267060 timings) |
   |---|---|---|---|---|
   | A | `rck_bulk` | – | 90 min (max 1908 s, licence waits included) | about 1.1 h at %1 |
   | B | `rck_mc` | – | 15 min (solves <= 4 s) | a few minutes |
   | C | `coral` | `aftercorr:A` (case i after bulk task i succeeds) | 90 min (max 2002 s) | about 1.9 h |
   | D | `ce` | `aftercorr:A` | 60 min (max 374 s; CE's own limit is 7200 s per MILP) | about 1 h |

   - **Gurobi licence (2 sessions):** A and B hold at most one each. C and D join as bulk tasks
     finish, so up to 4 may ask at once; a task waits up to 15 x 2 min.
   - **A failed bulk task** leaves that case's C / D pending (DependencyNeverSatisfied): cancel
     them, diagnose, and resubmit (finished runs are skipped).
   - **Expected total:** about 2–2.5 h wall time, about 4 h compute.
   - CoRAL / CE timings on the bulk route's graph may differ from these.
   - CE crashed on two bench_IW graphs; report crashes as method failures, do not patch CE.
6. **Done (2026-11-15): score and collect** (4 s; 144 rows, no `score_error.txt`, `strict_rule` and
   `baseline` carried). Original instructions:
   `python bench_iw.py --set iw_intpsi score --jobs 8`, then `python bench_iw.py --set iw_intpsi collect`.
   - Writes `bench_iw_intpsi_table.tsv`, `bench_iw_intpsi_pooled.tsv` and
     `bench_iw_intpsi_proportions.tsv` here.
   - Check for `score_error.txt` in any `score*/` folder (`collect` marks them `score_status=error`).
   - The pilot's 12 scored variants took 2 s, so all 144 rows (24 cases x 6 method variants)
     should take under a minute.
   - The table has one row per (case, `method_variant`); the pooled table groups by
     `method_variant` and gives `n` and `n_scored`.
   - Commit the TSVs with a short results README.
7. **Done (2026-11-15): `RESULTS_intpsi.md`** (pooled tables per method variant, class and S; findings).
   Original instructions:
   - a short `results/bench_IW/RESULTS_intpsi.md` with the pooled table per (method variant, class,
     S), committed with the TSVs;
   - lead with ecDNA copy-number recall / precision, proportions and the species matrix; give
     `total_cn_ecdna_*` before the whole-subgraph `total_cn_*` (purity-dominated);
   - failures;
   - CoRAL / CE behaviour, both walk variants side by side;
   - `rck_mc`'s within-block vs absolute counts.

**After step 7: open items for Yubi** (from `RESULTS_intpsi.md`; nothing started):
- **`rck_mc` on hard_S4_s16003:** copy-number recall 0.71 with psi within 1 of the truth. This is
  the case where the multi-clone model moved its given profiles on bench_IW too. Worth a look.
- **`rck_mc` extra cycles:** 4 simple cases have exact copy numbers but split a clone's copies
  over extra cycles (8 false species in all, proportion TV up to 0.40). Is this a decomposition
  tie-break or something to penalise?
- **Licence:** another holder used sessions during the arrays. For future arrays, consider a
  longer `rck_mc` time limit or fewer concurrent arrays.
- `elapsed_s` includes licence waits; the solvers' own runtimes are in each run's logs and
  `summary.json`, if timings are to be reported.

## Additions reviewers will likely ask for (manuscript; noted 2026-11-15, nothing started)

The manuscript's claim is that RCK-ecDNA reconstructs ecDNA better than CoRAL and CycleExtractor
(CE). The current evidence is `RESULTS.md` (figures 2–6, Table 1; `make_figures.py`) on
bench_IW_intpsi. Below is what a reviewer is likely to ask for, roughly in order of how much each
strengthens the paper, with the reasoning. Each needs Yubi's go-ahead before any new Slurm array.

1. **CoRAL and CE end to end, each building its own graph from the same reads.**
   - *Why:* today they only decompose RCK bulk's exported graph, which supports "better
     decomposition given the same graph", not "better method".
   - The obvious reviewer reply to figures 2 and 4 is that CoRAL was handed a total-copy graph
     with no ecDNA vs background separation, which is not its intended input. CoRAL builds its own
     breakpoint graph from long reads; CE normally runs on an AmpliconArchitect graph.
   - Running both from reads removes that objection, and keeps the current runs as a second,
     controlled comparison ("same graph").
   - *Where to start:* array 32267060 (bench_IW) already ran CoRAL / CE on the union subgraph
     exported directly. CoRAL's own graph-building step is in the `coral-3.1.1` env. AA for CE
     would need installing; check `module avail` first.
   - *Score with the same* `bench_iw.py score`, so the tables line up.
   - *How:* outlined step by step in "Addition 1 in detail" below.
2. **RCK multi-clone with estimated, not true, clone profiles.**
   - *Why:* `rck_mc` currently gets the true clone profiles (`--clone-profiles truth`), so its
     numbers are an upper bound and a reviewer will discount them.
   - Clone profiles from a standard tool (e.g. HATCHet, or whatever the CNP pipeline in
     amplicon-allele-analysis produces), or perturbed truth profiles (wrong purity, wrong clone
     proportions, ±1 copy noise), show how much of the advantage survives realistic input.
   - A curve of the metric against profile error is the most informative form.
   - Without this, RCK bulk is the only like-for-like comparison with CoRAL / CE. It still wins,
     but the clone-aware route is the paper's novelty.
3. **Robustness sweep: read depth, tumour purity, read length / error, breakpoint noise.**
   - *Why:* one simulator setting proves little about real data.
   - Reviewers expect metric-against-parameter curves (copy-number precision, J_seg, proportion
     TV), with all methods on the same cases.
   - The simulator already exposes depth and purity. Keep the S x class design small (e.g. S = 2
     and 4, nested only) to bound the run count.
   - Watch the 2-session Gurobi licence: this is the most expensive addition.
4. **At least one real tumour dataset.**
   - *Why:* an all-simulation methods paper is a common rejection reason.
   - Candidates: a cell line with well-characterised ecDNA and long reads (e.g. COLO320-DM,
     GBM39 or PC3 have published long-read / FISH / CRISPR-CATCH structures).
   - Ground truth will be partial, so the comparison is qualitative:
     - agreement with published structures;
     - long-read support for the reconstructed junctions;
     - whether the alleles RCK calls agree with heterozygous SNPs on the reads.
   - That last check is the evaluation unique to an allele-specific method; make it the centre
     of the real-data section.
5. **A larger, more varied benchmark.**
   - *Why:* 24 cases (4 per class and S, S <= 4, one shared set of clones) gives wide error bars
     and weak per-class tests (`table1_paired_tests.tsv`).
   - More replicates (8–10 per cell), S up to 6–8, several clone sets, and other structures
     (e.g. ecDNA with linear amplicons in the same region) would make the per-class claims
     testable.
   - Cheap for RCK; the baselines and the licence set the cost.
6. **Explain or fix the two RCK multi-clone weak spots** (already open, see "After step 7").
   - hard_S4_s16003: recall 0.71.
   - Extra cycles with exact copy numbers in 4 simple cases. These lower its species F1 below RCK
     bulk's (0.70 vs 0.76) and its simple-case exact-cycle rate (0.75 vs 0.96).
   - A reviewer will spot both in the per-case points of figures 2, 5 and 6. Either fix them
     (e.g. a parsimony term or tie-break that prefers fewer cycles per clone) or explain them in
     the text.
7. **Runtime and scaling.**
   - *Why:* methods papers are expected to report compute.
   - Current `elapsed_s` includes Gurobi licence waits, so it can't be used. Collect solver times
     from each run's logs / `summary.json`, and time a few cases with an exclusive licence.
   - Report time against graph size (segments, breakpoints), and note the solver time limits and
     the two runs that hit them.
8. **Ablations of RCK-ecDNA's own components.**
   - *Why:* reviewers ask which part of the model makes the difference.
   - Candidates:
     - allele-specific vs allele-blind graph;
     - with and without the long-read path constraints;
     - bulk vs multi-clone on the same input;
     - the opt-in model options from earlier work: `--depth-bias-*` (CHANGELOG) and
       `--unexplained-junction-penalty` (commit `23b5ecf`).
   - These are cheap local or small-array runs on bench_IW_intpsi.
9. **Fair-use statements for the baselines.**
   - *Why:* to pre-empt "the baselines were misconfigured".
   - Record versions (CoRAL 3.1.1, CycleExtractor 1.0.0), parameters, time limits, thread counts
     and the path-constraint handling in a supplementary table, all already in each `run.json`.
   - State that CE crashed on 2 of 24 bench_IW graphs inside CE itself (`NameError`, `CE.py`),
     and that CE was not patched.

### Addition 1 in detail: CoRAL and CE end to end on bench_IW_intpsi (outline, 2026-11-15)

> **Started 2026-11-15 evening:** route (B); status and next steps in
> [`results/bench_IW_e2e/NEXT_STEPS.md`](../bench_IW_e2e/NEXT_STEPS.md). Two refinements to the outline
> below: reads cover each ecDNA interval ± 5 Mb plus a 2 Mb window per chromosome arm (not the whole
> genome), and `--rck-junction-reads-per-copy` is dropped for the end-to-end RCK runs.

**Goal.** Score CoRAL and CE when each one builds its own graph from the same simulated reads.
Compare that with RCK-ecDNA's own end-to-end result on the same 24 cases, against the same truth,
with the same scorer. The current runs (CoRAL / CE decomposing RCK bulk's graph) stay in the
paper as the controlled "same graph" comparison.

**What blocks it today (checked 2026-11-15).**
- **No reads exist.** bench_IW_intpsi runs hold only `obs/` count tables and `truth/`. The
  simulator's read layer L5 is not written: `simulations/src/simulation.py::run_reads` raises
  `NotImplementedError`, and `emit_reads`, `l5_reads` and `--reference-fasta` are unused
  (amplicon-allele-analysis `simulations/README.md` §S10). RCK-ecDNA reads the count tables of
  L4 directly. CoRAL needs a long-read BAM; AA needs a short-read BAM.
- **No read tools on ionic.** `module avail` and `PATH` have no minimap2, samtools, NanoSim,
  CNVkit, ART / Mason or AmpliconSuite. The `coral-3.1.1` env has `coral` only. No hg38 FASTA is
  configured.
- **CE needs a graph.** It decomposes an AA (short read) or CoRAL (long read) `*_graph.txt`
  (`/n/fs/yubi-iw/tools/CycleExtractor/README.md`). On long reads, "CE end to end" therefore
  means CoRAL's graph followed by CE's decomposition. CE on an AA graph is a separate short-read
  arm (optional, step 8).

**Decision for Yubi: how far "end to end" goes for RCK-ecDNA.**
- **(A) Baselines from reads, RCK-ecDNA from the simulator's counts.**
  - Least work.
  - A reviewer will object: RCK-ecDNA gets error-free breakpoint calls, exact junction-read
    counts and the true SNP phasing, while CoRAL / CE pay for alignment and calling.
  - Not recommended as the main comparison.
- **(B) Everything from the same BAMs (recommended).**
  - RCK-ecDNA's inputs (`obs/`: depth windows, SNP BAF, adjacency / extremity / path tables,
    `#NOV_A/B`) are rebuilt from the BAMs. The tools: the callers in
    `prepare_adjacency_02/` (Severus, Sniffles, GRIDSS2; `run_phasing.sh`), `genotype_snps.sh`,
    and the UGP window layout that `scripts/read_ugp_outputs.py` reads.
  - Then `subgraph_construct.py` and `rck-ecdna` run unchanged.
  - Not yet checked: whether this real-data path writes `path.tsv` and the allele-split
    junction counts. If it does not, a converter from phased long reads is needed (step 4).
  - `rck_mc` still takes the true clone profiles (`--clone-profiles truth`). In (B) it should
    be reported with addition 2 (estimated profiles); otherwise mark it as an upper bound.
- Report (A) as well if (B) is done: the gap between (A) and (B) is RCK-ecDNA's own loss from
  read processing.

**Steps.**
1. **Tools and reference** (one-off; ask before installing).
   - A new conda env under `/n/fs/yubi-iw/conda/envs/` with minimap2, samtools, NanoSim (with
     ecSimulator's `promethion_10_4_1_guppy_model`), CNVkit and bedtools. For (B), also the
     `prepare_adjacency_02/env_*.yaml` envs.
   - A chr-prefixed GRCh38 FASTA with its minimap2 index (about 3 GB plus about 8 GB) in
     `/n/fs/yubi-iw/ref/`. The simulator's coordinates are GRCh38 (`simulation.yaml`
     `reference`), and ecSimulator ships `GRCh38_centromere.bed`, which `coral seed` needs.
2. **Write L5** (amplicon-allele-analysis `simulations/src/simulation.py::run_reads`, as planned
   in README §S10), as a new command that adds reads to an existing run without re-running
   L0–L4.
   - The truth tables stay byte-identical, so the current scores remain valid.
   - Seed it from the run's `random_seed`.
   - **Sequences.**
     - Apply the phased germline SNPs (`truth/phasing.tsv`, `obs/windows/snp.tsv.gz`) to the
       reference to get haplotypes A and B over the run's regions plus `flanking_length`.
     - Build each clone's linear karyotype from `truth/segments.tsv` (per-clone $(a_{i,j},
       b_{i,j})$) and each ecDNA species' circular sequence from `truth/walks.walk`, with the
       haplotype of every step. ecSimulator's `cycles_file_to_fasta.py` handles the walk part;
       the haplotype per step is new.
   - **Coverage.** Use the count layer's depth per copy (`l4_counts`; γ from `library_size`):
     clone $i$'s karyotype at $u_i \gamma$ per copy, and species $s$ at
     $\gamma \sum_i u_i \psi_{s,i}$, the same $\psi$ as `truth/proportions.tsv`. Generate the
     matched normal the same way (diploid, `normal_library_size`). This is how ecSimulator's
     `run_nanosim.py` scales amplicon against background (NanoSim `genome`, `-dna_type circular`
     for walks), done per sequence.
   - **Scope.** Simulate reads only over the run's regions plus flanks: the grid already has
     `restrict_to_regions: true`, and whole-genome coverage would be about 40 Gbp per case. Do
     align to the whole of GRCh38, so mapping ambiguity stays realistic.
   - **Check the reads against the counts on 2 pilot cases** (one simple_S1, one hard_S4):
     depth per copy, junction-spanning reads per adjacency against `obs/adjacency.tsv` `#NOV`,
     and SNP BAF. They should agree within sampling noise. If they do not, L5 and L4 disagree,
     and the comparison is not apples to apples.
3. **Align and call copy number** (per case, Slurm).
   - `minimap2 -ax map-ont` → `samtools sort/index` for the tumour and the normal BAM.
   - CNVkit with the matched normal as the reference, restricted to the simulated regions
     (`--targets` the region BED, `-m wgs`). The output `.cns` is CoRAL's `--cn-seg`.
   - Keep BAMs and intermediates in `/tmp` during the job, then write one compressed copy of
     each BAM under `/n/fs/yubi-iw/runs/bench_IW_intpsi/reads/<case>/`. That is a few large
     files per case, not many small ones.
4. **(B) only: rebuild RCK-ecDNA's inputs from the BAMs.**
   - Call SVs (Severus / Sniffles on the long reads), phase SNPs, and genotype the SNPs.
   - Write `obs/` in the simulator's schema (`build_adjacency_graph/schema.py`) into a new run
     folder, `runs/bench_IW_intpsi_reads/<case>/`, with `truth/` symlinked.
   - Then `bench_IW.py [--cnp] rebuild` builds the subgraphs and the DIAGs as now.
5. **CoRAL end to end** (method `coral_e2e`).
   - `coral seed --cn-seg <cns> --lr-bam <bam> --centromere-file GRCh38_centromere.bed`
     at defaults (`--gain 6`, `--min-seed-size 100000`).
   - Then `coral reconstruct --lr-bam --cnv-seed --cn-seg` with `--solver gurobi_direct`, the
     same 1800 s per MILP as today, and the default `--cycle-decomp-mode max_weight`.
   - Defaults are the main result, for fair use. The bench's smallest species are 300 kb with 4
     copies, below the `--gain 6` seed threshold. Record which truth species no seed covers
     (their recall is 0 for a seeding reason, not a decomposition reason), and add one
     sensitivity run with `--gain` lowered, reported separately.
6. **CE end to end** (method `ce_e2e`): `CE.py` at its defaults on each `*_graph.txt` from step
   5, path constraints included. Same Gurobi `gurobi.env` thread cap as `bench_iw.py run_ce`.
7. **RCK-ecDNA end to end** ((B) only; methods `rck_bulk_e2e`, `rck_mc_e2e`): today's
   `rck_argv` flags on step 4's subgraphs.
8. **Optional short-read arm for CE's native input.**
   - ART / Mason short reads (`l5_reads.short_simulator`), then AmpliconSuite-pipeline (AA, with
     its own CNVkit seeds), then CE on the AA graph.
   - Needs the AA data repo and AA's licence requirements; check them before committing to it.
     Only worth it if a reviewer is expected to say "CE is meant for AA graphs".
9. **Driver and scoring** (RCK `results/bench_IW/bench_iw.py`).
   - Add the new methods, with their own graph sources, to `METHODS` and `GRAPH_SOURCES`, and a
     set, e.g. `--set iw_intpsi_e2e`, that reuses bench_IW_intpsi's cases and truth.
   - CoRAL reconstruct can report several amplicons per case. Pool their cycles into one answer
     per case before scoring, as `results/graph_route/score_h2h.py` does per (case, amplicon),
     and say so in the methods text.
   - Score with `bench_iw.py score`, unchanged: `aa_eval` (species matrix, M1 / M2, strict and
     best-of) and `bench_metrics`; linear walks kept and dropped, as for `coral` / `ce` now.
   - **Add graph-level metrics.** Once the methods build their own graphs, build the truth mixture
     graph (`simulations/truth_mixture.py`; `aa_eval` reads it from `truth/aa/mixture_amplicon<N>_*` or `--mixture-dir`). `bp_accuracy`,
     novel-adjacency FNR / FDR and `segment_cn_*` then separate graph-construction errors from
     decomposition errors (`docs/pipeline/evaluation.md` §4.2–§4.3). Today's runs have none.
10. **Arrays.**
    - Per method, chained with `--dependency=afterok`, in this order:
      - reads (step 2) → align + CNVkit (step 3) → (B) inputs (step 4);
      - `coral_e2e` → `ce_e2e`;
      - `rck_*_e2e`.
    - Read simulation and alignment need no Gurobi and can run wide. The Gurobi steps share the
      2-session licence: keep `%1`–`%2` per array and the 15-try licence retry. Size memory from
      the 2 pilot cases (`seff`).
    - Yubi's go-ahead before each array.
11. **Report.**
    - Figures 2–6 and Table 1 rows for `coral_e2e`, `ce_e2e` and (B)'s `rck_*_e2e`, next to the
      current "same graph" rows.
    - A supplementary table of seeds missed by `coral seed`, graph-level metrics per method, and
      every tool version and flag (from each `run.json`).

**Cost, roughly.**
- Per case: one NanoSim run over a ≤ 8 Mb footprint plus flanks; one minimap2 alignment
  (minutes); CNVkit (minutes); CoRAL reconstruct (bounded by its time limits).
- The licence, not CPU, sets the wall time, as in the methods arrays of 2026-11-14.
- Writing L5 and checking it (step 2) is most of the human time.

## Earlier state (2026-11-14 afternoon, kept for the record)

> Superseded by "Resume here" above. Array 32267060 ran on bench_IW, whose multi-clone input and
> scoring truth disagree; its outputs are kept but are not the benchmark's results.

**State.**
- **Simulation side done and committed** (amplicon-allele-analysis, `yubi_CNP_ionic`): 24 accepted
  runs, union subgraphs (bulk route) and multi-clone DIAGs. All 24 share demo7's clones
  (`configs/bench_IW/shared_clones.tsv`).
- **RCK side committed:** the metrics (`rck/ecdna/bench_metrics.py`), the driver
  `results/bench_IW/bench_iw.py` and `bench_iw.sbatch`.
- **Methods array submitted:** Slurm job **32267060** (96 tasks = 24 cases x 4 methods, `%2`
  for the Gurobi WLS licence, 4 CPUs, 4G, 3 h each).
  - When this file was written (15:32, 2026-11-14), 13 of 96 runs were done: the 4 local pilot
    runs on `simple_S1_s11001` plus array tasks up to 12.
  - Median run time so far is 3 s, the longest 365 s. Expected total: 1–3 h, mostly CoRAL and
    the bulk pool runs on the S=4 cases.
- **Not done yet:** scoring and collecting over all 96 runs, the results write-up, and
  Yubi's decisions on CoRAL / CycleExtractor (below).
- Nothing pushed. Yubi pushes.

**Rules.** Small local runs are fine. New Slurm arrays, sweeps or long limits need Yubi's
approval (the licence allows 2 Gurobi sessions, so never run local Gurobi jobs while the array
holds both). Wait at least 60 s between `squeue` / `sacct` calls. Do not touch `~/RCK` or
`~/amplicon-allele-analysis` (`yubi_dev`).

**Next steps, in order:**
1. **Check the array.**
   ```bash
   cd ~/RCK-cnp/results/bench_IW
   export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
   squeue -u $USER                                   # anything left of 32267060?
   sacct -j 32267060 --format=JobID,State,Elapsed,MaxRSS -X | sort -k2 | uniq -c -f1 | head   # states
   python bench_iw.py status                         # done / missing / failed per (case, method)
   ls /n/fs/yubi-iw/runs/bench_IW/logs/slurm/        # bench-iw_32267060_<task>.out per task
   ```
   - A run is `done` when its outputs exist; `failed` when its `run.json` says so.
   - Look at the task log, and at `<run>/run.json` (command, attempts, elapsed). CoRAL / CE also
     write `<method>.stdout` in the run folder.
   - **If some failed or never ran:** diagnose first (licence error? time limit? a bug?). Only
     then resubmit with `sbatch --array=<idx,...>%2 --output=... bench_iw.sbatch`; finished runs
     are skipped, so the full array can also be resubmitted. Ask Yubi before resubmitting more
     than a handful.
   - `sacct` MaxRSS: if every task is far below 4G, note it for later arrays (`--mem`).
2. **Score and collect (offline, no Gurobi, local).**
   ```bash
   python bench_iw.py score --jobs 8    # export + aa_eval + bench_metrics per finished run
   python bench_iw.py collect           # bench_iw_table.tsv, bench_iw_pooled.tsv, bench_iw_proportions.tsv here
   ```
   - `score` rescores every finished run and overwrites `<run>/score/`. Runs that fail to score
     leave `<run>/score/score_error.txt`, a traceback; `collect` marks them `score_status=error`.
   - The three TSVs are written into `results/bench_IW/`. Commit them with a short README of the
     results (step 3).
3. **Report to Yubi:**
   - the per-(method, class, S) pooled table;
   - failures and rejections;
   - the CoRAL / CycleExtractor behaviour below, across all 24 cases.

   Do not change the metrics or methods before Yubi decides.
4. **Decisions waiting for the full results (Yubi):** *(resolved 2026-11-14 night: linear walks are
   reported both kept and dropped; see "Decisions that apply from here", 8)*
   - CoRAL and CycleExtractor answers include linear walks (`0+ ... 0-`) through non-ecDNA
     segments. They count as ecDNA mass, which lowers copy-number precision and proportion
     scores. Should precision (or everything) drop linear walks, as `score_h2h.normalise_cycles(...,
     drop_paths=True)` can? Recall is unaffected.
   - Whether the bulk route should also be scored with its `lambda_0` answer (not run here; only
     the pool run with lambda 0.2 re-ranking is).

**Pilot (`simple_S1_s11001`, one species, 4.97 copies; run locally before submitting):**

| method | time | cn_recall | cn_precision | proportion_tv (strict / all) | missed / false |
|---|---|---|---|---|---|
| rck_bulk | 1 s | 1.00 | 0.93 | 0 / 0 | 0 / 0 |
| rck_mc | 5 s | 1.00 | 0.96 | 0 / 0 | 0 / 0 |
| coral | 58 s | 1.00 | 0.32 | 0.11 / 0.22 | 0 / 1 |
| ce | 1 s | 1.00 | 0.16 | 0.58 / 0.58 | 0 / 2 |

What the answer files show (no further tests were run):
- **CoRAL** finds the species as a closed cycle (3.96 copies). It then adds three low-copy linear
  walks (0.50, 0.35, 0.27) over segments the species does not use, including the 10 kb flanks
  `graph_route.run_coral` adds.
  - `max_weight` mode keeps decomposing the remaining copy number, which is mostly chromosomal
    background: the exported graph's copy numbers are total copy number.
  - The long-read path constraints are credited to one of the linear walks, not to the cycle. Why
    is not clear from the files.
- **CycleExtractor** returns no closed cycle. It returns three linear walks (2.73, 2.57, 1.21):
  - two carry ecdna1's segments, each repeating part of it;
  - one is a single non-ecDNA segment.
  - Their sum, 6.51, is close to ecDNA plus about 2 copies of background. So it seems to spread
    the ecDNA copies over background-carrying walks (an inference, not checked against its code).

## What bench_IW is

- **Classes.**
  - `simple_S{1,2,4}`: species share no reference base. `dup` is kept, at most 2 traversals of a
    segment.
  - `hard_S{1,2,4}`: S >= 2 nests every species inside the first one's intervals (outer
    0.6–3 Mb). S = 1 is one ctrl_dup-like species with internal repeats or a fold-back; at least
    one segment traversed twice, at most 3.
  - 4 accepted seeds per scenario: 24 runs.
- **Shared ground truth.** Every run has `clone_seed: 7` and demo7's realised clone proportions
  `[0.2, 0.49981378091514816, 0.3001862190848519]`. These reproduce demo7's grid, phasing, clone
  tree, clone CNPs and proportions exactly. Also shared:
  - `ecdna_positive_clones` null (both tumour clones), `ecdna_clone_weights` [1, 0.5];
  - long reads, library 5e6 (about 6.5 junction reads per copy);
  - `hap_labels: inherit`, `p_cross` 0.2.
- **Drawn per seed** (`make_scenarios.iw_draw`, from `default_rng(random_seed)`, written into
  each config and its `bench_iw:` record block):
  - **size per species:** lognormal, median 1 Mb, sigma 1.3, within 0.3–3 Mb; simple sets
    redrawn until they sum to <= 8 Mb;
  - **hard nesting:** `nest_size_fraction` U(0.3, 0.7) per nested species;
  - **copies per species:** lognormal, median 8, sigma 1.2, within 4–60; sets redrawn until
    max/min <= 8 (`abundance_mode: per_sample`, `ecdna_burden` = sum, `ecdna_proportions` =
    shares);
  - **`num_intervals`:** simple U{1,2}, hard U{2,3};
  - **breakpoints per species:** U{2 .. min(8 or 10, ceil(size / 150 kb))}.

  The bounds are AmpliconRepository statistics (2,963 human ecDNA features); the reasons are in
  `make_scenarios.py`.
- **Seed checks** (`bench_nested.check_seed`, intent from `make_scenarios.IW_INTENT`):
  - bulk ratio within 10%;
  - every ecDNA junction has >= 20 reads;
  - traversal cap;
  - footprint <= 8 Mb;
  - every species placed with all its intervals;
  - separate / nested (>= 80% inside) / dup.

  Accepted 24 of 30 candidates. Rejections are in `configs/bench_IW/seeds.tsv`.
- **Subgraphs.** bench_nested's flags plus `--min-component-nodes 2 --min-sv-support 20`. The
  default (>= 5 breakpoint clusters) dropped a one-interval species in `simple_S4_s13004`. Every
  union subgraph has all truth segments and junctions (`subgraphs/bench_IW/subgraph_checks.tsv`).
- **Methods** (`bench_iw.py` docstring has the exact flags):
  - `rck_bulk`: pool run, `--block-flip unoriented`, answer re-ranked at lambda 0.2;
  - `rck_mc`: `--model multiclone` on the DIAG built with `--clone-profiles truth
    --diag-labeling-groups block`;
  - `coral`: CoRAL 3.1.1 on the exported union subgraph, with 10 kb flanks;
  - `ce`: CycleExtractor 1.0.0 on the same graph.
- **Metrics,** per (case, method), in `<run>/score/`:
  - **`aa_eval`:** `aa_species_matrix.tsv` (J_seg, J_adj; blind rows for every method, phased
    rows meaningful for `rck_bulk` only), M1 / M2 (`aa_partition_metrics.tsv`), strict and
    best-of species scores. No truth mixture graph was built, so there are no graph-level
    metrics.
  - **`bench_metrics`** (`rck/ecdna/bench_metrics.py`, approved definitions):
    - `cn_recall` = sum_x min(c, t) / sum_x t, over every base. Ecdna copies only,
      haplotype-blind, bp-weighted, strict cycles. t = sum_s PHI_s x traversals; c = sum of
      cycle bulk weight x traversals. `rck_mc` weights are sum_k u_k psi_k, with u from the DIAG.
    - `cn_precision` = the same over sum_x c.
    - `cn_within_tol`: the share of truth bases with |c - t| <= max(1, 0.1 t).
    - `_all` variants include low_weight cycles.
    - Proportions: truth PHI_s / sum vs inferred w / sum over the strict set and over all
      cycles, greedy blind-Jaccard matching. `proportion_tv` is the total-variation distance
      (unmatched mass counted); `proportion_mae_matched`.
- **Table columns** (`bench_iw_table.tsv`, one row per case x method):
  - case fields: `name`, `scenario`, `class`, `S`, `replicate`, `seed`;
  - run fields: `method`, `run_status`, `elapsed_s`;
  - the bench_metrics fields;
  - `n_inferred_strict`, `n_missed`, `n_false`, `mean_composite_h0`;
  - M1 / M2 headline precision / recall / F1 (`species_detection_*`, `adjacency_pairs_*`);
  - `{blind,phased}_jseg_best_mean`, `_jadj_best_mean`: the mean over truth species of the
    best strict cell; `_exact_cycles`.
- **`bench_iw_pooled.tsv`:** mean and median per (method, class, S), (method, class) and method;
  the `level` column says which.

## Where things are

| what | where |
|---|---|
| configs, seeds.tsv, shared_clones.tsv | amplicon-allele-analysis-cnp `simulations/configs/bench_IW/` (tracked) |
| simulator runs | `/n/fs/yubi-iw/runs/bench_IW/runs/` (linked as `simulations/runs/bench_IW`, ignored) |
| union subgraphs / multi-clone DIAGs | `simulations/subgraphs/bench_IW/`, `subgraphs/bench_IW_cnp/` (tracked) |
| simulation and subgraph logs | `/n/fs/yubi-iw/runs/bench_IW/logs/` |
| method runs and scores | `/n/fs/yubi-iw/runs/bench_IW/methods/<case>/<method>[/amplicon1]/{run.json,score/}` |
| exported graphs (CoRAL, CE input) | `/n/fs/yubi-iw/runs/bench_IW/methods/<case>/graph/` |
| Slurm logs | `/n/fs/yubi-iw/runs/bench_IW/logs/slurm/bench-iw_32267060_<task>.out` |
| driver, sbatch, results TSVs | RCK-cnp `results/bench_IW/` |
| **bench_IW_intpsi** (the benchmark from 2026-11-14 evening): configs, seeds.tsv | `simulations/configs/bench_IW_intpsi/` (tracked) |
| bench_IW_intpsi simulator runs | `/n/fs/yubi-iw/runs/bench_IW_intpsi/runs/` (linked as `simulations/runs/bench_IW_intpsi`, ignored) |
| bench_IW_intpsi subgraphs | `simulations/subgraphs/bench_IW_intpsi/`, `subgraphs/bench_IW_intpsi_cnp/` (tracked) |
| bench_IW_intpsi simulation and subgraph logs | `/n/fs/yubi-iw/runs/bench_IW_intpsi/logs/` |
| bench_IW_intpsi method runs (planned) | `/n/fs/yubi-iw/runs/bench_IW_intpsi/methods/` |

Commands on the simulation side (amplicon-allele-analysis-cnp `simulations/`; add `--set iw_intpsi`
as the first argument for bench_IW_intpsi):
- `python bench_IW.py seeds | report | check-shared`;
- `python bench_IW.py [--cnp] rebuild [RUN ...]`;
- `python bench_IW.py scenario <name> --accept 4 --candidates 20` (resumable; existing runs are
  reused).

**Session commits (2026-11-14), all on `yubi_CNP_ionic`, none pushed.**
- **amplicon-allele-analysis:**
  - `fead61a` `clone_seed`, `6283c0e` its test fix;
  - `0e05619` old bench configs and subgraphs removed;
  - `0c3a727` per-species structure settings;
  - `2c8f810` bench_IW code;
  - `21cce8c` bench_IW configs and subgraphs;
  - plus the docs commit that wrote this session's CHANGELOG / README entries.
- **RCK:**
  - `1b5bfd3` bench result data removed (scripts kept);
  - `eef70e6` bench_metrics, driver, sbatch;
  - plus the docs commit that wrote this file.
- **Local only:** in amplicon-allele-analysis-cnp, the `third_party/ecSimulator` (pinned
  `07a0e3c`) and `third_party/cnplot` (pinned `1e8bfef`) submodules were initialised. Without
  ecSimulator the simulator writes no `truth/aa/`, and the seed checker's traversal count
  reads `truth/aa/*_cycles.txt`.
