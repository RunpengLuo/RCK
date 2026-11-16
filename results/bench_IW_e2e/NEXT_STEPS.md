# bench_IW_e2e: status and next steps (branch `yubi_CNP_ionic`)

*Written 2026-11-15, 22:15, at the end of the session that built it. Plan of record:
[`results/bench_IW/NEXT_STEPS.md`](../bench_IW/NEXT_STEPS.md), "Addition 1 in detail" (reviewer item 1:
CoRAL and CE end to end). Folder layout and methods: [README.md](README.md). Read simulation and the
partial-genome design: amplicon-allele-analysis `simulations/e2e/README.md`.*

## Resume here

1. **Check the pilot array, Slurm job `32318772`** (tasks 12 = `hard_S1_s14001`, 21 = `hard_S4_s16003`;
   submitted 2026-11-15 22:11, 16 CPUs, 40 GB, 6 h each; `pilot_e2e.sbatch`).
   - Log: `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/logs/slurm/iwe-pilot_32318772_{12,21}.out`.
     Every step prints `== rc <code>, <seconds> s: <command>`, which gives the timings for sizing the arrays.
   - `sacct -j 32318772 --format=JobID,State,Elapsed,MaxRSS` and `seff 32318772_12`, `seff 32318772_21`.
   - `python bench_iw_e2e.py status` (gurobi-env python) shows reads / inputs / each method per case.
   - The chain per task: reads (`emit_reads.py` + `align_cnv.sh`) → `check_reads.py` → inputs
     (`obs_from_bam.py` + `build_subgraphs_e2e.py`, union subgraph and DIAG) → `coral_e2e` → `ce_e2e` →
     `rck_bulk_e2e` → `rck_mc_e2e` → `score --case IDX`. A failed method does not stop the others.
   - hard_S1 resumes at CNVkit (its reads and BAMs with sentinels were finished locally); hard_S4 starts over.
   - If a task failed: the step's own log is in the case folder (`reads/<case>/cnvkit/batch.log`,
     `subgraphs*/<case>/subgraph_construct.log`, `methods/<case>/<method>/run.json` and `*.stdout`).
     Fix, then resubmit the same command; finished steps are skipped.
2. **Check the pilots' outputs** (the questions the pilots were run to answer):
   - **Reads against the counts:** `reads/<case>/check_reads.json`. Expect depth / expectation ≈ 0.97–1.0
     (background and ecDNA, tumour and normal), ecDNA ALT-fraction correlation ≈ 0.998. Junction-spanning
     reads are expected at ≈ 2.4× the count layer's `#NOV` (see "Findings").
   - **Sentinels restore the purity / ploidy fit:** `subgraphs/<case>/union/scale.json`. Expect
     `ploidy_class` 2, purity ≈ 0.8, `depth_per_copy` ≈ 14.8 (the count-layer route's fit on hard_S1 was
     class 2, purity 0.80, 14.87; footprint-only reads gave class 4, purity 0.28, 10.57). If it is still
     wrong, the sentinels are not enough: discuss with Yubi (bigger sentinels, or whole-genome reads for
     the background) before going on. Also note `junction_reads_per_copy` and `junction_reads_method`
     there and in `methods/<case>/rck_bulk_e2e/summary.json` (θ is now estimated, not set; below).
   - **RCK-ecDNA end to end runs at all:** `methods/<case>/rck_{bulk,mc}_e2e/run.json` status `ok`; this is
     the first time `rck_bulk_e2e` / `rck_mc_e2e` run (the BAM-derived subgraph had not been built with
     sentinels before the session ended).
   - **Scores:** `python bench_iw_e2e.py collect`, then compare with the earlier rows of the same cases in
     `results/bench_IW/bench_iw_intpsi_table.tsv` (the snippet used on 2026-11-15 is under "Pilot
     results so far").
3. **Resize the arrays from the pilot timings, then ask Yubi before submitting** (`submit_e2e.sh` prints
   the commands; `--submit` submits). The values in it (reads 1 h / 16 GB, coral 3 h / 8 GB, ce 2 h / 4 GB)
   were set before the sentinel windows, which made the reads step about 5× longer; they are too small.
   The script also has no arrays yet for `inputs`, `rck_bulk_e2e` and `rck_mc_e2e`: add them
   (`inputs` after reads, no Gurobi; `rck_*_e2e` after inputs, `%1` each for the 2-session licence; the
   `rck_bulk_e2e` limits as `results/bench_IW/submit_intpsi.sh`'s `rck_bulk`, 90 min). Run
   `python bench_iw_e2e.py mixture` once (truth mixture graphs for all 24 cases, no Gurobi; the two
   pilots have theirs) before scoring.
4. After the arrays: `score`, `collect`, then a write-up (`RESULTS_e2e.md` here, figures from
   `results/bench_IW/make_figures.py` pointed at `bench_iw_e2e_table.tsv`, as new figures; the existing
   figures stay as they are).

## Decisions so far (with whom)

- **Route (B): everything from the same BAMs** (Yubi: "proceed with the steps in the outline"; outline
  recommended B). CoRAL, CE and RCK-ecDNA all read the simulated reads; RCK-ecDNA's `obs/` is rebuilt from
  the BAMs (`obs_from_bam.py`).
- **Keep previous results separate** (Yubi). Nothing in `results/bench_IW/`, the bench_IW_intpsi simulator
  runs, subgraphs or method runs is written; see README "Kept apart from the earlier results".
- **Partial-genome reads, kept for now** (Yubi, 2026-11-15): each ecDNA interval ± 5 Mb plus a 2 Mb
  sentinel per chromosome arm (median 102 Mb per case). Reasons, sources, limitations and the open
  whole-genome check: amplicon-allele-analysis `simulations/e2e/README.md`, "Partial-genome design".
- **`--rck-junction-reads-per-copy 6.5` dropped** (Yubi, 2026-11-15). bench_nested passes the count
  layer's θ to `subgraph_construct.py`, which only writes it into `scale.json`. It matters only for the bulk
  route (junction copy number = `#NOV` / θ, `rck/ecdna/scale.py`); the multi-clone route reads θ only with
  `--mc-junction-weight` > 0 (default 0, never set by the bench). The reads carry ≈ 2.4× the count layer's
  junction reads, so 6.5 would be wrong; without the flag, `rck-ecdna` estimates θ itself
  (`scale.estimate_junction_scale`). Note for the write-up: the earlier `rck_bulk` was given the true θ,
  `rck_bulk_e2e` is not.
- **`rck_mc_e2e` keeps the true clone profiles** (`--clone-profiles truth`), as `rck_mc` did; its segment
  copy numbers come from the truth, its graph (breakpoints, adjacencies, `#NOV`, haplotype labels) from the
  BAMs. Estimated profiles are reviewer item 2, not this one.
- **CoRAL and CE at their defaults**, except CoRAL's `--solver-time-limit 1800` (as the earlier runs). CE runs
  on each CoRAL amplicon graph; both report one cycles file per amplicon, pooled into one answer per case.

## Choices made without a decision from Yubi (check these)

- **Own split-read clustering instead of an SV caller** in `obs_from_bam.py`: it gives `#NOV`, per-haplotype
  counts per end and read paths in the exact form `subgraph_construct.py` reads. A reviewer may ask why not
  Sniffles (installed in `e2e-align`, version 2.8.1, unused).
- **SNPs**: het SNVs called by bcftools on the normal (QUAL ≥ 20, DP ≥ 8), phased by WhatsHap on both BAMs.
- **Haplotype per read alignment**: ≥ 2 phased SNPs, ≥ 75% agreeing; an adjacency read counts as A or B
  only when its known ends agree.
- **NanoSim**: ecSimulator's PromethION r10.4.1 model (trained with NanoSim 3.1.0) in NanoSim 3.2.3; no
  base-quality profile, so reads get a constant Phred 20; read counts raised by the model's unaligned-read
  rate (14%), which otherwise left depth 15% low.
- **CNVkit in its own env** (0.9.14): the shared env resolved 0.9.8, which fails on pandas 2.

## Findings so far (2026-11-15 pilots, before the sentinel windows)

- **Reads match the counts** except junction reads: depth / expectation 0.90–1.01; ecDNA ALT fraction
  correlation 0.997–0.998 (MAE 0.02). Junction-spanning reads are ≈ 2.4× the count layer's `#NOV`
  (hard_S1: 1,222 vs 466 over 6 ecDNA adjacencies; hard_S4: 3,494 vs 1,539 over 19). Cause: the count
  layer's `gamma_frag` = 0.0005 fragments per base per copy, about half the read-start rate its own γ
  implies (γ / mean read length ≈ 14.8 / 15,000 ≈ 0.001). This is a property of the simulator config,
  worth a line in the simulator's docs; the reads are consistent with the depth.
- **CoRAL end to end ran**: seeds cover every truth interval; reconstruct 90–105 s, optimal.
- **Pilot results so far** (strict, `mean_composite_h0`; earlier rows from `bench_iw_intpsi_table.tsv`):

  | case | coral_e2e | ce_e2e | coral (RCK bulk's graph) | ce (same) | rck_bulk (counts) | rck_mc (counts) |
  |---|---|---|---|---|---|---|
  | hard_S1_s14001 | 0.89 | 0.89 | 0.96 | 0.96 | 1.00 | 0.70 |
  | hard_S4_s16003 | 0.28 | 0.28 | 0.05 | 0.17 | 0.48 | 0.73 |

  CoRAL's own graph on hard_S4 is nearly right (breakpoint accuracy 1.0, novel-adjacency FDR 0.05); its
  copy count on hard_S1 is 6.6 against a true 11.7 copies per cell (cn_recall 0.56). These were before the
  sentinel windows and will be replaced by the pilot array's numbers.
  Comparison snippet (run in this folder, gurobi-env python):

  ```python
  import pandas as pd
  cols = ["name", "method_variant", "n_inferred_strict", "n_missed", "n_false", "mean_composite_h0",
          "cn_recall", "cn_precision", "blind_jseg_best_mean", "blind_jadj_best_mean", "blind_exact_cycles"]
  e = pd.read_csv("bench_iw_e2e_table.tsv", sep="\t"); e = e[e["run_status"].notna()]
  o = pd.read_csv("../bench_IW/bench_iw_intpsi_table.tsv", sep="\t"); o = o[o["name"].isin(e["name"])]
  print(pd.concat([o[cols], e[cols]]).sort_values(["name", "method_variant"]).round(3).to_string(index=False))
  ```

## Open items

- **Whole-genome check** of the partial-genome design on the two pilot cases (amplicon-allele-analysis
  `simulations/e2e/README.md`, "Open"); about 3–4 h and tens of GB per case. Yubi to decide when.
- **Sniffles instead of the split-read clustering**, if Yubi wants a standard SV caller on RCK-ecDNA's side.
- **Optional short-read arm** (CE on AmpliconArchitect graphs; outline step 8): not started.
- **Nothing is committed.** New, uncommitted: RCK `results/bench_IW_e2e/` (this folder);
  amplicon-allele-analysis `simulations/e2e/`. Also uncommitted from the same day: RCK
  `docs/pipeline/evaluation.md` (audit) and `results/bench_IW/NEXT_STEPS.md` (the outline). Commit when
  Yubi asks.

## Where things are

| what | where |
|---|---|
| driver, sbatch, submit script, tables | RCK `results/bench_IW_e2e/` |
| read simulation, BAM → obs, subgraph build | amplicon-allele-analysis `simulations/e2e/` |
| reads, BAMs, CNVkit (`<case>/`) | `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/reads/` |
| BAM-derived obs runs (`obs/` + `truth` link) | `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/obs_runs/` |
| union subgraphs / DIAGs from the BAMs | `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/subgraphs/`, `subgraphs_cnp/` |
| truth mixture graphs | `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/mixture/` |
| method runs and scores | `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/methods/<case>/<method>/` |
| logs (local pilots, Slurm) | `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/logs/`, `logs/slurm/` (`submitted.txt`) |
| GRCh38 no-alt FASTA + .fai, minimap2 map-ont index, NanoSim model | `/n/fs/yubi-iw/ref/` |
| envs | `/n/fs/yubi-iw/conda/envs/e2e-align`, `e2e-cnvkit`, `e2e-nanosim` (plus the existing `gurobi-env`, `coral-3.1.1`) |
