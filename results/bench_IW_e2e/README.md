# bench_IW_e2e: CoRAL and CycleExtractor end to end on bench_IW_intpsi

*Started 2026-11-15 (Yubi: "proceed with the steps in the outline"; plan in
[`results/bench_IW/NEXT_STEPS.md`](../bench_IW/NEXT_STEPS.md), "Addition 1 in detail").*

The finished bench_IW_intpsi comparison in [`results/bench_IW/`](../bench_IW/) gives CoRAL and CE
RCK bulk's graph to decompose ("same graph", decomposition only). Here each baseline builds its own
graph from simulated long reads of the same 24 cases and is scored against the same truth with the
same scorer.

## Kept apart from the earlier results

- **Nothing in `results/bench_IW/` is written or edited.** `bench_iw_e2e.py` imports `bench_iw.py` (case
  list, scorer, table code) and `results/graph_route/graph_route.py` (CoRAL / CE settings) and points
  copies of their paths here; the TSVs it writes are `bench_iw_e2e_*.tsv` in this folder.
- **Data** live under `/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e/` (`reads/`, `mixture/`, `methods/`, `logs/`).
  The simulator runs in `/n/fs/yubi-iw/runs/bench_IW_intpsi/runs/` and the method runs in
  `/n/fs/yubi-iw/runs/bench_IW_intpsi/methods/` are only read.
- **Read simulation code** is new, in amplicon-allele-analysis `simulations/e2e/` (its README); the
  simulator's `src/` is not edited.

## Steps

| Step | Command | Where |
|---|---|---|
| reads, BAMs, CNVkit | `bench_iw_e2e.py reads IDX` (= `e2e/emit_reads.py` + `e2e/align_cnv.sh`) | `reads/<case>/` |
| read check (pilots) | `e2e/check_reads.py` (e2e-align python) | `reads/<case>/check_reads.json` |
| CoRAL end to end | `bench_iw_e2e.py run --method coral_e2e IDX` | `methods/<case>/coral_e2e/` |
| CE end to end | `bench_iw_e2e.py run --method ce_e2e IDX` (on CoRAL's graphs) | `methods/<case>/ce_e2e/amplicon<N>/` |
| truth mixture graphs | `bench_iw_e2e.py mixture` (`truth_mixture.py --outdir`) | `mixture/<case>/` |
| scoring | `bench_iw_e2e.py score` | `methods/<case>/<method>/score[_nolinear]/` |
| table | `bench_iw_e2e.py collect` | `bench_iw_e2e_{table,pooled,proportions}.tsv` here |
| arrays | `submit_e2e.sh` (prints; `--submit` submits) | `logs/slurm/` |

## Reads cover part of the genome

Reads are simulated for each ecDNA interval ± 5 Mb plus a 2 Mb window at the middle of every
chromosome arm (median 102 Mb per case, about 3.7% of the genome), because reads for the whole
simulated genome would be 150 Gbp per case. The design, its sources (ecSimulator's background
regions, CoRAL's 10 Mb baseline requirement, the pilot purity / ploidy fit) and its limitations are in
amplicon-allele-analysis
[`simulations/e2e/README.md`, "Partial-genome design"](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_CNP_ionic/simulations/e2e/README.md#partial-genome-design).
Every method reads the same BAMs; report the results as regional WGS.

## Methods (fair-use record)

- **coral_e2e**: CoRAL 3.1.1 `seed` at its defaults (`--gain 6`, `--min-seed-size 100000`,
  `--max-seg-gap 300000`, GRCh38 centromeres) on CNVkit's `.cns`, then `reconstruct` at its defaults
  except `--solver-time-limit 1800` and `--solver-threads` (as the earlier CoRAL runs).
- **ce_e2e**: CycleExtractor 1.0.0 at its defaults on each CoRAL amplicon graph (path constraints included).
- Both report one cycles file per amplicon; they are pooled into one answer per case for scoring.

Status, decisions and open items: [NEXT_STEPS.md](NEXT_STEPS.md).
