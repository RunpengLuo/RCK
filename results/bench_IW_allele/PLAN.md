# bench_IW_allele: analyses that show RCK-ecDNA's allele-specific reconstruction (branch `yubi_CNP_ionic`)

*Started 2026-11-16 (Yubi: "write down your ideas in a document in the repo that can be referred to in
future sessions. then, start items 2 and 3 ... make sure to describe how [item 1] can be implemented").
Revised 2026-11-16 (later the same day) after reading the RCK paper's evaluation: the plan now follows
RCK's own evaluation format first (Part A, Part B), and item 1 became the **allelic-sibling** benchmark
instead of allelic twins (reasons under "Alignment with the literature"). The end-to-end comparison
([`../bench_IW_e2e/NEXT_STEPS.md`](../bench_IW_e2e/NEXT_STEPS.md)) is paused by Yubi on 2026-11-16,
blocked on the purity / ploidy fit. Revised again 2026-11-16 (Yubi's decisions): Part A + Part B mimic RCK's
2020 evaluation as closely as possible, on the multi-clone route with the true clone profiles and
proportions as input (RCK's simulation also used the true copy numbers), P = 0.9 only, no homogeneous S′
cases, and the bulk route kept only as the graph CoRAL / CE decompose.*

## Resume here

**State at hand-off (2026-11-16, ~21:20 EST).** Part B's 72-task array **32364902** (`0-71%1`, submitted
21:10 EST with Yubi's go-ahead) was running: 11 of 72 runs done, all optimal (the 3 pilot runs, array
32364762, included), about one task per minute, expected to finish ~22:20–22:30 EST. Everything else from the
session is committed on `yubi_CNP_ionic` (RCK `ec5284a` … this commit; amplicon-allele-analysis unchanged),
not pushed. Partial Part B tables / Figure 10 were deliberately not committed; step 4 below regenerates them.

**Pick up here, in `results/bench_IW_allele/`** (`/n/fs/yubi-iw/conda/envs/gurobi-env` python):

1. **Check the array.** `python rck_partB.py status` reads the run folders (no Slurm calls): every variant
   should show `'done (optimal)': 24`. If needed, `sacct -j 32364902 --format=JobID,State,Elapsed,MaxRSS`
   (≥ 60 s between `squeue` / `sacct` calls). Logs: `/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/logs/slurm/
   iwa-partB_32364902_<task>.out` (task index → variant / case: `python rck_partB.py list`).
2. **If any task failed or is `done (status ≠ 2)`:** read its log and `<case>/<variant>/rck_mc/run.json`,
   diagnose before anything else (licence waits are retried 15 × 2 min inside a task; a Slurm TIMEOUT at
   30 min means the MILP needed longer than the 600 s limit allows end to end). Resubmitting
   (`./submit_partB.sh --submit`; finished runs are skipped) is a new Gurobi submission: **ask Yubi first**.
3. **Export the answers (no Gurobi, seconds):** `python rck_partB.py export --jobs 8` → `<case>/<variant>/
   rck_mc/score/answer/<case>_rck_mc_cycles.txt`.
4. **Score, collect, plot (local, seconds):** `python score_partB.py run --jobs 16 && python score_partB.py
   collect && python make_partB_figures.py` → `partB_table.tsv`, `figures/fig10_partB.{pdf,png}`,
   `figures/table4.{md,tsv}`, `figures/figure_data/fig10_partB.tsv`. Finished rows are cached in
   `<run>/score/partB.json` (the labelled reference's in `<case>/labelled/score/partB.json`); delete those to
   rescore after a scorer change. 30 rows (24 labelled + 6 Part B) were already scored at hand-off.
5. **Look before writing:** open Figure 10 (layout and values); check every run optimal (Table 4 "runs
   (optimal)"); RCK's FPR over spurious inputs (Table 4 footer); and the pilot observation below.
6. **Write "Part B results"** under "Part B status" (RCK's Fig. 2 reading: what CHIL alone, averaging, and
   the full errors each cost against the labelled reference, single vs mixed cases; descriptive, n = 6
   mixed), update this list, commit `partB_table.tsv`, Figure 10 and Table 4 with the text.
7. **Then, waiting for Yubi:** item 1 (allelic siblings): the simulator option first (local, cheap; amp repo
   `simulations/src/structures.py`), method runs need arrays and Yubi's go-ahead.

**Observation to check across all cases (from the pilot case only):** on hard_S4_s16003 (the known weak
multi-clone case, item 8) the label-free `chil` input gave near-exact clone copy numbers (Eq. 11 over clones,
ecDNA bases, relative ~0 vs 0.25 labelled; allele CN recall 0.96 vs 0.71) but worse junction phasing
(labelled adjacency FNR 0.37 vs 0.26, FDR 0.25 vs 0.00). The hard junction labels may push that case to a
wrong copy-number solution. If it holds on other cases it belongs in the results (and in item 8).

**Done before this point:** items 2 and 3 (figures 7, 8, Table 2), Part A (Figure 9, Table 3), Part B steps
1–2 and the pilot (sections below). To redo any of them from scratch, see the Files table; every driver
skips finished tasks.

Rules for this line of work (Yubi's, carried over): keep everything in `results/bench_IW_allele/` and
separate from `results/bench_IW/` and the bench_IW_intpsi runs (only read); new Slurm arrays and new Gurobi
work need Yubi's go-ahead; Gurobi licence 2 sessions (`%1` per Gurobi array, no local Gurobi while arrays
hold both); scripts that would take more than a few minutes locally go to Slurm on ionic instead; wait ≥ 60 s between `squeue` / `sacct`
calls; don't touch `~/RCK` or `~/amplicon-allele-analysis` (`yubi_dev`); don't patch CoRAL or CE; commit on
`yubi_CNP_ionic`, don't push.

## Why these analyses

The current manuscript results (`results/bench_IW/RESULTS.md`, figures 2–6 and Table 1) mostly show
**structure and copy-number** accuracy, where RCK-ecDNA wins clearly. The allele-specific claim rests
on Figure 6 alone, which has no baseline, and on it phased and allele-blind scores agree in 37 of 48 runs.
A reviewer can read that as: "in this benchmark, getting the alleles right comes for free once the
structure is right; a haplotype-blind method could assign alleles afterwards." The analyses below test
that directly, and they are laid out so that a reviewer who knows RCK recognises the evaluation.

## Alignment with the literature

**How RCK evaluated haplotype-specific reconstruction** (Aganezov & Raphael 2020, read in full on
2026-11-16 from Europe PMC, PMC7545144; section names are the paper's):

- *Evaluation of RCK on simulated data* / *Methods: Simulating rearranged cancer samples*: 100 two-clone
  samples (~600 shared and 100 / 200 clone-specific rearrangements, simple and complex). The input
  hides the truth by **clone and haplotype information loss (CHIL)**: every segment's allele-specific copy
  numbers are flipped A/B with probability 0.5, and novel adjacencies lose their A/B labels. A second input
  adds **errors**: ±50 bp on breakpoints, 10% spurious novel adjacencies, copy numbers averaged over 50 kb
  fragments, and ±1 on 5% of fragments.
- Scoring is **haplotype-specific "up to haplotype symmetries"**: false negative and false positive
  rates of novel adjacencies (Fig. 2A), and a **length-weighted allele-specific segment copy-number
  distance** (their Eq. 11: per fragment, the minimum over the A/B swap of
  Σ_segments Σ_clones (|a − d| + |b − d′|) · length; Fig. 2B).
- *Evaluation of RCK on prostate cancer* (17 samples): compared with ReMixT (McPherson et al. 2017), which
  is clone- and allele-specific for segments but "does not infer haplotype A/B labels for the extremities
  that are involved in each novel adjacency". The haplotype argument there is **plausibility**: ReMixT's
  karyotypes violate the **generalized infinite-sites (IS) constraints** (Fig. 4): a novel adjacency that
  must be labelled both A and B ("at least two large-scale somatic rearrangements occurred independently
  at the same genomic position on both homologous chromosomes, which is highly unlikely"), the same across
  clones, and reciprocal extremities used on opposite homologs. RCK's karyotypes have none by construction.
- There is **no haplotype-blind or post-hoc-phasing baseline and no ablation** in RCK; its claim is
  carried by the simulation scores and the IS-violation comparison.

**What we do, in that frame:**

| RCK (2020) | here | where |
|---|---|---|
| simulated two-clone samples, known haplotype-labelled genomes | bench_IW_intpsi: 24 cases, 2 tumour clones, haplotype-labelled ecDNA walks | done (bench_IW) |
| true clone- and allele-specific copy numbers as input (not HATCHet: that is the real-data section) | the subgraph labelled with the true clone profiles and proportions (`--clone-profiles truth`), multi-clone route | done (bench_IW) |
| CHIL on segments (A/B flipped per segment) | today one unknown A/B per phase block; per-segment flips in `chil` | **Part B** |
| CHIL on adjacencies (no labels) | our junction ends carry long-read haplotype labels; removed in `chil` | **Part B** |
| errors: 50 kb averaging, ±1 on 5%, ±50 bp, 10% spurious adjacencies | read-sampling noise only; RCK's errors in `chil_avg` / `chil_err` | **Part B** |
| P = 0.9 | `--run-nas-fp 0.1` (the default) | Part B |
| Eq. 11 summed over clones; copy-number state agreement | bulk-level Eq. 11 (Part A); per clone in Part B | **Part B** |
| scoring up to haplotype symmetries | within-phase-block scores (one A/B swap per phase block) | done (`_block` columns) |
| novel-adjacency FNR / FPR, haplotype-specific | phased J_adj only | **Part A** |
| Eq. 11 allele-specific CN distance | allele CN recall / MAE (items 2-3) | **Part A** (Eq. 11 itself) |
| generalized-IS violations, vs a method without adjacency labels | CoRAL / CE phased post-hoc, RCK-ecDNA, and the truth as a control | **Part A** |

**Beyond RCK** (our own designs; reasonable but without a published precedent to cite): the post-hoc
phasing baseline (item 2), allele copy numbers per segment (item 3), and the allelic-sibling benchmark
(item 1). The sibling benchmark has biological grounding in allele-specific copy-number work:
**mirrored-subclonal allelic imbalance**, different cells amplifying the two alleles of the same region,
is reported by CHISEL (Zaccaria & Raphael 2021, mirrored-subclonal CNAs in breast cancer single cells) and
Alleloscope (Wu et al. 2021, mirrored subclones in gastric, colorectal and breast cancer); the term goes
back to TRACERx (Jamal-Hanjani et al. 2017; the claim there was not re-read in this session). A literature
search on 2026-11-16 found **no report of ecDNA species derived from both parental alleles of one
locus**, so the sibling scenario is presented as a modelled case motivated by mirrored allelic imbalance,
not as an observed ecDNA phenomenon, unless a citation turns up.

**Why siblings and not twins** (the earlier item 1 design, superseded): a twin copies the source
species' junctions onto the other allele, i.e. the same novel adjacency on both homologs. That is exactly
what RCK's generalized IS forbids as "highly unlikely", and `rck/ecdna/ilp_circular.py` (`one-labeling`
constraint, one labelling per novel adjacency, inherited from RCK) cannot output it, nor can the
multi-clone route (RCK's `OptModelMultiClone`). Also, the subgraph labels a junction read from both AA and
BB copies `mixed`, which `subgraph_io._junction_phasings` turns into AB/BA only. A twin benchmark would
measure the model's assumption, not its accuracy, and contradicts the model it extends. Siblings (same
locus, other allele, **independent breakpoints**) respect the constraints.

Sources (DOIs checked against Crossref on 2026-11-16; Aganezov 2020 also read in full):
- Aganezov S, Raphael BJ. Reconstruction of clone- and haplotype-specific cancer genome karyotypes from
  bulk tumor samples. *Genome Res* 30:1274–1290 (2020). doi:10.1101/gr.256701.119
- McPherson AW, et al. ReMixT: clone-specific genomic structure estimation in cancer. *Genome Biol* 18 (2017).
  doi:10.1186/s13059-017-1267-2
- Zaccaria S, Raphael BJ. Characterizing allele- and haplotype-specific copy numbers in single cells with
  CHISEL. *Nat Biotechnol* 39:207–214 (2021). doi:10.1038/s41587-020-0661-6
- Wu CY, et al. Integrative single-cell analysis of allele-specific copy number alterations and chromatin
  accessibility in cancer (Alleloscope). *Nat Biotechnol* 39:1259–1269 (2021). doi:10.1038/s41587-021-00911-w
- Jamal-Hanjani M, et al. Tracking the evolution of non-small-cell lung cancer. *N Engl J Med*
  376:2109–2121 (2017). doi:10.1056/NEJMoa1616288
- Context already in `docs/pipeline/references.md`: CoRAL (Zhu et al. 2024), Decoil (Giurgiu et al. 2024;
  co-occurring ecDNA with overlapping footprints), HATCHet (Zaccaria & Raphael 2020).

## The analyses, ranked

| # | analysis | status | cost |
|---|---|---|---|
| A | **RCK-style evaluation of the existing runs**: haplotype-specific adjacency FNR / FDR / FPR, Eq. 11 CN distance, generalized-IS violations | **done 2026-11-16** | local, seconds |
| B | **RCK's simulation inputs on the multi-clone route**: CHIL, 50 kb averaging, full errors; P = 0.9; all 24 cases | steps 1–2 + pilot done; array 32364902 submitted 2026-11-16; scoring scripts ready; results pending ("Resume here") | 72 `rck_mc` runs (array, `%1`) |
| 1 | **allelic-sibling cases**: a second species on the other allele of the same locus, own breakpoints | designed (below), not started | simulator option + seeds; method arrays |
| 2 | **post-hoc phasing baseline**: CoRAL / CE (and RCK's own structures) phased afterwards by a rule | **done 2026-11-16** | local, minutes |
| 3 | **allele-specific copy numbers** per segment: true vs inferred CN_A / CN_B, BAF error | **done 2026-11-16** | local, minutes |
| 4 | ablation: RCK with the alleles collapsed (allele-blind graph) | idea | small array |
| 5 | clone-level results of the multi-clone route (ψ per clone; Eq. 11 per clone), and estimated clone profiles | idea (bench_IW NEXT_STEPS item 2) | local scoring; profile perturbation runs |
| 6 | when phasing fails: phased − blind scores vs SNP density / phase-block length | idea | local on existing runs; a SNP-density sweep needs seeds |
| 7 | one worked example figure: truth vs each method, segments coloured by allele | idea | local |
| 8 | explain or fix the multi-clone weak spots (hard_S4_s16003; extra split cycles) | open (bench_IW NEXT_STEPS item 6) | local |

Order: A (in RCK's own terms, on what exists), then B (completes RCK's simulation design on the
multi-clone route), then 1 (the case RCK's evaluation never isolates: both alleles of a locus amplified).

Items 4–8 in a sentence each:
- **4.** Run `rck-ecdna` on the 24 cases with A and B merged in the input (or the model's allele split
  switched off, if it has a switch); the drop on nested S = 4 cases shows what the allele model buys.
- **5.** Score ψ_{i,s} (copies of species s per cell of clone i) against the truth as a scatter or heatmap
  (per-clone Eq. 11 moved to Part B); estimated clone profiles (`--clone-profiles hatchet`, the real-data
  analogue) so `rck_mc` is not only an upper bound.
- **6.** Plot (phased − blind J_seg) and the within-block vs absolute counts against heterozygous SNP
  density and phase-block length; turns "alleles are identified within a phase block" into a measured
  behaviour.
- **7.** A nested S = 4 case drawn as truth / RCK bulk / RCK multi-clone / CoRAL, each cycle's segments
  coloured by allele; candidate Figure 1b.
- **8.** Either a parsimony tie-break that prefers fewer cycles per clone, or a short explanation in the text.

## Part A in detail: RCK-style evaluation of the existing runs (local, no Gurobi)

Scores every answer `allele_bench.py` already phased and scored (24 cases × rck_bulk / rck_mc native,
and every method × post-hoc rule), plus the truth itself as a control, in RCK's terms. Module
`rck/ecdna/rck_eval.py` (definitions in its docstring), tests `tests/ecdna/test_rck_eval.py`, driver
`rck_style.py` here, figure / table script `make_rck_style_figures.py`.

1. **Novel adjacencies, haplotype-specific** (RCK Fig. 2A): the answer's distinct labelled novel
   adjacencies matched one to one to the truth's within 100 bp (`aa_eval.adjacency_items`, phased).
   FNR = missed / truth; FDR = unmatched / answer's. RCK's FPR is over the input candidates
   (spurious inputs used / spurious inputs); here: graph candidates that match no truth junction,
   used by the answer / all such, haplotype-blind (candidates carry no labels). Unlabelled versions of
   FNR / FDR alongside. Within phase blocks: one A/B swap per block, the same swaps the copy-number
   score picks, applied to the cycles (`evaluate.relabel`).
2. **Eq. 11 copy-number distance** (RCK Fig. 2B), bulk level (n = 1 genome: the bulk ecDNA copies
   Σ_k u_k ψ_k, as items 2–3 score): Σ over phase blocks of min over the A/B swap of
   Σ_bases (|c_A − t_A| + |c_B − t_B| + c_U), in copies × bp, and relative to Σ t (so cases compare).
   Phase blocks are RCK's fragments. Per-clone Eq. 11 for `rck_mc` belongs to item 5.
3. **Generalized-IS violations** (RCK Fig. 4), over the union of the answer's labelled novel adjacencies
   (RCK imposes them on the whole sample, so intra- and inter-clone violations are both counted):
   - `is_ext`: an (unlabelled) extremity in two or more distinct novel adjacencies (extremity
     exclusivity, as `aa_eval.is_violations`);
   - `is_hom_same`: one novel adjacency used with two or more labellings (RCK Fig. 4A / 4C);
   - `is_hom`: an extremity used on A in one novel adjacency and on B in any (homologous-extremity
     exclusivity; includes `is_hom_same`);
   - `is_recip`: reference-adjacent extremities (end of p on "+", start of p + 1 on "−") used on
     opposite haplotypes (homologous-reciprocal exclusivity, RCK Fig. 4E).
   Within phase blocks after the same swaps. The truth is scored too: violations there are properties of
   the simulator (e.g. two species breaking the same bin), to be reported, not hidden.

Within-block frame, two scoring details (both in `rck_eval`'s docstring, both with unit tests): the block
swap is applied per junction end (the block holding that breakend), because a bulk-graph segment may cross
a block boundary; and a join of reference-contiguous steps is dropped when it lies on a block boundary or
carries one haplotype after the swaps. Without these, within-block scoring alone manufactured
homologous-reciprocal "violations" in RCK's own answers (6 multi-clone and 1 bulk answer, all at block
boundaries, 0 in the absolute frame). With them the within-block counts never exceed the absolute ones
(+1 at most, 2 of 336 answers, the other direction).

### Part A results (done 2026-11-16)

Note (decision after these results): the method under test is the multi-clone route; RCK bulk's rows are
kept as a reference and because CoRAL / CE decompose its graph, not as a second RCK-ecDNA result.

`python rck_style.py run && python rck_style.py collect && python make_rck_style_figures.py`: 360 rows
(336 answers + 24 truths), all ok, seconds with 16 processes; Figure 9 (`figures/fig9_rck_style.{pdf,png}`)
and Table 3 (`figures/table3.md`, `table3.tsv`). Means over cases, within phase blocks (Eq. 11: median):

| case group | method, phasing | labelled adj. FNR | labelled adj. FDR | Eq. 11 (rel.) | homologous | homologous-reciprocal |
|---|---|---|---|---|---|---|
| both alleles (6) | truth | – | – | – | 0.50 | 0.00 |
| | RCK bulk, own | **0.00** | **0.01** | **0.068** | 0.50 | **0.00** |
| | RCK bulk, oracle post-hoc | 0.35 | 0.48 | 0.340 | 0.00 | 4.00 |
| | RCK multi-clone\*, own | 0.21 | 0.05 | 0.000 | 0.17 | **0.00** |
| | RCK multi-clone\*, oracle post-hoc | 0.44 | 0.48 | 0.305 | 0.00 | 4.00 |
| | CoRAL†, BAF post-hoc | 0.40 | 0.41 | 2.651 | 0.50 | 4.33 |
| | CE†, BAF post-hoc | 0.41 | 0.40 | 1.587 | 0.00 | 3.50 |
| single allele (18) | RCK bulk, own / oracle | 0.01 / 0.00 | 0.01 / 0.00 | 0.049 / 0.049 | 0 / 0 | 0 / 0 |
| | CoRAL†, BAF post-hoc | 0.08 | 0.21 | 2.383 | 0.06 | 1.56 |

Reading, in RCK's terms:
- **Fig. 2A analogue.** Where both alleles are amplified, RCK-ecDNA bulk recovers every labelled novel
  adjacency (FNR 0.00, FDR 0.01). Phasing the *same* structures after the fact, even with the oracle rule,
  misses a third of them and a half of what it reports is wrongly labelled (0.35 / 0.48); CoRAL / CE land at
  ~0.40 / 0.40. Unlabelled FNR is 0.00 for RCK bulk in both cases, so the whole difference is phasing.
- **Fig. 2B analogue.** Eq. 11 distance rises five-fold when RCK bulk's structures are phased post hoc
  (0.068 → 0.340); the multi-clone route is exact (0.000) with its own phasing and 0.305 post hoc. CoRAL /
  CE are dominated by over-counted copies (median 1.6–2.7 in both groups, the same in single-allele cases),
  so Eq. 11 separates them from RCK by structure, not by phasing.
- **Fig. 4 analogue (the ReMixT argument, carried over).** Post-hoc phasing produces 3–4.5 homologous-
  reciprocal violations per mixed case, on every method's structures: it gives two reference-contiguous
  segments of one walk different alleles, i.e. an A-to-B join at a single position (e.g. hard_S2_s15003,
  chr16:85,415,999 A → 85,416,000 B), which needs breaks at homologous positions on both chromosomes, what
  RCK calls "highly unlikely". RCK-ecDNA has **0** in all 24 cases although it runs without `--enforce-is`;
  its homologous violations (0.50) are exactly the truth's own (below). In single-allele cases the BAF rule
  still produces 0.7–1.6 per case (wrong major-allele calls on noisy segments), the oracle rule none.
- My earlier expectation (this plan, first version of Part A: "post-hoc phasing ... creates IS violations
  only where neighbouring-step filling or a segment boundary splits differently") was wrong for mixed
  cases: a segment boundary where the major allele changes is common there, and each one is a violation.
- **Input-level FPR is not measurable on these inputs**: the union subgraphs hold 9 candidate junctions per
  case (median) and **none is spurious** in any of the 24 (the subgraph keeps junctions with ≥ 20 reads;
  the clone background's have ≤ 8). RCK's 10%-spurious condition needs Part B.
- **The truth is not IS-clean in 10 of 24 cases**, a simulator property to state in the methods:
  extremity reuse in 10 (within one species, from ecSimulator's duplication / translocation loop, e.g.
  hard_S1, and between nested species that break at the same 1 kb bin boundary), homologous reuse in 2
  (hard_S4_s16003: 2, s16005: 1; nested species on opposite alleles breaking at the same bin). RCK bulk's
  own homologous violations are exactly these two cases. Item 1's `is_clean` seed check prevents this for
  siblings; for bench_IW it is a known limit, not a scoring error.

Limits: n = 6 mixed cases, descriptive only; the multi-clone route has the true clone profiles; CoRAL / CE
are decomposition-only on RCK bulk's graph; Eq. 11 is bulk-level (per clone: item 5).

## Part B in detail: RCK's simulation inputs on the multi-clone route (decided 2026-11-16)

**Decisions (Yubi, 2026-11-16).** Part A + Part B together mimic, as closely as our setting allows, how
Aganezov & Raphael (2020) tested RCK. The method under test is the **multi-clone route** (`rck_mc`) on the
subgraph labelled with the tumour clones' copy-number profiles and clone proportions
(`subgraph_construct.py --clone-profiles truth`). Only **P = 0.9**; **no homogeneous S′ = (G1) cases**; the
**bulk route is kept only as the graph CoRAL / CE decompose** (its own scores stay in the tables for
reference, not as a second RCK-ecDNA result).

**Why the true clone profiles are faithful, not a shortcut.** HATCHet and Battenberg appear only in RCK's
real-data section. In the simulation study RCK was given the **true** clone- and allele-specific copy
numbers and the true novel adjacencies, then degraded on purpose ("Evaluation of RCK on simulated data";
Methods "Simulating rearranged cancer samples"): clone and haplotype information loss (CHIL), and errors.
Our labelled subgraph is that true input; Part B derives RCK's degraded versions from it. Running HATCHet
on the simulated reads (`--clone-profiles hatchet`) would be the real-data analogue; it stays out of
scope (item 5).

**RCK's own knobs already exist in our model.** `rck-ecdna --model multiclone` runs RCK's
`OptModelMultiClone` extended to ecDNA: `rck.fragments.tsv` defines RCK's fragments, the units whose A/B
the model may flip (`fragment_flipping`); today they are the phase blocks. `--run-nas-fp` (default 0.1) is
RCK's input-adjacency parameter, P = 1 − nas_fp, so the default is the paper's P = 0.9. RCK's molecule
groups (Methods "Third-generation sequencing technologies and novel adjacency groups") are our long-read
path groups (`rck.adjgroups.tsv`, `agt=M`).

### Input variants (RCK's Fig. 2, Supplemental Figs. S2–S3)

| variant | how it is made from the labelled input | RCK (2020) |
|---|---|---|
| `labelled` (reference) | the existing runs: true profiles in the observed phase frame (one unknown A/B per phase block), long-read junction haplotype tags, labelling groups | – (our setting: what long reads add) |
| `chil` | **CHIL**: one fragment per DIAG segment, each segment's A/B swapped with probability 0.5 (the same swap in every clone); junction haplotype tags (`junction_hap`, `hap_obs_end1/2`) and labelling groups (`agt=L`) removed; molecule groups (`agt=M`) kept | error-free input with CHIL |
| `chil_avg` | each clone's allele-specific profile averaged over 50 kb fragments (fragments = 50 kb windows), then CHIL | fragment averaging only (Supplemental Fig. S3) |
| `chil_err` | 50 kb averaging; ±1 copy on a random 5% of fragments (per clone and allele); ±50 bp on the breakends of a random half of the novel adjacencies; 10% spurious novel adjacencies added; then CHIL | the full error model (Fig. 2) |

3 new variants × 24 cases × P = 0.9 = **72 `rck_mc` runs**, one Slurm array with `%1` (Gurobi licence),
after Yubi's go-ahead. The `labelled` runs exist (bench_IW_intpsi).

**How each transform is made** (a seeded script, `make_rck_inputs.py` here; outputs under
`/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/<case>/<variant>/rck_input/`; the bench_IW_intpsi subgraphs are
only read):
- Copy-number transforms (averaging, ±1, CHIL flips) edit a copy of `rck.scnt.tsv` / `rck.fragments.tsv`
  (and `diag_segments.tsv` for the record), in RCK's order: average, perturb, flip. Averaged values are
  rounded to integers if `read_scnt` requires it (to check); scnb bounds are kept.
- Adjacency transforms (±50 bp, spurious) move breakpoints, so segments must be re-aligned to them
  (`multiclone.py` checks `positions_aligned`). They are applied **upstream**: to a copy of the case's
  observed adjacency table, then the DIAG is rebuilt with `subgraph_construct.py` and the same flags
  (`bench_IW.py` `--cnp` flags, `--clone-profiles truth`), which aligns segments exactly as for real data.
  Spurious adjacencies get read support drawn from the true junctions' (so they pass the ≥ 20-read filter;
  RCK's are "in the input"), with breakends inside the amplicon region. The rebuild is the one step that is
  not a few seconds; if all 24 together take more than a few minutes they run as a Slurm job.
- Every transform gets a unit test against the paper's definition (flip probability, 5% share, ±1, ±50 bp,
  10% count, 50 kb windows), and a check that the `chil` input with flips undone equals the labelled
  input's copy numbers.

### Scores (Part A's, plus what RCK reports per clone)

- Haplotype-specific novel-adjacency FNR and FPR (RCK Fig. 2A). FPR over the input: spurious input
  adjacencies used / spurious input adjacencies; defined in `chil_err` only (the clean inputs hold none,
  Part A). FDR alongside.
- **Eq. 11 summed over clones**, as RCK defines it: per clone i, the inferred allele-specific segment copy
  numbers against the truth's clone karyotype (`truth/rck/`), over the input's fragments, minimum over the
  A/B swap per fragment. New in `rck_eval` (where `rck_mc` writes its per-clone segment copy numbers is to
  be checked). Bulk-level Eq. 11 (Part A) kept alongside for CoRAL / CE.
- Copy-number state agreement (RCK Supplemental Fig. S2B,C): per clone and allele, loss / neutral /
  gain, share of base pairs in the same state as the truth.
- Generalized-IS violations (Part A), and the truth's own as the control.
- RCK's ReMixT section, carried to CoRAL / CE (as the comparator, on their existing runs): IS violations
  (done, Part A), share of input adjacencies used (RCK Fig. 3B), and linear walks (open paths) as the
  ecDNA analogue of novel telomeres (RCK Supplemental Figs. S6, S9).

**Differences from RCK, to state in the methods:** ecDNA cycles in the amplicon region instead of whole
karyotypes; 24 cases instead of 100; RCK's simulation compared RCK with no other method (ReMixT only on
real data), so CoRAL / CE stand in for ReMixT's role on simulated data; long-read junction tags (`labelled`)
have no RCK counterpart and are reported as our setting next to RCK's `chil`.

### Part B status (2026-11-16)

**Step 1 done** (`make_rck_inputs.py`, transforms in `rck/ecdna/rck_inputs.py`, 8 unit tests): 72 inputs under
`/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/<case>/<variant>/rck_input/` (+ `transforms.json`), all loading
with `multiclone.load_rck_input` (RCK's readers and alignment checks); counts in `rck_inputs_summary.tsv`.
Built in 3 s locally. Choices made while building, each recorded in the module docstring:
- **No subgraph rebuild.** Moved and spurious breakends are aligned by splitting the segment they fall in
  (copy numbers copied), as RCK aligns segment extremities to novel adjacency positions (Methods "Deriving
  extremities and novel adjacencies from data"). Simpler than the rebuild planned above and stays in RCK's
  input format; the amp repo is untouched.
- **50 kb windows only where they change something.** The DIAGs span up to ~150 Mb (median segment 1 Mb;
  99.7% of the bases in segments > 50 kb), so cutting everything into 50 kb fragments would make the model
  ~15x larger (hard_S4_s16003: 214 → 3,138 segments) with no change in copy numbers. A window overlapping two
  or more segments is cut and averaged; one segment alone over whole windows stays one fragment. Result:
  214 → 242 (`chil_avg`) / 267 (`chil_err`) segments for hard_S4_s16003; medians 71 → 92 / 111.
- **Averaged copy numbers stay real-valued** (3 decimals; RCK's readers accept them; whether the solve and
  the decomposition do is what the pilot checks).
- **CHIL flips are per fragment** (= per segment in `chil`), one swap shared by both clones; RCK flips the
  fragment-level profile the same way.
- **Spurious adjacencies: max(1, round(10% of n))**, so 1–2 per case (3–19 true junctions); above 10% in the
  smallest cases. Breakends at random positions of the DIAG's non-flank segments (mostly non-amplified
  background, like RCK's random spurious adjacencies), ≥ 1 kb from other breakends, read support copied from
  a true junction.
- **Labels:** the multi-clone ILP takes junction haplotype information only through labelling groups
  (`agt=L`, 1–9 per case, all dropped); the `junction_hap` / `hap_obs_end*` tags are dropped too. Molecule
  groups (`agt=M`, 0–13 per case) are kept.

**Step 2 done** (`rck_eval`: `clone_cn_scores`, `multiclone_input_use`, `graph_input_use`, `linear_walks`; 4 new
unit tests; `rck_style.py` re-run, seconds). On the existing (`labelled`) multi-clone runs: per-clone Eq. 11
is 0 in all 18 single-allele cases and 0.042 relative on the ecDNA bases in the 6 mixed ones (max 0.25),
state agreement ≥ 0.875 on ecDNA bases (mean 0.98), every input junction used, no linear walks. That is
RCK's clean-input finding ("karyotypes that use all input novel adjacencies and whose haplotype-labeled copy
numbers agree with the simulated copy numbers") and Part B's baseline. CoRAL / CE use 77–96% of the
junctions in the graph they decompose and leave 0.7 (CoRAL) / 1.2 (CE) linear walks per case when walks
are kept.

**Steps 3–4 prepared, not run:** `rck_partB.py` (list / run / status / export), `rck_partB.sbatch`,
`submit_partB.sh` (prints the command unless `--submit`). Settings as bench_IW_intpsi's rck_mc
(`bench_iw.rck_argv`) except 600 s per MILP instead of 60 (the reference solved in ≤ 5 s; label-free inputs
may need longer; each run records its solver status). Array `0-71%1`, 30 min, 2G. Proposed first: a pilot
of 3 tasks (21, 45, 69 = hard_S4_s16003, the largest DIAG, in each variant) to check real-valued copy
numbers and solve times before the full array.

**Pilot done 2026-11-16** (array 32364762, tasks 21 / 45 / 69 = hard_S4_s16003, the largest DIAG, in each
variant; submitted 21:07 EST with no other array holding the licence, finished by ~21:09): all three optimal
(status 2, gap 0) on the first attempt, 8.8 s (`chil`), 16.6 s (`chil_avg`), 34.0 s (`chil_err`) per task
(model solve 5.0 / 10.6 / 26.3 s; reference 4.9 s), far under the 600 s limit. Real-valued averaged copy
numbers pass through the solve, the decomposition and `rck_partB.py export`. Every true junction used; the 2
spurious ones in `chil_err` not used. Cycles in the decomposition: 3 / 7 / 8 (reference 7). Not scored yet.
**Array submitted 2026-11-16 21:10 EST** with Yubi's go-ahead: 32364902, `0-71%1` (the 3 pilot runs skip as
done). **Step 4 scripts written** (`score_partB.py`, `make_partB_figures.py`; committed `ba4a6f5`), run end
to end on the 30 rows available at hand-off. Scoring frame: within phase blocks for every variant (CHIL
flips A/B per fragment, so swapping per fragment to match the truth would forgive exactly the phasing errors
CHIL tests; phase across blocks is unknowable from any data); Eq. 11 also in each variant's own fragments
(RCK's definition as written). Remaining: "Resume here" steps 1–6.

### Steps
1. **Local, seconds:** `make_rck_inputs.py` (copy-number transforms; adjacency transforms + DIAG rebuild,
   as a Slurm job if over a few minutes) and its tests; check one case of each variant by hand.
2. **Local, seconds:** per-clone Eq. 11, state agreement, input-adjacency use and linear-walk counts in
   `rck_eval`; score the existing `labelled` runs and CoRAL / CE with them.
3. **Slurm array (Gurobi, Yubi's go-ahead):** 72 `rck_mc` runs, `%1`, sized as bench_IW_intpsi's rck_mc
   array; check at ≥ 60 s intervals.
4. **Scoring:** `rck_style.py` gets the Part B variants; Figure 10 laid out as RCK's Fig. 2
   (adjacency FNR / FPR and Eq. 11 per variant), Table 4 per variant; the CoRAL / CE comparison as
   RCK's Figs. 3–4.

## Item 1 in detail: the allelic-sibling benchmark (design, not started)

**Goal.** Cases where both alleles of one locus are amplified on different ecDNA species, so the
allele-specific model is necessary by construction, while respecting RCK's generalized IS: species 1
plus a **sibling** on the other allele of the same reference intervals, with its **own** breakpoints and
rearrangement. A haplotype-blind method sees overlapping species with summed copy numbers on shared
segments; post-hoc phasing can give each base only one allele, so it recovers at most the major share.
Varying the minor sibling's share gives a curve (prediction: post-hoc ≈ 1 − minor share on sibling
bases; RCK flat if it separates them).

**What the simulator does today** (amplicon-allele-analysis `simulations/src/structures.py`,
`build_structures`): every step of a species walk carries a haplotype. With `hap_labels: inherit`
(bench_IW's setting) `inherit_haplotypes` draws one haplotype per source interval, switching between
intervals with probability `p_cross` (0.2). Nested species (`nest_species: true`) are placed inside
species 1's intervals (`select_nested_intervals`) with their own rearrangement and haplotypes;
`nest_same_haplotype` forces the nested species onto species 1's starting haplotype.

**Simulator change** (defaults leave every existing config byte-identical; a test checks it):
1. `l2_structures.allelic_sibling`: a per-species list, `null` or the 0-based index of the source species.
   A sibling takes the source's `intervals` (or, with `sibling_footprint` < 1, nested intervals inside them,
   as `select_nested_intervals`), then draws its own breakpoints (`assign_breakpoints`) and rearrangement.
2. Haplotypes: each sibling segment gets the opposite of the haplotype the source carries on the source
   interval that contains it (inherit mode gives one per interval), so the sibling follows the source's own
   crossings on the other allele.
3. Per-species clone weights (`ecdna_clone_weights` is per clone and shared by every species today): an
   optional per-species override, so the A-sibling can sit mostly in clone 1 and the B-sibling mostly in
   clone 2, the ecDNA form of mirrored-subclonal allelic imbalance.

**Scenarios** (make_scenarios.py set `iw_allelic`, `integer_psi: true`, otherwise `iw_intpsi`):
- `sib_S2`: species 1 + a full-footprint sibling; minor share on a grid 0.1 / 0.2 / 0.3 / 0.5 on the same
  structures (the structures layer has its own random stream, `layer_rng(cfg, "structures")`, so
  changing `ecdna_proportions` keeps the walks); total copies fixed (~40) so the minor sibling keeps
  ≥ 4 copies and ≥ 20 reads per junction.
- `sib_S2c`: as `sib_S2` with the clone split (A-sibling clone 1, B-sibling clone 2).
- 4 structure seeds × 4 shares per scenario (32 cases, or 16 for `sib_S2` alone to start).

**Seed checks** (bench_IW.py `--set iw_allelic`): `sibling` (the sibling's bases lie in the source's;
opposite allele on ≥ 90% of them), **`is_clean` (no breakpoint of the sibling within 1 bin of a source
breakpoint, so the truth obeys RCK's generalized IS; Part A's `is_*` counts on the truth must be 0)**,
`placed`, `footprint`, ≥ 20 reads per junction, bulk ratio. Enough heterozygous SNPs on the intervals
(~1 per 1.8 kb in truth/phasing.tsv; intervals ≥ 100 kb are fine).

**Runs and scoring:** `bench_iw.py` `SETS` entry `iw_allelic`; arrays for `rck_mc` (the method under
test) and `rck_bulk` (only to build the graph CoRAL / CE decompose), then CoRAL / CE on that graph (`%1` per Gurobi array), **after Yubi's go-ahead**, preceded by a local pilot on
2–3 cases (also Gurobi, so also needs the go-ahead and a free licence). Scored with `allele_bench.py`
(`--set`), Part A's metrics, and phased species matching for proportions (`evaluate._match(...,
with_hap=True)`; `bench_metrics.proportion_scores` matches blind).

**Expected figure:** x = minor-sibling share; y = allele-specific copy-number recall, Eq. 11 distance and
phased species F1; lines for RCK bulk / multi-clone (own), post-hoc oracle, CoRAL / CE (BAF). If RCK does
*not* separate siblings (it may prefer one species carrying both alleles where the step penalty favours
fewer cycles), that is a finding too; the pilot answers it before any array.

## Results (items 2 and 3, 2026-11-16)

**Key fact first.** Only **6 of the 24 cases** (all nested: hard_S2_s15003, s15005, hard_S4_s16002,
s16003, s16005, s16006) have ecDNA bases where **both alleles are amplified** (5–70% of their ecDNA bases,
`mixed_bp_frac`). In the other 18, every ecDNA base carries copies of one allele only. On those, any rule
that gives each segment its major allele is exactly right, so phasing after the decomposition is as good
as phasing inside it. The allele-specific model can only show its value where alleles are mixed, which is
why item 1 matters.

Post-hoc phasing rules (`rck.ecdna.allele_metrics.rephase`), applied to CoRAL / CE answers and, as a
control, to RCK's own structures:
- **BAF**: each segment gets the major allele of the subgraph's observed CN_A / CN_B (`segment_blocks.tsv`,
  from depth and phased SNP counts, fitted before any decomposition: what a CoRAL user with phased SNPs has);
- **oracle**: the major allele of the true ecDNA copies; the best any one-allele-per-segment rule can do;
- segments with no informative SNPs take the nearest called step's allele in the same cycle (8,995 steps
  filled, 3 left unphased over all 336 answers).

Table 2 (`figures/table2.md`; within phase blocks, mean over cases):

| | single-allele cases (18) | | both alleles amplified (6) | | |
|---|---|---|---|---|---|
| method, phasing | allele CN recall | phased J_adj | allele CN recall | ecDNA BAF error | phased J_adj |
| RCK bulk, own | 0.97 | 0.95 | **0.95** | **0.012** | 0.65 |
| RCK bulk, oracle post-hoc | 0.97 | 0.95 | 0.82 | 0.130 | 0.40 |
| RCK multi-clone\*, own | 1.00 | 0.89 | **0.95** | 0.028 | **0.80** |
| RCK multi-clone\*, oracle post-hoc | 1.00 | 0.89 | 0.81 | 0.130 | 0.47 |
| CoRAL†, BAF post-hoc | 0.89 | 0.79 | 0.69 | 0.129 | 0.26 |
| CE†, BAF post-hoc | 0.99 | 0.82 | 0.76 | 0.130 | 0.30 |

Reading:
- **Where both alleles are amplified, phasing after the decomposition loses ~15–20% of the ecDNA copies
  to the wrong allele, even with the oracle rule and even on RCK's own structures** (0.95 → 0.82). The ecDNA
  BAF error grows tenfold (0.01 → 0.13) and phased J_adj drops (0.65 → 0.40 bulk; 0.80 → 0.47 multi-clone).
  So the allele assignment has to happen inside the decomposition; it cannot be added afterwards.
- With the same post-hoc phasing, CoRAL and CE are below RCK's structures (allele CN recall 0.69 / 0.76 vs
  0.82), i.e. their structural errors add to the phasing loss.
- **Single-allele cases: post-hoc phasing matches RCK's own phasing** (oracle = own on every column; BAF
  only loses junction-level detail, J_adj 0.84 vs 0.95). Say this in the paper rather than leave a
  reviewer to find it: the allele model matters where alleles are mixed.
- Figure 8: RCK's per-segment allele copy numbers sit on the diagonal, including the bases where both
  alleles are amplified (open markers); CoRAL / CE phased post-hoc put the minor allele at 0 (the row at
  the bottom) and carry their structural copy-number errors.

Limits: n = 6 mixed cases (descriptive only, no tests); the multi-clone route has the true clone profiles;
CoRAL / CE are decomposition-only on RCK bulk's graph, as in figures 2–6; within-block scores choose the
best swap per phase block for every method alike (absolute scores are in `allele_table.tsv`).

Checks that the new scorer is right: native RCK rows reproduce `bench_iw_intpsi_table.tsv` exactly
(phased J_seg, exact cycles, composite, species counts; rck_bulk absolute, rck_mc within-block), and
rephasing leaves every haplotype-blind CoRAL / CE score unchanged (difference 0 on all 24 cases);
`blind_cn_recall` equals the table's `cn_recall` on all 24. Unit tests: `tests/ecdna/test_allele_metrics.py`.

## Files

| what | where |
|---|---|
| allele-specific scorer and post-hoc phasing | `rck/ecdna/allele_metrics.py` (docstring has the definitions); tests `tests/ecdna/test_allele_metrics.py` |
| RCK-style metrics (Part A: adjacency FNR / FDR / FPR, Eq. 11, generalized-IS violations) | `rck/ecdna/rck_eval.py`; tests `tests/ecdna/test_rck_eval.py` |
| Part B inputs (CHIL, averaging, errors) | `rck/ecdna/rck_inputs.py` (tests `tests/ecdna/test_rck_inputs.py`); driver `make_rck_inputs.py` → `/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/<case>/<variant>/rck_input/`; counts `rck_inputs_summary.tsv` |
| Part B runs (multi-clone on the degraded inputs) | `rck_partB.py`, `rck_partB.sbatch`, `submit_partB.sh` → `<case>/<variant>/rck_mc/`; logs `/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/logs/slurm/` (`submitted.txt` lists the arrays) |
| Part B scoring, Figure 10, Table 4 | `score_partB.py` (→ `<run>/score/partB.json`, `partB_table.tsv`), `make_partB_figures.py` (→ `figures/fig10_partB.*`, `figures/table4.*`) |
| driver, items 2–3 (phase, score, collect) | `allele_bench.py` |
| driver, Part A (score the item 2–3 answers and the truth, collect) | `rck_style.py` |
| figures 7, 8, Table 2 | `make_allele_figures.py` → `figures/` (`figure_data/` has the plotted numbers) |
| Part A figure and Table 3 | `make_rck_style_figures.py` → `figures/` |
| per-run outputs (cycles, aa_eval, aa_eval_block, allele/, rck_style/) | `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele/<case>/<method_variant>/<phasing>/`; truth control in `<case>/truth/rck_style/` |
| collected tables | `allele_table.tsv` (one row per case × method variant × phasing; absolute and `_block` columns), `allele_segments.tsv` (per-interval t_A, t_B, c_A, c_B, c_U), `rck_style_table.tsv` (Part A, one row per case × answer, truth rows included) |
