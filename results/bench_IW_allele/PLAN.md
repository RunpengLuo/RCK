# bench_IW_allele: analyses that show RCK-ecDNA's allele-specific reconstruction (branch `yubi_CNP_ionic`)

*Started 2026-11-16 (Yubi: "write down your ideas in a document in the repo that can be referred to in
future sessions. then, start items 2 and 3 ... make sure to describe how [item 1] can be implemented").
Revised 2026-11-16 (later the same day) after reading the RCK paper's evaluation: the plan now follows
RCK's own evaluation format first (Part A, Part B), and item 1 became the **allelic-sibling** benchmark
instead of allelic twins (reasons under "Alignment with the literature"). The end-to-end comparison
([`../bench_IW_e2e/NEXT_STEPS.md`](../bench_IW_e2e/NEXT_STEPS.md)) is paused by Yubi on 2026-11-16,
blocked on the purity / ploidy fit.*

## Resume here

1. **Items 2 and 3 are done** (below, "Results"): `allele_bench.py run` / `collect`, then
   `make_allele_figures.py` (figures 7, 8 and Table 2 in `figures/`). All local, no Gurobi, no Slurm.
   To redo from scratch: delete `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele/` (finished tasks are skipped).
2. **Part A is done** (below, "Part A results"; Figure 9, Table 3): `rck_style.py run` / `collect`, then
   `make_rck_style_figures.py`. Local, no Gurobi; seconds. To redo: delete the `rck_style/` folders under
   `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele/` (finished tasks are skipped).
3. **Next, waiting for Yubi:** Part B (RCK-style perturbed and label-free inputs; RCK runs, so Gurobi and
   an array, which need Yubi's go-ahead), then item 1 (allelic siblings: simulator option first, local and
   cheap; method runs need arrays and Yubi's go-ahead).
4. Write-up: the paragraph under "Results" is ready to adapt; RESULTS.md-style captions are in the
   figure scripts' docstrings.

Rules for this line of work (Yubi's, carried over): keep everything in `results/bench_IW_allele/` and
separate from `results/bench_IW/` and the bench_IW_intpsi runs (only read); new Slurm arrays and new Gurobi
work need Yubi's go-ahead; Gurobi licence 2 sessions (`%1` per Gurobi array, no local Gurobi while arrays
hold both); local jobs longer than ~30 min go to Slurm on ionic; wait ≥ 60 s between `squeue` / `sacct`
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
| CHIL on segments (A/B hidden per segment) | allele copy numbers from phased SNPs in phase blocks of unknown orientation (A/B hidden per block) | done (simulator / subgraph) |
| CHIL on adjacencies (no labels) | **not yet**: our junction ends carry long-read haplotype labels, used as hard constraints | **Part B** (label-free variant) |
| errors in breakpoints, spurious adjacencies, copy numbers | read-sampling noise only | **Part B** (perturbed variants) |
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
| B | **RCK-style inputs**: label-free junctions (CHIL) and perturbed inputs, on the 6 mixed-allele cases | planned, needs go-ahead | 36 new RCK runs (array, `%1`) |
| 1 | **allelic-sibling cases**: a second species on the other allele of the same locus, own breakpoints | designed (below), not started | simulator option + seeds; method arrays |
| 2 | **post-hoc phasing baseline**: CoRAL / CE (and RCK's own structures) phased afterwards by a rule | **done 2026-11-16** | local, minutes |
| 3 | **allele-specific copy numbers** per segment: true vs inferred CN_A / CN_B, BAF error | **done 2026-11-16** | local, minutes |
| 4 | ablation: RCK with the alleles collapsed (allele-blind graph) | idea | small array |
| 5 | clone-level results of the multi-clone route (ψ per clone; Eq. 11 per clone), and estimated clone profiles | idea (bench_IW NEXT_STEPS item 2) | local scoring; profile perturbation runs |
| 6 | when phasing fails: phased − blind scores vs SNP density / phase-block length | idea | local on existing runs; a SNP-density sweep needs seeds |
| 7 | one worked example figure: truth vs each method, segments coloured by allele | idea | local |
| 8 | explain or fix the multi-clone weak spots (hard_S4_s16003; extra split cycles) | open (bench_IW NEXT_STEPS item 6) | local |

Order: A (in RCK's own terms, on what exists), then B (completes RCK's simulation design), then 1 (the
case RCK's evaluation never isolates: both alleles of a locus amplified).

Items 4–8 in a sentence each:
- **4.** Run `rck-ecdna` on the 24 cases with A and B merged in the input (or the model's allele split
  switched off, if it has a switch); the drop on nested S = 4 cases shows what the allele model buys.
- **5.** Score ψ_{i,s} (copies of species s per cell of clone i) against the truth as a scatter or heatmap,
  and Eq. 11 summed over clones as RCK does; pair with perturbed truth profiles (±1 copy, wrong u) so
  `rck_mc` is not only an upper bound.
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

## Part B in detail: RCK-style inputs (needs Yubi's go-ahead; Gurobi)

On the 6 cases with both alleles amplified (hard_S2_s15003, s15005, hard_S4_s16002, s16003, s16005,
s16006), 4 input variants × 2 routes (rck_bulk, rck_mc) = 48 runs, as an array with `%1`:

| variant | what changes | RCK analogue |
|---|---|---|
| clean | nothing (re-uses the existing runs; no new run) | error-free input |
| nolabel | junction haplotype labels removed (every labelling allowed; `junction_hap` emptied in a copy of edges.tsv, not `--soft-junction-hap`, which only turns labels into penalties) | CHIL on adjacencies |
| noisy | ±1 copy on 5% of segment pieces (CN_A or CN_B), and 10% spurious novel adjacencies (random breakend pairs inside the region, with read support drawn from the observed junctions' lower quartile) | errors in input |
| nolabel+noisy | both | CHIL + errors |

So 3 new variants × 6 cases × 2 routes = 36 new runs. Perturbations are written to copies of the
subgraph under `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele_partB/` by a seeded script; scored with the
Part A metrics plus items 2–3. Open questions before building: how `rck_mc`'s DIAG input
(`rck_input/`) carries junction labels (labelling groups, `--diag-labeling-groups block`) and how to strip
them there; whether spurious adjacencies should be ecDNA-region-only (as here) or also reach the flanks.
Part B does not need Part A's results to run, only its scorer to be scored (done). Part A showed
the clean inputs have no spurious candidates at all, so the noisy variant is the only place RCK's FPR
can be measured.

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

**Runs and scoring:** `bench_iw.py` `SETS` entry `iw_allelic`; arrays for `rck_bulk`, `rck_mc`, then CoRAL /
CE on RCK bulk's graph (`%1` per Gurobi array), **after Yubi's go-ahead**, preceded by a local pilot on
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
| driver, items 2–3 (phase, score, collect) | `allele_bench.py` |
| driver, Part A (score the item 2–3 answers and the truth, collect) | `rck_style.py` |
| figures 7, 8, Table 2 | `make_allele_figures.py` → `figures/` (`figure_data/` has the plotted numbers) |
| Part A figure and Table 3 | `make_rck_style_figures.py` → `figures/` |
| per-run outputs (cycles, aa_eval, aa_eval_block, allele/, rck_style/) | `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele/<case>/<method_variant>/<phasing>/`; truth control in `<case>/truth/rck_style/` |
| collected tables | `allele_table.tsv` (one row per case × method variant × phasing; absolute and `_block` columns), `allele_segments.tsv` (per-interval t_A, t_B, c_A, c_B, c_U), `rck_style_table.tsv` (Part A, one row per case × answer, truth rows included) |
