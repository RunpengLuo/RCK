# bench_IW_allele: analyses that show RCK-ecDNA's allele-specific reconstruction (branch `yubi_CNP_ionic`)

*Started 2026-11-16 (Yubi: "write down your ideas in a document in the repo that can be referred to in
future sessions. then, start items 2 and 3 ... make sure to describe how [item 1] can be implemented").
The end-to-end comparison ([`../bench_IW_e2e/NEXT_STEPS.md`](../bench_IW_e2e/NEXT_STEPS.md)) is paused
by Yubi on 2026-11-16, blocked on the purity / ploidy fit.*

## Resume here

1. **Items 2 and 3 are done** (below, "Results"): `allele_bench.py run` / `collect`, then
   `make_allele_figures.py` (figures 7, 8 and Table 2 in `figures/`). All local, no Gurobi, no Slurm.
   To redo from scratch: delete `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele/` (finished tasks are skipped).
2. **Next, waiting for Yubi:** item 1 (the allelic-twin benchmark, design under "Item 1 in detail").
   Its first step is local and cheap (simulator option + a few seeds); the method runs need arrays,
   which need Yubi's go-ahead.
3. Write-up: the paragraph under "Results" is ready to adapt; RESULTS.md-style captions are in the
   figure script's docstrings.

Rules for this line of work (Yubi's, carried over): keep everything separate from `results/bench_IW/`
and the bench_IW_intpsi runs (only read); new Slurm arrays need Yubi's go-ahead; Gurobi licence 2 sessions
(`%1` per Gurobi array); commit on `yubi_CNP_ionic`, don't push.

## Why these analyses

The current manuscript results (`results/bench_IW/RESULTS.md`, figures 2–6 and Table 1) mostly show
**structure and copy-number** accuracy, where RCK-ecDNA wins clearly. The allele-specific claim rests
on Figure 6 alone, which has no baseline, and on it phased and allele-blind scores agree in 37 of 48 runs.
A reviewer can read that as: "in this benchmark, getting the alleles right comes for free once the
structure is right; a haplotype-blind method could assign alleles afterwards." The analyses below test
that directly. Ranked by how much they strengthen the allele-specific claim:

| # | analysis | status | cost |
|---|---|---|---|
| 1 | **allelic-twin cases**: species that differ only (or partly) by allele | designed (below), not started | simulator option + new seeds; method arrays |
| 2 | **post-hoc phasing baseline**: CoRAL / CE (and RCK's own structures) phased afterwards by a rule | **done 2026-11-16** | local, minutes |
| 3 | **allele-specific copy numbers** per segment: true vs inferred CN_A / CN_B, BAF error | **done 2026-11-16** | local, minutes |
| 4 | ablation: RCK with the alleles collapsed (allele-blind graph) | idea | small array |
| 5 | clone-level results of the multi-clone route (ψ per clone), and estimated clone profiles | idea (bench_IW NEXT_STEPS item 2) | local scoring; profile perturbation runs |
| 6 | when phasing fails: phased − blind scores vs SNP density / phase-block length | idea | local on existing runs; a SNP-density sweep needs seeds |
| 7 | one worked example figure: truth vs each method, segments coloured by allele | idea | local |
| 8 | explain or fix the multi-clone weak spots (hard_S4_s16003; extra split cycles) | open (bench_IW NEXT_STEPS item 6) | local |

Items 4–8 in a sentence each:
- **4.** Run `rck-ecdna` on the 24 cases with A and B merged in the input (or the model's allele split
  switched off, if it has a switch); the drop on nested S = 4 cases shows what the allele model buys.
- **5.** Score ψ_{i,s} (copies of species s per cell of clone i) against the truth as a scatter or heatmap;
  pair with perturbed truth profiles (±1 copy, wrong u) so `rck_mc` is not only an upper bound.
- **6.** Plot (phased − blind J_seg) and the within-block vs absolute counts against heterozygous SNP
  density and phase-block length; turns "alleles are identified within a phase block" into a measured
  behaviour.
- **7.** A nested S = 4 case drawn as truth / RCK bulk / RCK multi-clone / CoRAL, each cycle's segments
  coloured by allele; candidate Figure 1b.
- **8.** Either a parsimony tie-break that prefers fewer cycles per clone, or a short explanation in the text.

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

## Item 1 in detail: the allelic-twin benchmark (design, not started)

**Goal.** Cases where ecDNA species differ only, or partly, in which parental allele they carry, so the
allele-specific model is necessary by construction. A haplotype-blind method sees one species at the
summed copy number; RCK should separate them. Varying the minor twin's share gives a curve: any post-hoc
rule recovers at most the major share on twinned bases, so its allele-specific recall falls as the minor
share grows toward 0.5, and RCK's should stay flat.

**What the simulator does today** (amplicon-allele-analysis `simulations/src/structures.py`,
`build_structures`): every step of a species walk carries a haplotype. With `hap_labels: inherit`
(bench_IW's setting) `inherit_haplotypes` draws one haplotype per source interval, switching between
intervals with probability `p_cross` (0.2). Nested species (`nest_species: true`) are placed inside
species 1's intervals (`select_nested_intervals`) with their own rearrangement and haplotypes;
`nest_same_haplotype` forces the nested species onto species 1's starting haplotype. Nothing yet makes a
species that copies another's structure on the other allele.

**Simulator change** (one new option, defaults leave every existing config byte-identical):
1. `l2_structures.allelic_twin`: a per-species list (like `nest_size_fraction`), `null` or the index of
   the source species. A twin copies the source's `intervals`, `segments` and rearranged `steps`, and
   flips the haplotype (A ↔ B) of the steps on a chosen subset of the source's intervals.
2. `l2_structures.twin_flip_fraction` (per species, default 1.0): the share of the source's intervals
   flipped. 1.0 = full twin (same structure, other allele everywhere); 0.5 = partial twin (structures equal,
   alleles differ on half the intervals, so haplotype-crossing junctions differ too).
3. Make no RNG draws for a twin (or draw the subset from a separate `layer_rng(cfg, "twins")`), so species,
   clones and counts after it keep their random stream; the file already keeps this rule ("NB: the same
   draws ...") and it keeps old seeds reproducible.
4. Copies: the twin gets its own entry in `ecdna_proportions` (iw_draw draws one per species), so the
   minor share is set per seed; add a bench variant that fixes it on a grid (0.1, 0.2, 0.3, 0.5).
   Optionally give twins different clone weights (`ecdna_clone_weights`): A-twin mostly in clone 1, B-twin
   mostly in clone 2, the case the multi-clone route is built for.

**Check before building** (each could block the design):
- Two novel adjacencies at the same coordinates with different phasings (AA in the source, BB in the
  twin): `novel_multiset` keys junctions with the haplotype (`left[2]`), but check `emit_counts` /
  `obs/adjacency*.tsv` and `subgraph_construct.py` keep them as one adjacency with per-haplotype counts
  (`NOV_A` / `NOV_B`, `adjacency_allele_ends.tsv`) rather than merging or dropping one.
- Species matching for proportions (`bench_metrics.proportion_scores`) uses **blind** Jaccard, so a full
  twin and its source are indistinguishable there. Add a phased matching option (`ev._match(...,
  with_hap=True)`) and report proportions both ways; the species matrix (`aa_species_matrix.tsv`) already
  has phased rows.
- The seed checks (`bench_IW.py` / `IW_INTENT`): a new `twin` check (the twin's bases are the source's;
  every twinned junction has ≥ 20 reads on each haplotype), `placed`, `footprint`.
- Enough heterozygous SNPs on the twinned intervals to phase them (truth/phasing.tsv has ~1.73 M heterozygous SNPs over the genome, ~1 per 1.8 kb;
  intervals ≥ 100 kb are fine).

**Benchmark** (make_scenarios.py, a new set `iw_allelic` next to `iw_intpsi`, `integer_psi: true`):
- `twin_S2`: species 1 + a full twin (2 species).
- `twin_S2p`: species 1 + a partial twin (flip fraction 0.5).
- `twin_S3`: species 1 + a full twin + a nested species (combines nesting and twinning).
- 6 seeds each = 18 cases, minor-twin share spread over 0.1–0.5. Each case is cheap to simulate (counts
  layer only; no reads). Subgraphs (bulk and multi-clone) with `bench_IW.py --set iw_allelic`.

**Runs and scoring:**
- `bench_iw.py` gets a `SETS` entry `iw_allelic`; arrays for `rck_bulk`, `rck_mc`, then CoRAL / CE on
  RCK bulk's graph, sized as `submit_intpsi.sh` (`%1` per Gurobi array). **Needs Yubi's go-ahead.**
- `allele_bench.py` gets `--set` and runs unchanged (post-hoc phasing, allele CN, phased species scores).

**Expected figure:** x = minor-twin share; y = allele-specific copy-number recall and phased species F1;
lines for RCK bulk / multi-clone (own), post-hoc oracle, CoRAL / CE (BAF). Prediction: post-hoc recall
≈ 1 − (minor share) on twinned bases, RCK near 1; blind methods report one species per twin pair.
If RCK does *not* separate full twins, that is a finding too (the bulk model may tie between "two
species on A and B" and "one species with both alleles"; the step penalty and the haplotype-labelled
junction counts decide it). Find out on 2–3 local cases before any array.

**Cost estimate:** simulator change and tests ~half a day; 18 cases simulate in minutes; RCK runs are
minutes per case (rck_bulk up to its 30-min limit); the CoRAL / CE arrays are the slow part, as in
bench_IW_intpsi.

## Files

| what | where |
|---|---|
| allele-specific scorer and post-hoc phasing | `rck/ecdna/allele_metrics.py` (docstring has the definitions); tests `tests/ecdna/test_allele_metrics.py` |
| driver (phase, score, collect) | `allele_bench.py` |
| figures 7, 8, Table 2 | `make_allele_figures.py` → `figures/` (`figure_data/` has the plotted numbers) |
| per-run outputs (cycles, aa_eval, aa_eval_block, allele/) | `/n/fs/yubi-iw/runs/bench_IW_intpsi_allele/<case>/<method_variant>/<phasing>/` |
| collected tables | `allele_table.tsv` (one row per case × method variant × phasing; absolute and `_block` columns), `allele_segments.tsv` (per-interval t_A, t_B, c_A, c_B, c_U) |
