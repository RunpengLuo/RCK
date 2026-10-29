# Allele-specific ecDNA reconstruction (`rck-ecdna`)

`rck-ecdna` extends RCK from linear, clone-specific karyotypes to circular, haplotype-specific
ecDNA. It reads an ecDNA breakpoint subgraph built by
amplicon-allele-analysis `subgraph_construct.py` and reports every circular species the
data support, each as an ordered walk over phased segments, with alternatives.

```
subgraph (nodes.tsv, edges.tsv) ──► normalized DIAG ──► copy-number observations
                                                       │
                     two-layer flow MILP (ecDNA + background) ◄┘
                                                       │
            minimum-species cycle decomposition + alternatives
                                                       │
         alternating Euler walks + alternative orderings ranked by long-read paths
                                                       │
       recovery matrix against the truth (optional) + one figure page per species
```

Code lives in `rck/ecdna/`; the command-line entry points are `rck/ecdna_run.py`
(`rck-ecdna`), `rck/ecdna/evaluate.py` (`rck-ecdna-eval`) and `rck/ecdna/plot.py`
(`rck-ecdna-plot`). RCK's own pipeline (`rck`, `rck/core/*`) is unchanged.

## 1. Install

`rck-ecdna` needs Python 3.8 or newer, `numpy`, `pandas`, `networkx` (already RCK
dependencies) and `gurobipy`.

```bash
cd RCK
pip install -e .          # installs rck, rck-ecdna, rck-ecdna-eval
pip install gurobipy      # or: conda install -c gurobi gurobi
```

The figures (section 5) also need `matplotlib` and `cnplot`, the copy-number plotting package
in amplicon-allele-analysis `third_party/cnplot` (Python 3.10 or newer). Both are optional:
without them the run finishes and logs that the plots were skipped.

```bash
pip install -e /path/to/amplicon-allele-analysis/third_party/cnplot   # also installs matplotlib and seaborn
```

Instead of installing cnplot you can point `--cnplot-path` (or `$CNPLOT_PATH`) at that checkout.

`pip install gurobipy` ships a size-limited license (2,000 variables and 2,000 constraints).
It is enough for the unit tests but not for the demo3 decomposition, which needs your
academic license (see section 7).

## 2. Run on the demo3 simulation

From the root of amplicon-allele-analysis, with RCK installed:

```bash
rck-ecdna \
    --subgraph-dir subgraphs/demo3/amplicon_0 \
    --obs-dir simulations/runs/demo3_clone1_ecdna/obs \
    --truth-dir simulations/runs/demo3_clone1_ecdna \
    --outdir rck_ecdna_out/demo3/amplicon_0 \
    --threads 4
```

`--obs-dir` points at `windows/`, `extremity.tsv` and `path.tsv` in one go. Leave out
`--truth-dir` on real data. To skip the scale estimation, pass
`--depth-per-copy 15 --junction-reads-per-copy 6` (the demo3 values).

### Run from the subgraph folder only

Subgraphs built by amplicon-allele-analysis `simulations/subgraph_construct.py` (current
`yubi_dev`) carry everything RCK needs, computed from the observations without truth
(`simulations/rck_inputs.py`):

| File in the subgraph folder | Content |
|---|---|
| `scale.json` | depth per copy and junction reads per copy, the purity / tumor ploidy solution they come from, every ploidy class considered, the sample, and the sample the junction reads come from |
| `edges.tsv` columns `CN_TOTAL`, `CN_A`, `CN_B` | cell-averaged segment copy number and its allele split |
| `edges.tsv` column `CN` | junction copy number (`#REF` or `#NOV` / junction reads per copy) |
| `cn_profile.tsv` | copy number in 100 kb bins along each chromosome of the subgraph, for the region flanks (background layer) |
| `paths.tsv` | long-read paths whose extremities are all in the subgraph |
| `phase_blocks.tsv` | one row per phase block of the subgraph's segments and junction flanks: `votes_same`, `votes_flip`, `orientation` (1 as observed, −1 flipped, 0 not oriented), `n_snps`, `status` (`voted` / `no_votes`). Older folders list only the voted blocks |
| `segment_blocks.tsv` | one row per (segment edge, phase block it overlaps): the parent edge's `NODE_1 NODE_2 CHROM_1 POS_1 POS_2` and `edge_index`, `phase_block`, `piece_start`, `piece_end` (half-open; the pieces tile the edge), `piece_length`, `n_snps`, `A_count_obs` / `B_count_obs` (observed phase), `orientation`, `A_count` / `B_count` / `BAF` and `CN_TOTAL` / `CN_A` / `CN_B` (oriented; `CN_A` / `CN_B` empty below 5 SNPs). Read with `--block-flip` (section 4) |

```bash
python -m rck.ecdna_run --subgraph-dir subgraphs/demo4/amplicon_0 --outdir rck_out/amplicon_0 \
    --truth-dir simulations/runs/demo4_clone3_ecdna        # optional, for scoring only
```

`--inputs auto` (default) uses these files whenever `scale.json` exists; `--inputs windows`
returns to estimating from `--windows-dir` / `--obs-dir`. `--depth-per-copy`,
`--junction-reads-per-copy` and `--baseline-cn` still override. A subgraph describes one
sample (`subgraph_construct.py --sample`, default the first tumor); build one subgraph folder
per sample. `--sample` is ignored with a warning when it differs from `scale.json`.

**How the depth scale is fitted (in `rck_inputs.py`).** Depth only gives copy number up to a
scale, and dividing by 2 is wrong for any non-diploid tumor: a whole-genome doubling made the
benchmark copy numbers about 0.55x the truth. The fit follows ASCAT / HATCHet: 1 Mb blocks of
depth and allelic imbalance (phased SNP counts in 100 kb sub-blocks, binomial variance
removed), segmented per chromosome, then a grid over purity and the copy number at median
depth that puts every segment's tumor allele copies closest to whole numbers, plus a BAF
term that checks each segment's predicted imbalance (copy-neutral LOH reads BAF 0 in a pure
diploid tumor, (1 - purity) / 2 after a doubling). The best solution per ploidy class (about
2, 4, 6 copies) is kept and the class is chosen jointly over the tumor samples, because a
doubling is truncal. On the runs checked: demo3 13.76 (truth about 13.8, purity 0.80
recovered), demo4 13.68 / 13.62 (truth about 13.7), a run with a whole-genome doubling 7.29 /
7.27 at purity 0.80 (truth 7.29, 0.80). `scale.json` lists every class's solution;
`subgraph_construct.py --rck-ploidy 2|4|6` or `--rck-depth-per-copy` override the call.

**Limitation (simulator).** `adjacency.tsv`, `extremity.tsv` and `path.tsv` are written for
the first tumor sample only, so junction copy numbers always describe that sample
(`scale.json` `junction_sample`). For other samples the segment copy numbers are theirs but
the junction evidence is not.

## 3. Inputs

| Option | Content |
|---|---|
| `--subgraph-dir` | `nodes.tsv`, `edges.tsv` (and optionally `phase_blocks.tsv`, `segment_blocks.tsv`) of one amplicon; with `scale.json` etc. (above) it is the only input needed |
| `--inputs` | `auto` (subgraph files when `scale.json` exists), `subgraph`, or `windows` |
| `--windows-dir` | UGP windows: `window.tsv.gz`, `window.depth.npz`, `sample_ids.tsv`, optional `snp.tsv.gz` + `snp.[AB]allele.npz` |
| `--extremities` | `extremity.tsv` with genome-wide `#REF`, used to estimate junction reads per copy |
| `--paths` | `path.tsv` long-read paths, used to rank alternative orderings |
| `--hatchet-seg` | HATCHet `best.seg.ucn`, used to estimate depth per copy when there is no matched normal (cell lines such as COLO320) |
| `--truth-dir` / `--truth-cycles` | ground truth for the optional evaluation |

### Graph normalization

`subgraph_construct.py` edge types map onto RCK's interval adjacency graph as follows:
`reference` → segment, `extremity` → reference adjacency (`#REF`), `variant` → observed
novel adjacency (`#NOV`), `imputed_variant` → candidate novel adjacency. When the graph is
read, three things happen (`rck/ecdna/subgraph_io.py`):

1. **Phantom segments are contracted.** Simulated breakpoints have two extremities,
   `x+` and `x+1-`. Between them sits a 1 bp segment (`seg_hap = no_span`) that only ever
   carries reference flow. It is removed and its two reference joins are merged into one.
   Any other edge attached to the phantom moves to the true extremity on the same side.
   `--phantom-max-bp 0` turns this off.
2. **Parallel imputed edges are merged.** An imputed edge that duplicates a reference join
   (`inter` edges at imputed extremities, and anchors that loop back onto their own
   extremity) is merged into that reference adjacency.
3. **Terminals are marked.** These are the extremities with no segment (the region flanks)
   and chromosome-boundary nodes. Only the background layer may start or end at a terminal.

**Haplotype labelings.** Each adjacency gets a set of allowed RCK `Phasing` labelings:

- reference adjacencies: AA or BB;
- observed junctions follow `junction_hap`: `A` → AA, `B` → BB, `mixed` → AB or BA
  (a haplotype crossing). With `--soft-junction-hap` these become penalties instead of
  constraints;
- imputed candidates: AA or BB;
- novel self-loops can never be AB or BA, as in RCK.

### Copy-number scale

Two scales convert the data to copy numbers (`rck/ecdna/scale.py`):

- **depth per copy** (γ): segment copy number = mean tumor depth / γ;
- **junction reads per copy** (θ): junction copy number = reads / θ.

**γ** is taken from the first available source:

1. `--depth-per-copy`;
2. a HATCHet fit (`--hatchet-seg`): the length-weighted median of tumor depth ÷ fractional
   copy number over long background segments;
3. the windows: the mode of tumor depth over allelically balanced background blocks (mirrored
   BAF ≥ 0.4, away from the amplicon), divided by `--baseline-cn` (default 2).

`--depth-mode ratio` divides tumor depth by matched-normal depth window by window first.
That removes coverage bias the two libraries share (GC content, mappability), which is the
usual case for real data. The amplicon-allele-analysis simulator draws the bias field
independently for each sample, so on simulated data the ratio adds noise; that is why
`tumor` is the default.

**θ** comes from `--junction-reads-per-copy`, from the `#REF` of background extremities
(when `--extremities` and windows are given), or, as a logged fallback, from the median
`#REF` of the subgraph's own reference joins.

**Allele-specific observations.** These use the phase-oriented segment BAF from the
subgraph: c_A = c(1 − BAF), c_B = c × BAF. Segments with fewer than `--min-snps-allele`
SNPs get a total-only observation.

**Per-phase-block observations (`--block-flip unoriented`).** The segment BAF is only on the
true haplotype scale where `subgraph_construct.py` could orient every phase block the segment
touches. With `--block-flip unoriented`, every segment with rows in `segment_blocks.tsv` gets
one observation per *piece* (its part inside one phase block) instead:

- pieces are matched to segments by their parent edge (`NODE_1 NODE_2 CHROM_1 POS_1 POS_2`),
  following phantom contraction (a contracted phantom's pieces go to the segment that absorbed
  its base);
- the piece's total is its `CN_TOTAL`; its A / B are `CN_A` / `CN_B` (or, in windows mode or
  with `--depth-per-copy`, the piece total split by its BAF). A piece with fewer than
  `--min-snps-allele` SNPs is total-only;
- a piece in a block with orientation 1 or −1 in `phase_blocks.tsv` is *fixed* and taken in the
  oriented frame. A piece in a block with orientation 0, or a block missing from
  `phase_blocks.tsv`, is taken in the observed frame and depends on that block's flip (below).

Without `segment_blocks.tsv` the run logs a warning and uses the segment observations, as with
`--block-flip off` (the default).

## 4. Model

### Flow MILP (`rck/ecdna/ilp_circular.py`)

`OptModelCircular` keeps from RCK's `OptModelMultiClone`:

- integer copy numbers per segment haplotype and per adjacency labeling;
- big-M presence indicators;
- at most one labeling per novel adjacency;
- AA/BB-only reference adjacencies;
- the haplotype-aware balance constraint, with self-loops counted twice.

It changes the following:

| | RCK | `rck-ecdna` |
|---|---|---|
| telomeres | allowed at reference / listed telomeres | ecDNA-layer copies may start or stop only at *linear ends* (`--linear-ends`), each copy-end costing `--linear-end-penalty`; every other extremity is strictly balanced |
| chromosomal copies | part of the one karyotype | separate background layer that may end only at terminals, cannot use imputed edges, and is softly tied to the flank copy number |
| novel adjacency use | at least a fraction (1 − `--run-nas-fp`) must be used | no requirement; each one used costs `--novel-penalty` |
| missing junctions | n/a | imputed candidates, each costing `--imputed-penalty` |
| objective | length-weighted copy-number distance | length-capped, sqrt-normalized allele-specific segment distance + junction read support (`--alpha-novel`, `--alpha-ref`) + penalties |
| infinite sites | on | optional, `--enforce-is` (chromothripsis-derived ecDNA violates it) |

#### Linear amplicons

Not every focal amplification is circular. AmpliconArchitect and CoRAL also report linear
(BFB, heavily rearranged or tandem-duplicated) amplicons and only call an amplicon ecDNA when
cycles explain most of the amplified copy number. The ecDNA layer therefore allows a copy to
start or stop at selected extremities, using an integer end variable per extremity haplotype:
the balance becomes `segment flow = junction flow + ends`. `--linear-ends` picks the extremities:

- `boundaries` (default): an extremity with no reference neighbour, at a region flank, or where
  total copy number changes by at least max(`--boundary-min-copies`, `--boundary-frac` × the
  larger side). This is where the end of a linear amplicon is visible in the depth;
- `all`: every extremity of a segment;
- `none`: the strictly circular model of the previous versions.

Each copy-end costs `--linear-end-penalty` (0.25), so a linear walk costs 0.5 per copy. An
imputed closing edge costs `--imputed-penalty` (2.0) once, so with ≥ 4 copies an affordable
imputed closure still wins, while a long, otherwise unsupported gap is explained as linear.
A tiny `bg_balance_weight` term prefers a heterozygous background (1A + 1B); without it,
linear copies can swap haplotype at no cost.

#### Phase-block flips (`--block-flip unoriented`)

`subgraph_construct.py` orients a phase block only from votes of pure-haplotype junctions at
segments that are clearly one haplotype. A block with no votes keeps an arbitrary A/B labelling,
so its allele observations can contradict the junction haplotype labels (demo7: blocks chr4:185
and chr4:186, flow objective 39.3 with two imputed edges). With `--block-flip unoriented` the
flow MILP chooses the labelling:

- one binary f<sub>b</sub> per phase block the subgraph left unoriented that has an
  allele-specific piece, shared by every piece in that block, across segments;
- a piece with observed copy numbers a, b gets the targets t<sub>A</sub> = a + f<sub>b</sub>(b − a)
  and t<sub>B</sub> = b + f<sub>b</sub>(a − b). These are linear in f<sub>b</sub>, so the usual
  absolute-deviation constraints apply. The model value is the segment's inferred A / B copy
  number (ecDNA + background), which all pieces of a segment share; total-only pieces fit the
  segment total;
- the segment's data weight is split over its pieces in proportion to piece length (the
  segment weight is length-based), so each segment's total weight and the objective scale are
  unchanged. The 1/√c normalizer uses the segment's total copy number, so it is the same for A
  and B and does not depend on f<sub>b</sub>;
- each flipped block costs 10<sup>−6</sup>, so a block with no evidence either way stays as
  observed. The junction haplotype labels, through the balance constraints, then decide the
  orientation;
- after solving, `delta_if_reversed` is the objective change from reversing a block's flip with
  every copy number fixed. A block not oriented by the subgraph whose change is at most
  `--flip-tol` (0.01) is `unresolved` and kept as observed;
- `--flip-oriented-penalty P` (off by default) also gives the oriented blocks a flip, at cost P
  each, for when a block's votes are wrong.

The flips are fixed after the flow MILP: the copy-number refit, `flow_segments.tsv` and the
figures use the observations with the chosen flips applied.

On the demo7 union subgraph rebuilt with `segment_blocks.tsv` (amplicon-allele-analysis
`simulations/runs/demo7_nestedcycles_ecdna/rck_out_union_blockflip/`):

- the flow flips chr4:185 and chr4:186 and keeps chr4:171, chr4:172 and chr4:187 as oriented.
  The 12 blocks of the unamplified 172.4–185.5 Mb gap are `unresolved`;
- the flow objective is 9.13 with no imputed edge (8.8 with total-only observations, 39.3
  without flips);
- the rank-1 decomposition has the two nested species as separate circles, at 16.43 and 9.08
  copies against a true 16.90 and 9.10, with composites 0.994 and 0.9996.

Without flips, the total-only run merged the two cycles into one walk. Correctly phased, the
species share 9 reference segments but only one segment haplotype (`ecdna1` is on B and
`ecdna2` on A where they overlap), so a merged walk can no longer reproduce the
haplotype-specific flow. Nested species on the *same* haplotype may still merge; see the known
limitations in `CHANGELOG.md`.

### Cycle decomposition (`rck/ecdna/decompose.py`)

The ecDNA-layer flow is split into at most `--max-species` slots. Each slot has an integer
copy weight w and a per-edge multiplicity m between 0 and `--max-multiplicity`. Constraints:

- the products w × m reproduce the flow; any residual is penalized;
- each slot is balanced with alternation;
- each slot is connected (enforced with lazy cuts), so it forms a single closed walk;
- each slot uses either no linear end (a circle) or exactly two copy-ends (one linear walk).

The objective minimizes the residual plus `--species-penalty` × the number of species. A tiny
penalty on total edge traversals breaks ties in favour of the simplest walk, and a slot that
goes round the same circle g times is reported as that circle at g times the weight. Species
lighter than `--min-species-weight` (default `--min-species-weight-frac` = 10% of the
heaviest species) are kept but flagged `low_weight`: they usually absorb small unevenness in
the flow rather than describing a molecule, and the evaluation ignores them unless
`--include-low-weight` is passed to `rck-ecdna-eval`. Decompositions are ranked by residual, then by a structure
cost (`--species-penalty` per species that is not `low_weight`, plus `--step-penalty` per step
of those species, below), then by the number of such species, total species and objective. One
molecule plus small leftovers therefore ranks above the same molecule split in two.
Each species has a `topology` (`circular` or `linear`). Each decomposition gets an amplicon
call from its significant species: the circular fraction is the copy-weighted share of
circular species with weight ≥ `--min-ecdna-copies` (2). The call is `ecDNA` when that share is
≥ `--ecdna-call-fraction` (0.8), `linear` when it is ≤ 1 − that fraction, `mixed` in between,
and `none` when there are no species.
#### Walk complexity (`--step-penalty`) and ambiguity

Species that share reference sequence, such as nested ecDNA, can often be explained in several
ways that reproduce the flow, and so the copy numbers, exactly. Take demo8, an outer cycle at
19.5 copies with a nested cycle at 6.5 on the same haplotype. "Outer circle at 18 + nested
circle at 6" and "outer circle at 12 + one walk through outer and nested at 6" predict the same
copy number on every segment and junction: outer-only segments 18 = 12 + 6, shared segments
18 + 6 = 12 + 2 × 6. Neither a better fit nor more samples can separate them (unless the nested
species outnumbers the outer one somewhere); only reads that span the
difference, or a prior, can. demo7 with total-only observations shows the other form: one walk
round the outer cycle twice and the nested one once needs one species fewer, and the species
penalty chose it.

`--step-penalty λ` (default 0) adds λ per *step* of a significant species to the ranking's
structure cost. A step is one segment traversal of one species (a walk's length in segments, not
weighted by copies). This is a minimum-description-length prior:

- splitting a species into pieces keeps the total step count about the same, so it does not
  push towards splits;
- a merged walk re-traverses the shared sequence, so it costs the extra steps (demo8: 47 steps
  against 28 for the separate circles);
- merging saves one species and costs the extra steps, so it is undone once λ × extra steps
  exceeds `--species-penalty`;
- a true internal duplication (a tandem duplication traverses its segment twice) has the same
  steps as its split into two species, so the species penalty keeps it whole.

**Where the penalty acts.** It is applied in the ranking, to species that are not `low_weight`,
and not in the decomposition MILP. In the MILP it would also charge the light leftover species
that absorb unevenness in the integer flow. Their walks are often long (45 steps on demo3), so
the solver reshaped real species to shorten them: it split demo3's duplicated section and
demo7's tandem duplication into separate species. A weight threshold inside the MILP was gamed
too: it split a merged walk into copies too light to pay. With λ > 0, the solution pool keeps
every solution within `--ambiguity-window` of the best objective (never less than `--pool-gap`)
instead of within 10%. The alternatives the prior chooses between, which can need more species,
then reach the ranking. Solutions the prior cannot prefer are dropped again, so the listed
alternatives stay comparable. The wider pool costs time: the demo8 run went from 8 s to 31 s.

With λ = 0 only the 10<sup>−4</sup> tie-break on edge traversals separates such answers, as on
demo8, where the right one won by 0.004.

**Effect on the demos** (rank-1 composite per truth species; λ 0.05, 0.1 and 0.2 gave identical
results):

| run | λ = 0 | λ = 0.02 | λ ≥ 0.05 |
|---|---|---|---|
| demo3 | 0.960 | 0.959 | 0.959 (weight 12.77 → 12.70, an equivalent solution) |
| demo4 amplicon_0 / amplicon_1 | 0.983 / 0.997 | same | same |
| demo5 | 0.990 | same | same |
| demo6 union | 0.984, 0.994 | same | same |
| demo7 union, `--block-flip unoriented` | 0.994, 1.000 | same | same |
| demo7 union, total-only (`--min-snps-allele 100000000`) | 0.51, 0.34 (merged walk) | same | **0.98, 0.97** (two species, no leftover) |
| demo8 union (3:1, one haplotype) | 0.985, 0.989 (tie-break, gap 0.004) | same | same, decided by `step_penalty` (gap 0.95 at 0.05, 1.9 at 0.1) |

demo7's total-only merge needs λ above about 1/36: merging saves one species and adds 36 steps
(67 against 31). The default stays 0, so earlier results are unchanged; `--step-penalty 0.1`
is recommended for nested amplicons.

**Ambiguity report** (`--ambiguity-window`, default 2.0; 0 turns it off). After ranking, one more
MILP looks for alternatives that the data cannot tell apart from the top decomposition. It fixes
the ecDNA flow covered by every edge to the top decomposition's, so the fit is identical. It
forbids two species with the same support, which excludes one walk split into copies. Then it
enumerates solutions within the window of the objective. An alternative whose significant species
structures differ makes the top decomposition `ambiguous`. `summary.json` reports, under
`ambiguity` (top level, and in the first decomposition):

- `ambiguous` and `n_alternatives` (distinct structure sets);
- `gap`: the closest alternative's objective minus the top one's. Negative means the ranking
  preferred fewer species over a lower objective;
- `decided_by`: what separates them: `species_penalty`, `step_penalty`, `tie-break`, or
  `reads or tie-break`;
- the first 10 alternatives: species count, steps, integer weights, and their rank if listed.

A run whose top decomposition is ambiguous also logs a warning. Ambiguity is expected with
nested species. The report says whether the choice came from the data or from a prior, and how
strongly the prior decided. The extra MILP has the size of the decomposition MILP plus K(K−1)/2
binaries per edge (K = `--max-species`).

Alternatives come from Gurobi's solution pool (`--alternatives-mode pool`) or from repeated
solves with a no-good cut (`--alternatives-mode nogood`), within `--pool-gap` of the
optimum. Every alternative is re-checked and deduplicated.

**Re-ranking a kept pool** (`--keep-pool`, `rck/ecdna/rerank.py`). λ acts only in the ranking,
and every λ > 0 widens the pool to the same `--ambiguity-window`, so one run at any λ > 0 holds
every answer a λ > 0 can give. With `--keep-pool` (it needs `--step-penalty` > 0) the run also
writes `pool/`: every pool solution before the λ-dependent pruning, annotated, refitted and, with
a truth, evaluated like a ranked decomposition; and, for each solution that is rank 1 at some λ of
`--pool-lambdas` (default 0 0.02 0.05 0.1 0.2 0.5), its same-fit alternatives (the ambiguity
MILP does not depend on λ either). The run's own outputs are unchanged.

```bash
rck-ecdna --subgraph-dir SUB --truth-dir RUN --outdir OUT --step-penalty 0.1 --keep-pool
python -m rck.ecdna.rerank OUT --lambdas 0.02 0.05 0.1 0.2 0.5    # one row per λ
```

`rerank` applies the run's own pruning (`pool_keeps`) and ranking (`ranking_key`) to the recorded
pool and gives, per λ, the rank-1 solution, its refitted significant species, steps, residual,
composite per truth species and `ambiguous`, `decided_by`, `gap`, `n_alternatives` (missing for a
λ whose rank-1 solution was not rank 1 at any `--pool-lambdas` value). λ = 0 is only emulated (the
solutions within `--pool-gap` of the best): a λ = 0 run solves a narrower pool, whose ties can
come out differently (demo3: 12.77 copies at λ = 0, an equivalent 12.70 at λ > 0), so run it
instead. `results/bench_nested/` uses this to sweep λ with two runs per case.

#### Noise-tolerant residual (opt-in)

The flow MILP rounds noisy observations to an integer flow, and the decomposition then covers
that flow. By default one copy of residual costs 1 on a segment of 50 kb or more (less on shorter
ones, down to 0.05) and `--adj-residual-weight` (1.0) on a junction, the same as one species.
A dip of 1–3 copies from depth bias is therefore cheaper to cover exactly with a leftover or a
merged species than to leave as residual. Four options, all off by default, change this:

- `--seg-residual-scale S` multiplies the segment weights by S;
- `--adj-residual-weight W` sets the junction weight;
- `--residual-mode relative` divides every edge's weight, junctions included, by the expected noise
  `--residual-sigma` (0.1) × max(flow, 1), so the residual counts noise units. Junctions are
  included because the flow is balanced: a dip on a segment comes with dips on its junctions.
  Below 1 / σ copies (10) the relative weight is above the absolute one;
- `--rank-residual sum` ranks decompositions on residual + structure cost instead of residual
  first. Lower weights alone only change what reaches the pool (and what `--pool-gap` keeps);
  with the default ranking, an exact fit in the pool still ranks first.

For example, `--adj-residual-weight 0.25 --seg-residual-scale 0.25 --rank-residual sum` lets a
residual of up to 4 copies on long segments or junctions cost less than one species. `--keep-pool`
records these settings, and `rck.ecdna.rerank` ranks with the run's `--rank-residual`.

#### Long reads in the decomposition

With `--paths`, every long-read path that spans at least `--min-path-junctions` (2) junctions,
has at least `--min-path-support` reads, and whose junctions all carry ecDNA flow becomes a
read constraint (identical paths, in either direction, are merged and their reads summed).
For read p and slot k, a binary z<sub>kp</sub> is 1 exactly when slot k uses every junction of
p, and x<sub>p</sub> = OR<sub>k</sub> z<sub>kp</sub> marks p as explained. Two ways to use them, as in CoRAL:

- **reward** (default): unexplained reads cost `--path-weight` (0.5) × their share of the read
  support. The total is capped at 0.5, below one species (`--species-penalty` 1.0), so reads
  choose between decompositions that fit the flow equally well (e.g. two circles sharing a
  segment against one merged circle) but never add a species on their own;
- **requirement**: `--path-fraction 0.9` requires that share of the support to be explained
  (CoRAL's default). If no decomposition can, the run logs a warning and solves without it.

z and x are fully determined by the slot supports, so they do not multiply the solution pool.
Containing every junction is necessary, not sufficient: after solving, each read is checked
against the actual orderings, and `reads_explained` reports the share of support found in
order in some species' top ordering.

#### Fractional copy numbers and uneven flow (`rck/ecdna/refit.py`)

The structure search uses integer weights and reproduces the integer flow. That flow is
fitted edge by edge to noisy observations, so it can be uneven (15 copies on most edges of a
molecule, 28 on one segment). Matching it exactly then needs leftover species, and an integer
leftover can take copies from the real molecule (demo4 tumor2: 6 + 1 against a true 7.44).

After the decomposition, every decomposition keeps its structures and its copy numbers are
refitted as real numbers, directly against the observations, by a small linear program. It
uses the flow MILP's own data terms (allele-specific segment copy numbers and junction reads,
same weights), with the background layer held at the flow solution. A species refitted below
`--refit-min-copies` (0.5) is not supported by the data: the lightest one is dropped and the
fit repeated. The ranking is unchanged (it is made on the integer solution); `low_weight`
flags and the amplicon call are recomputed from the new weights. `--no-refit` reports the
integer weights. The LP is solved with scipy (HiGHS), so it needs no Gurobi licence.

### Walks (`rck/ecdna/walk.py`)

For each species, Hierholzer's algorithm with alternation gives one walk, and a bounded
depth-first search enumerates up to `--max-orderings` distinct orderings, up to rotation and
reversal (only reversal for a linear walk, which starts at one end and finishes at the other).
Orderings of the same multiset have identical copy numbers. They are ranked by the
support of long-read paths (`--paths`) that span two or more junctions and appear in that
order.

## 5. Outputs

| File | Content |
|---|---|
| `decompositions.tsv` | one row per decomposition, species and ordering: `amplicon_call`, `circular_fraction`, `reads_explained`, `data_objective` (observation misfit after the refit), `topology`, `status` (`species` or `low_weight`), `weight` (copies per cell, bulk average, fractional after the refit), `weight_int` (integer weight of the structure search), `length_bp`, `path_support`, `uses_imputed`, `junctions`, `walk` (`chrom:start-end:hap:strand`, comma-separated) and `joins` (adjacency after each step, `kind:aid:phasing`; a linear walk's last join is `end:<extremity>:<hap>`) |
| `flow_segments.tsv` | per segment: observed total and allele copy numbers, plus inferred ecDNA and background copy numbers for A and B. With `--block-flip`, a segment with pieces reports its total split by the BAF of its piece counts summed after the chosen flips |
| `flow_adjacencies.tsv` | per adjacency labeling: read support, observed copy number, ecDNA and background flow; `linear_end` rows for copy-ends |
| `aa/<name>_graph.txt`, `aa/<name>_decomp<k>_cycles.txt` | AmpliconArchitect-style graph and cycles (total copy number; cycle lines also carry `Haplotypes=`; linear walks start and end with AA's `0+`/`0-` source) |
| `observations.json` | scales, how each was estimated, and every observation; with `--block-flip`, also `pieces`: each piece's segment, block, span, totals, A / B in its base frame, `fixed`, `flipped` and the A / B after the flip (`a_final`, `b_final`) |
| `block_flips.tsv` | with `--block-flip` only: one row per phase block with a piece: `phase_block`, `orientation` (from the subgraph), `flip` (1 when A and B were swapped), `status` (`fixed`, `flipped`, `kept` or `unresolved`), `n_pieces`, `n_snps` and `delta_if_reversed` |
| `summary.json` | graph normalization log, flow status (with `linear_end_copies`), top-level `amplicon_call` and `circular_fraction`, decompositions (with `reads_total`, `reads_explained`, `data_objective` and `data_objective_int`, `pruned`: species the refit dropped, and `steps`), top-level `ambiguity` (section 4, walk complexity); with `--block-flip`, `block_flip` (mode, number of pieces and blocks, counts per status, the flipped blocks) |
| `pool/pool.tsv`, `pool/pool.json`, `pool/decompositions.tsv`, `pool/evaluation.tsv` | with `--keep-pool` only: every solution of the decomposition pool before the λ-dependent pruning (`pool_index`, objective, residual, species, significant species and their steps, integer weights), the same solutions as decompositions (`decomposition_rank` = pool index + 1) and their evaluation; `pool.json` adds the settings, the rank-1 solution per `--pool-lambdas` value and those solutions' same-fit alternatives (section 4, re-ranking) |
| `gurobi_flow.log`, `gurobi_decompose.log` | solver logs |
| `evaluation.tsv`, `evaluation.json`, `truth_cycles.tsv`, `recovery_matrix.tsv` | only with `--truth-dir` or `--truth-cycles` (section 6) |
| `plots/<name>_ecdna.pdf`, `plots/<name>_rank<r>_<species>.png`, `plots/<name>_recovery_matrix.png` | figures: one page per inferred species, then the recovery matrix (below) |

### Figures (`rck/ecdna/plot.py`)

At the end of a run, `rck-ecdna` draws one page for each species of the top decomposition. Species
flagged `low_weight` are left out unless you pass `--plot-low-weight`, and `--plot-ranks N` also
draws the next N − 1 alternatives. Every track shares one cnplot `GenomeAxis` that covers only
the segments the species and its matched truth species use, padded by 5% (at least 5 kb), with the
gaps between them collapsed. The tracks are:

- **topology**: one row per step of the top ordering, numbered. Each step is a block at its
  reference position, blue for haplotype A and orange for B, with an arrow for the strand.
  Connectors follow the traversal order and are dotted for an imputed join and for the join that
  closes the circle. This is CoRAL's cycles figure, the same track amplicon-allele-analysis draws
  for simulated truth (`simulations/src/plots.py`), so the two can be compared side by side;
- **truth**: the matched truth species drawn the same way, when the run was evaluated;
- **junctions**: one arc per non-reference join of the walk, coloured by the haplotype it joins.
  Grey marks a haplotype crossing (AB / BA) and a dotted arc an imputed junction;
- **CN-A / CN-B**: allele-specific copy number, B mirrored under A. Points are the observations,
  drawn with cnplot `plot_scatter_1d` along each segment (with `--block-flip`, along each piece,
  after the chosen flips). The step line is the inferred ecDNA +
  background copy number. Shaded blocks are this species' own copies (weight × the number of
  times its walk uses the segment).

Beside the tracks, the walk is drawn as a ring: clockwise from the top in traversal order, one arc
per step coloured by haplotype, with a tick at every novel junction (dotted for an imputed one). A
linear species is drawn as an open ring. Arc lengths are to scale, except that every step gets at
least 1.5% of the ring, which the page notes when it applies. The title gives the decomposition's
call, the species' weight (refitted and integer), its length and, when the run was evaluated, the
matched truth species and its scores. After an evaluation, the last page is a heatmap of the
recovery matrix.

`rck-ecdna-plot --outdir <run>` redraws an existing output folder (`--plot-ranks`,
`--plot-low-weight`, `--cnplot-path`, `--no-png`), and `--no-plots` turns the figures off.
A plotting error is logged as a warning and never fails the run.

## 6. Evaluation

`evaluation.tsv` reports `cn_rel_err` for the refitted weight and `cn_rel_err_int` for the
integer weight, so the effect of the refit is visible per species.

```bash
rck-ecdna-eval --decompositions rck_ecdna_out/demo3/amplicon_0/decompositions.tsv \
               --truth-dir simulations/runs/demo3_clone1_ecdna
```

or pass `--truth-dir` (and optionally `--truth-sample`) to `rck-ecdna` itself.

**How truth is read from a simulator run (`--truth-dir`).** For each species listed in
`truth/ecdna.tsv`:

1. **Walk and coordinates:** taken from `truth/aa/<species>_cycles.txt` when it exists (all
   current simulator runs). This works whatever the species went through (deletions, shared
   chromosomes, several species). Older runs without `truth/aa/` fall back to rebuilding the
   walk from `truth/walks.walk` + `truth/segments.tsv`.
2. **Haplotypes:** from the species' junction phasings (`truth/adjacency_clone.tsv` +
   `obs/adjacency.tsv`), carried across reference joins. Steps no junction labels take the
   species' majority haplotype. `PHASING` keeps one phasing per junction, its most common.
   Runs simulated with per-step haplotype labels (amplicon-allele-analysis
   `hap_labels: per_step`, every run before that key) can traverse one junction with several
   phasings, such as AA on one copy and BB on its duplicate; `truth/rck/rck.acnt.tsv` lists them
   all. Each traversal of such a junction takes one of them: the assignment with the fewest
   haplotype conflicts along the walk, then the one that uses the most of the junction's
   phasings. A crossing between reference neighbours may also be traversed on one haplotype,
   as a plain reference join that the truth does not list. Applying one `PHASING` to every
   traversal had mislabelled 0.995 Mb of the 5.37 Mb bench_nested `ctrl_dup_s702` truth walk.
3. **Copy number:** the species' copies per cell in the evaluated sample (`ecdna_burden` rows
   of `truth/proportions.tsv`). The sample is `--truth-sample`, else `--sample`
   (`sim_tumor1` maps to `tumor1`), else the first tumor. Without proportions it falls back to
   `ecdna.tsv` `PHI`.

The rebuilt truth is written to `truth_cycles.tsv`. For real data, write that table yourself
and pass it as `--truth-cycles` (columns `species order chrom start end strand hap weight`,
0-based inclusive coordinates).

**Linear inferred species.** Truth species are circular. A linear inferred walk is scored
without its wrap-around junction and never counts as an order match; `evaluation.tsv` has
`inferred_topology` and `topology_match` columns.

**Which truth species are scored.** A run covers one amplicon subgraph, but a simulation can
hold several species. Only truth species with at least `--min-truth-fraction` (default 0.5) of
their bases inside the subgraph's regions are scored. The others are listed in
`evaluation.json` under `truth_species` with `in_scope: false` and their
`fraction_in_subgraph`, and do not count as missed. `--all-truth-species` scores every species.
Subgraph regions split wherever segments on one chromosome are more than 1 Mb apart.

**Metrics.** Truth and inferred cycles are compared after merging reference-contiguous steps
into runs, so segmentation differences do not count against the inference:

- base-pair Jaccard, with and without haplotype;
- junction recall and precision, unphased and phased;
- order match up to rotation and reversal;
- relative copy-number error.

All coordinates are compared within `--eval-tol-bp` (default 10 kb; imputed breakpoints sit
on a 10 kb grid).

### Recovery matrix (`recovery_matrix.tsv`)

`recovery_matrix.tsv` has one row per truth species and shows how well the run recovered it. The
scores describe the rank-1 decomposition, which is RCK's answer; the `best_*` columns describe the
best alternative. Each truth species is matched one to one to an inferred species by haplotype
Jaccard, as in `evaluation.tsv`. The table is also stored under `recovery_matrix` in
`evaluation.json`, with its in-scope mean as `mean_composite`.

| Column | Meaning |
|---|---|
| `truth`, `in_scope`, `fraction_in_subgraph` | truth species and whether it is scored in this subgraph (out-of-scope rows have no scores) |
| `truth_length_bp`, `truth_n_junctions`, `truth_weight` | its size, its distinct junctions and its copies per cell |
| `inferred`, `inferred_topology`, `inferred_weight` | the matched rank-1 species (empty when it was missed) |
| `segment_recall`, `segment_precision`, `segment_jaccard` | bases covered by both (with multiplicity, haplotype ignored) over the truth bases, over the inferred bases, and over their union |
| `haplotype_jaccard` | the same Jaccard with haplotypes |
| `junction_recall`, `junction_precision`, `junction_f1` | distinct junctions, unphased |
| `phased_junction_recall`, `phased_junction_precision`, `phased_junction_f1` | distinct junctions with both haplotypes right |
| `phasing_accuracy` | of the bases both cover, the share on the right haplotype |
| `phasing_accuracy_swap_invariant` | the same, with A/B as given or swapped, whichever is better (for real data, where the phase labels need not match the truth's) |
| `order_match`, `topology_match` | as in `evaluation.tsv` |
| `cn_rel_err` | \|weight − truth weight\| / truth weight |
| `composite` | equal-weight mean of `segment_jaccard`, `junction_f1`, `phasing_accuracy`, `order_match` (0 or 1) and max(0, 1 − `cn_rel_err`); 0 for a missed species |
| `best_rank`, `best_inferred`, `best_composite` | the decomposition (and species) with the highest composite for this truth species, over every reported alternative: an upper bound on what the ranking could have picked |

`evaluation.tsv` carries the same per-pair scores (`segment_*`, `junction_f1`,
`phased_junction_f1`, `phasing_accuracy*`, `composite`) for every decomposition.

**Failures.** If the truth files cannot be read, `rck-ecdna` still finishes: the
reconstruction outputs are written, and `evaluation.json` records the error and traceback.

## 7. Running on Princeton ionic

ionic schedules jobs with Slurm. `docs/examples/ecdna_ionic.sbatch` is a template.

- **Gurobi license:** a named-user academic license is tied to the machine where you ran
  `grbgetkey`, so it does not work on compute nodes. For batch jobs, request a Gurobi
  academic WLS license and put its `gurobi.lic` (which holds `WLSACCESSID`, `WLSSECRET`,
  `LICENSEID`) in your home directory, or point `GRB_LICENSE_FILE` at it. If ionic provides
  a Gurobi module (`module avail gurobi`), its floating license works as well.
- **Environment:** create a conda environment once, on the login node:

  ```bash
  conda create -n rck-ecdna python=3.11 numpy pandas networkx scipy -y
  conda activate rck-ecdna
  pip install gurobipy
  pip install -e /path/to/RCK
  pip install -e /path/to/amplicon-allele-analysis/third_party/cnplot   # optional: figures
  ```

- **Batch jobs:** put the environment's `bin/` on `PATH`
  (`export PATH=$HOME/.conda/envs/<env>/bin:$PATH`), as `results/step_sweep/sweep.sbatch` does.
  `source ~/.bashrc; conda activate` can fail on a compute node, where conda's shell setup is
  not loaded. Keep array jobs to 2 concurrent tasks (`--array=0-N%2`) for the academic WLS licence.
- **Nested-ecDNA benchmark** (`results/bench_nested/`): 14 simulations from amplicon-allele-analysis
  `simulations/make_scenarios.py --set nested` and the 8 step-sweep cases, for choosing
  `--step-penalty`. Simulations and subgraphs run as one array without Gurobi
  (amplicon-allele-analysis `simulations/bench_nested.sbatch`), RCK as a second
  (`results/bench_nested/bench_nested.sbatch`), and `collect.py` summarizes; the order is in
  `results/bench_nested/README.md`.

## 8. Tests

The tests are plain `unittest` modules under `tests/ecdna/`. Run them from the RCK root,
either with the standard library or with pytest:

```bash
python -m unittest discover -s tests/ecdna -t . -v
# or
python -m pytest tests/ecdna -v
```

Gurobi uses `--threads`, else `SLURM_CPUS_PER_TASK`, else the CPUs the process may run on,
so it never starts more threads than the job was given. Its log goes to `gurobi_*.log`; pass
`-v` to also see it on the console.

Test outputs go to `$RCK_ECDNA_TEST_TMP` if set, else the system temp directory if it has at
least 200 MB free, else `~/.cache/rck-ecdna-tests`, and temp-directory cleanup errors (such as
`.nfsXXXX` placeholders on NFS homes) never fail a test. Every Gurobi model is disposed, and its
log file closed, as soon as its solution has been read. `docs/examples/ecdna_demo3_tests.sbatch`
runs the whole suite, both demo3 tests included, as a Slurm job on ionic.

Tests that solve a MILP are skipped when `gurobipy` cannot be imported. `test_demo3.py`
also needs a full Gurobi license. Setting `RCK_ECDNA_DEMO3_OBS` to
`.../simulations/runs/demo3_clone1_ecdna/obs` enables its second test, which estimates the
scales from the windows.

| File | Needs Gurobi | Covers |
|---|---|---|
| `test_subgraph_io.py` | no | synthetic chains; phantom contraction and edge relocation; parallel imputed merging; junction labelings; the `amplicon_0` fixture (18 segments; 20 reference, 14 novel and 12 imputed adjacencies; 4 terminals; crossings AB/BA) |
| `test_scale.py` | no | depth mode from balanced blocks (rejecting LOH and CN 3 blocks); allele split; flank background; override and fallback; ratio mode cancelling shared bias; θ from extremities; sample selection |
| `test_walk.py` | no | alternating Hierholzer; canonical form under rotation and reversal; a tandem-duplicated segment; a hub with three loops (2 orderings); balance checks; path parsing and ranking; linear walks from end to end (round trip, canonical form, one ordering, an unpaired end rejected) |
| `test_evaluate.py` | no | Jaccard with multiplicity; junction scores; order match with rotation, reversal, finer segmentation and tolerance; a linear inferred walk never matches a circular truth (`topology_match`); truth table I/O; rebuilding the demo3 truth from the simulator (16 steps, 2,343,000 bp, seg4 on A and B, seg10 on B); demo4 truth from `truth/aa` (two species, per-sample weights); scoping species to a subgraph; truth walks with a junction traversed as AA and BB, and with a crossing between reference neighbours traversed back on one haplotype, labelled without conflict; every phasing read from `rck.acnt.tsv` |
| `test_ilp_circular.py` | yes (small) | single circle (10 copies on A plus a 1+1 background); haplotype crossings; imputed closing edge used when cheap, linear flow when expensive, no flow with `--linear-ends none`; linear ends only at boundaries and flanks; a linear amplicon; `--enforce-is` |
| `test_decompose.py` | yes (small) | repeated segment (one species, multiplicity 2); two species sharing a segment, with every optimal answer reported; disconnected circles; pool and no-good alternatives; orderings; a linear species; a circle plus a linear piece (called ecDNA); end-to-end `rck-ecdna` for a circle (`ecDNA`) and for an amplicon with no closing junction (`linear`, and `none` with `--linear-ends none`); read constraints (merging a read and its reverse, dropping single-junction, weak and off-flow reads); a read through a>x>b choosing the merged circle (reward and `--path-fraction`); a read never adding a species, and the order post-check; noise tolerance: residual weights of `--seg-residual-scale`, `--adj-residual-weight` and `--residual-mode relative` (no Gurobi), the `--rank-residual` key, and a walk with a 2-copy dip on one loop, an exact two-species fit by default and one walk with cheaper residuals (absolute or relative) and the sum ranking |
| `test_refit.py` | LP: no; CLI: yes (small) | a fractional weight (10.4); an unsupported leftover dropped; a supported one kept; the `--refit-min-copies` threshold; the last species never dropped; observation rows use the flow background; end-to-end fractional copy number (10.4, integer 10) and `--no-refit` |
| `test_subgraph_inputs.py` | observations: no; run: yes (small) | scales, segment and junction copy numbers and flank background read from `scale.json`, `CN_*` columns and `cn_profile.tsv`; CN columns win over the scale; `--depth-per-copy` recomputes from `WEIGHT`; `--inputs auto/subgraph/windows`; a full run from the subgraph folder alone |
| `test_recovery.py` | no | segment recall and precision, junction F1, phasing accuracy (one haplotype wrong; every haplotype swapped, which only the swap-invariant score forgives), the composite and its copy-number term; the recovery matrix: a missed species scores 0, `best_rank` finds a better alternative, out-of-scope rows are listed but not scored or averaged, the TSV columns, strict JSON, and a run with no decompositions |
| `test_plot.py` | no (figures: matplotlib and cnplot) | axis regions padded and merged; figures skipped when cnplot cannot be imported; a circular species with truth and recovery-matrix pages and their PNGs; a linear species without truth |
| `test_step_penalty.py` | yes (small) | a hub segment with two 3-segment loops (nested-like). At 2:1 the merged walk wins without `--step-penalty` (ambiguous, decided by the species penalty) and the separate circles win with it. At 3:1 the answer is decided by the tie-break and, with `--step-penalty`, by the steps; the gaps match λ × extra steps. Default ranking unchanged; a single circle is not ambiguous; `--ambiguity-window 0` turns the report off. Pool ranking (no Gurobi): λ reorders and prunes recorded entries, residual first. `--keep-pool`: the pool kept at λ = 0.1 re-ranks to the rank-1 decomposition and ambiguity fields of direct runs at λ = 0.02–0.5 (2:1 and 3:1); the command line refuses it at λ = 0, writes `pool/`, and `rerank` reproduces the run's rank 1 |
| `test_block_flip.py` | inputs: no; flow and run: yes (small) | pieces read and mapped through phantom contraction; orientation −1 swapped and fixed; a block missing from `phase_blocks.tsv` is unoriented; low-SNP pieces total-only; a missing `segment_blocks.tsv` falls back with a warning; a circle on A whose unoriented block has swapped observations: the block is flipped, all ecDNA on A, no imputed edge, the objective of the unswapped input; `--block-flip off` unchanged; an oriented block never flipped (but flips with `--flip-oriented-penalty`); a block without evidence `unresolved`; a segment over an oriented and an unoriented block; the refit uses flipped observations; `block_flips.tsv`, `summary.json`, `observations.json` pieces, `flow_segments.tsv` and the plot points |
| `test_gurobi_env.py` | no | thread default (explicit, `SLURM_CPUS_PER_TASK`, CPU affinity); gurobipy logger kept off the console |
| `test_demo4.py` | yes (full license) | both demo4 subgraphs: evaluation never errors, each subgraph is scored only against its own species, per-sample truth weights; accuracy is printed, not asserted (`RCK_ECDNA_DEMO4_OBS` enables estimated scales for both samples) |
| `test_demo7.py` | yes (full license) | `--block-flip unoriented` on a rebuilt demo7 union subgraph (`RCK_ECDNA_DEMO7_SUBGRAPH`, with `segment_blocks.tsv`): chr4:185 and chr4:186 flipped, chr4:171, chr4:172 and chr4:187 fixed, no imputed edge; the flow objective is printed. Skipped without a full licence |
| `test_demo3.py` | yes (full license) | `amplicon_0` against the demo3 truth, rank-1 decomposition only: cycle order matches, haplotype Jaccard ≥ 0.95, junction recall ≥ 0.9, copy-number error ≤ 30%; the same with the fixture's long-read paths in the decomposition (≥ 90% of read support explained) |
