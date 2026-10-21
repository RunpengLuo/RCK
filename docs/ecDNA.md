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
```

Code lives in `rck/ecdna/`; the command-line entry points are `rck/ecdna_run.py`
(`rck-ecdna`) and `rck/ecdna/evaluate.py` (`rck-ecdna-eval`). RCK's own pipeline
(`rck`, `rck/core/*`) is unchanged.

## 1. Install

`rck-ecdna` needs Python 3.8 or newer, `numpy`, `pandas`, `networkx` (already RCK
dependencies) and `gurobipy`.

```bash
cd RCK
pip install -e .          # installs rck, rck-ecdna, rck-ecdna-eval
pip install gurobipy      # or: conda install -c gurobi gurobi
```

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

## 3. Inputs

| Option | Content |
|---|---|
| `--subgraph-dir` | `nodes.tsv`, `edges.tsv` (and optionally `phase_blocks.tsv`) of one amplicon |
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
`--include-low-weight` is passed to `rck-ecdna-eval`. Decompositions are ranked by residual, then by the number of
species that are not `low_weight`, then by total species and objective, so one molecule plus
small leftovers ranks above the same molecule split in two.
Each species has a `topology` (`circular` or `linear`). Each decomposition gets an amplicon
call from its significant species: the circular fraction is the copy-weighted share of
circular species with weight ≥ `--min-ecdna-copies` (2). The call is `ecDNA` when that share is
≥ `--ecdna-call-fraction` (0.8), `linear` when it is ≤ 1 − that fraction, `mixed` in between,
and `none` when there are no species.
Alternatives come from Gurobi's solution pool (`--alternatives-mode pool`) or from repeated
solves with a no-good cut (`--alternatives-mode nogood`), within `--pool-gap` of the
optimum. Every alternative is re-checked and deduplicated.

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
| `decompositions.tsv` | one row per decomposition, species and ordering: `amplicon_call`, `circular_fraction`, `topology`, `status` (`species` or `low_weight`), `weight` (copies per cell, bulk average), `length_bp`, `path_support`, `uses_imputed`, `junctions`, `walk` (`chrom:start-end:hap:strand`, comma-separated) and `joins` (adjacency after each step, `kind:aid:phasing`; a linear walk's last join is `end:<extremity>:<hap>`) |
| `flow_segments.tsv` | per segment: observed total and allele copy numbers, plus inferred ecDNA and background copy numbers for A and B |
| `flow_adjacencies.tsv` | per adjacency labeling: read support, observed copy number, ecDNA and background flow; `linear_end` rows for copy-ends |
| `aa/<name>_graph.txt`, `aa/<name>_decomp<k>_cycles.txt` | AmpliconArchitect-style graph and cycles (total copy number; cycle lines also carry `Haplotypes=`; linear walks start and end with AA's `0+`/`0-` source) |
| `observations.json` | scales, how each was estimated, and every observation |
| `summary.json` | graph normalization log, flow status (with `linear_end_copies`), top-level `amplicon_call` and `circular_fraction`, decompositions |
| `gurobi_flow.log`, `gurobi_decompose.log` | solver logs |
| `evaluation.tsv`, `evaluation.json`, `truth_cycles.tsv` | only with `--truth-dir` or `--truth-cycles` |

## 6. Evaluation

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
2. **Haplotypes:** from the species' junction `PHASING` (`truth/adjacency_clone.tsv` +
   `obs/adjacency.tsv`), carried across reference joins. Steps no junction labels take the
   species' majority haplotype.
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
  ```

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
| `test_evaluate.py` | no | Jaccard with multiplicity; junction scores; order match with rotation, reversal, finer segmentation and tolerance; a linear inferred walk never matches a circular truth (`topology_match`); truth table I/O; rebuilding the demo3 truth from the simulator (16 steps, 2,343,000 bp, seg4 on A and B, seg10 on B); demo4 truth from `truth/aa` (two species, per-sample weights); scoping species to a subgraph |
| `test_ilp_circular.py` | yes (small) | single circle (10 copies on A plus a 1+1 background); haplotype crossings; imputed closing edge used when cheap, linear flow when expensive, no flow with `--linear-ends none`; linear ends only at boundaries and flanks; a linear amplicon; `--enforce-is` |
| `test_decompose.py` | yes (small) | repeated segment (one species, multiplicity 2); two species sharing a segment, with every optimal answer reported; disconnected circles; pool and no-good alternatives; orderings; a linear species; a circle plus a linear piece (called ecDNA); end-to-end `rck-ecdna` for a circle (`ecDNA`) and for an amplicon with no closing junction (`linear`, and `none` with `--linear-ends none`) |
| `test_gurobi_env.py` | no | thread default (explicit, `SLURM_CPUS_PER_TASK`, CPU affinity); gurobipy logger kept off the console |
| `test_demo4.py` | yes (full license) | both demo4 subgraphs: evaluation never errors, each subgraph is scored only against its own species, per-sample truth weights; accuracy is printed, not asserted (`RCK_ECDNA_DEMO4_OBS` enables estimated scales for both samples) |
| `test_demo3.py` | yes (full license) | `amplicon_0` against the demo3 truth, rank-1 decomposition only: cycle order matches, haplotype Jaccard ≥ 0.95, junction recall ≥ 0.9, copy-number error ≤ 30% |
