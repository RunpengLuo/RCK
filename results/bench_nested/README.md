# bench_nested: nested-ecDNA benchmark and the choice of `--step-penalty`

bench_nested is a small end-to-end benchmark of the pipeline's last two stages. Simulated tumours
with nested ecDNA species pass through the subgraph stage (amplicon-allele-analysis
`simulations/subgraph_construct.py` and `rck_inputs.py`) and then through `rck-ecdna`. Each
answer is scored against the simulator's truth. The benchmark has two uses:

1. **Choose the walk-complexity penalty λ (`--step-penalty`).** The rule is the λ selection rule
   below.
2. **Measure how well nested species and haplotype-crossing novel adjacencies (junctions) are
   recovered, end to end.** When a case fails, the failure-stage rule says which stage explains
   it.

The pipeline is described in [docs/pipeline/README.md](../../docs/pipeline/README.md):
simulation → subgraph → `rck-ecdna` → evaluation. The decomposition and its ranking (λ, the
solution pool, the ambiguity report) are in [docs/ecDNA.md](../../docs/ecDNA.md). Every metric
used here (the composite, its parts and the method-agnostic metrics) is defined in
[docs/pipeline/evaluation.md](../../docs/pipeline/evaluation.md). This README defines only the
rules that belong to the benchmark.

**Current state (2026-11-01).** The benchmark is the *grown* bench_nested_v2:

- 8 scenarios, 32 bench cases and 8 regression cases, 40 cases in all;
- inherited haplotype labels;
- subgraphs rebuilt with the long-read path fix at ploidy class 2.

At the chosen λ = 0.2 the mean composite is 0.8531, and 17 of the 32 bench cases fail. The plan
for these failures is in [NEXT_STEPS.md](NEXT_STEPS.md), "Grown bench: failure analysis and
revised plan (2026-11-01)". The earlier versions of the benchmark are kept under
[Earlier runs](#appendix-b-earlier-runs).

> **Provenance:** (c) new, RCK yubi_dev@949f796 (harness, `merged` flag), @2fdf73d (λ plateau
> rule and `scale` stage), @e4c2409 (`--aa-metrics`); scenarios in amplicon-allele-analysis
> yubi_dev@d1a1dce, @618697f (inherited labels) and @f1b09b6 (`nest_samehap_6to1`).
> **Code:** `results/bench_nested/run_case.py`, `results/bench_nested/collect.py`;
> amplicon-allele-analysis `simulations/make_scenarios.py`, `simulations/bench_nested.py`.

The design follows the simulation benchmarks of AmpliconArchitect [Deshpande et al. 2019](../../docs/pipeline/references.md#deshpande2019) and CoRAL
[Zhu et al. 2024](../../docs/pipeline/references.md#zhu2024): known structures are simulated, reconstructed and scored against the truth. What
is specific here is the focus on species that share reference sequence (nested species), which
neither of those benchmarks varies on purpose.

## 1. Cases

### 1.1 Bench scenarios

> **Provenance:** (c) new, amplicon-allele-analysis yubi_dev@d1a1dce, then @618697f and
> @f1b09b6. **Code:** `simulations/make_scenarios.py::NESTED_SCENARIOS`, `NESTED_INTENT`,
> `NESTED_SEEDS`.

Every scenario is demo7's simulation config with one or two keys changed. Common to all:

- 3 clones and 1 tumour sample;
- long-read counts, no simulated reads;
- library size, read model and bin size as in demo7, so junction reads per copy stay 6.5.

All scenarios set `hap_labels: inherit` with `p_cross` 0.5; the same-haplotype scenarios set
`p_cross` 0 and `nest_same_haplotype: true` (`NESTED_BASE`, `SAME_HAPLOTYPE`). The generators,
the inherited labels and the nesting options are described in amplicon-allele-analysis
[simulations/README.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/README.md).

| scenario | species | nesting | haplotypes | bulk ratio (`ecdna_proportions`) | intent (`NESTED_INTENT`) | purpose |
|---|---|---|---|---|---|---|
| `nest_cross_2to1` | 2 | nested | `p_cross` 0.5 | 2:1 | nested, cross | demo7-like with new seeds |
| `nest_samehap_2to1` | 2 | nested | one haplotype | 2:1 | nested, samehap | a merged walk saves a species |
| `nest_samehap_3to1` | 2 | nested | one haplotype | 3:1 | nested, samehap | demo8-like: the separate and merged answers tie on species count, so λ decides |
| `nest_cross_3sp` | 3 | 2 inside 1 | `p_cross` 0.5 | 0.5 / 0.3 / 0.2 | nested, cross | several nested species with crossings |
| `nest_cross_1to1` | 2 | nested | `p_cross` 0.5 | 0.55 / 0.45 (1.22:1) | nested, cross | close, non-integer copy numbers |
| `ctrl_separate` | 2 | none (`allow_interval_reuse` false) | `p_cross` 0.5 | 1:1 | separate, cross | demo6-like control: λ must not hurt |
| `ctrl_dup` | 1 | – | `p_cross` 0.5 | – | dup, cross | true duplications and fold-backs (`sv_probs` del 0, dup 0.8, inv 0.3, trans 0.2, fback 0.3): λ must not split the species |
| `nest_samehap_6to1` | 2 | nested | one haplotype | 6:1, `num_ecdna_units` 60 | nested, samehap | a light true nested species (about 6–7 copies against 36–43), which a leftover-type prior would drop (added 2026-11-01) |

The bulk copy ratio of two species equals their `ecdna_proportions` ratio whatever the clone
weights: $\psi_{s,k} = \text{scale} \cdot \text{share}_s \cdot w_k$. Here $\psi_{s,k}$ is the
copies of species $s$ per cell of clone $k$, $\text{share}_s$ is the species' normalised
proportion and $w_k$ is the clone's ecDNA weight.

### 1.2 Seeds and the seed-acceptance checks

> **Provenance:** (c) new, amplicon-allele-analysis yubi_dev@d1a1dce.
> **Code:** `simulations/bench_nested.py::check_seed`, `cmd_scenario`.

Each scenario draws its seeds from its own block of 100: 101–, 201–, …, 801–. The walks depend
on the seed alone, so with shared seeds every nested scenario would have the same geometry.

A seed is accepted only if its truth meets the scenario's intent. The checks below are applied
from `truth/`; the thresholds are `MIN_INSIDE` 0.80, `MIN_READS` 20 and `MAX_RATIO_ERR` 0.10:

- **nested:** every nested species has at least 80% of its bases inside the first (outer)
  species;
- **cross:** at least one haplotype-crossing ecDNA junction (phasing AB or BA);
- **samehap:** no crossing, and every species on one haplotype;
- **separate:** the species share no base;
- **dup:** some segment is traversed at least twice by the species;
- **every scenario:** the bulk copy ratio is within 10% of the target, and every ecDNA junction
  has at least 20 reads.

A rejected seed is recorded and the next one tried, at most `--max-tries` (5) per accepted seed.
Every seed is recorded in amplicon-allele-analysis `simulations/configs/bench_nested/seeds.tsv`.

The scenarios have different numbers of accepted seeds:

- nest_samehap_3to1, nest_cross_1to1, ctrl_dup and nest_samehap_6to1 have 6 accepted seeds each
  (grown on 2026-11-01). They were built with `--accept 6 --candidates 30`, so they could draw
  from the first 30 seeds of their block. All 18 new seeds came from the first 10: s305–s308,
  s506–s509, s706–s710 and s801–s806.
- The other four scenarios have 2 accepted seeds each.

### 1.3 Regression cases

The 8 cases of [`results/step_sweep/`](../step_sweep/README.md), with the same subgraphs and
options (`run_case.py::REGRESSION`):

- demo3;
- demo4 `amplicon_0` and `amplicon_1`;
- demo5;
- demo6 union;
- demo7 with `--block-flip unoriented`;
- demo7 total-only;
- demo8.

They count as controls in the λ rule.

## 2. Runs

> **Provenance:** (c) new, RCK yubi_dev@949f796. **Code:** `results/bench_nested/run_case.py::run_case`.

Every case gets two `rck-ecdna` runs. Each starts from the subgraph folder alone, with
`--truth-dir`, `--no-plots`, `--threads 4` and `--time-limit 1800`. Bench cases add
`--block-flip unoriented`.

- `runs/<case>/lambda_0/`: λ = 0, RCK-ecDNA's default.
- `runs/<case>/pool/`: λ = 0.1 with `--keep-pool --pool-lambdas 0 0.02 0.05 0.1 0.2 0.5`.

λ enters only the ranking of decompositions, and every λ > 0 widens the decomposition pool to
the same `--ambiguity-window`. So this one pool holds every answer any λ > 0 can give
([docs/ecDNA.md](../../docs/ecDNA.md), "Offline re-ranking"). `python -m rck.ecdna.rerank
runs/<case>/pool` re-ranks it for any λ with the run's own pruning and ranking functions. It
reports the composites (from `pool/evaluation.tsv`) and the ambiguity fields (from the same-fit
alternatives recorded for each rank-1 solution of the λ grid).

λ = 0 is run rather than re-ranked. A λ = 0 run solves a narrower pool (10% relative gap, 40
solutions), where ties can be broken differently. demo3 in `results/step_sweep/` shows it: at
λ = 0, rank 1 has 12.77 copies; at every λ > 0, an equivalent solution with 12.70 copies. The
λ = 0 runs cost about a quarter of the pool runs (demo8: 8 s against 31 s).

λ grid: $\Lambda = \{0, 0.02, 0.05, 0.1, 0.2, 0.5\}$ (`run_case.py::LAMBDAS`).

**Run infrastructure** — (b) patch, RCK yubi_dev@949f796 and @89a5456 (paths). None of it
changes a result:

- the Slurm array runs at most 2 tasks at once (`%2`), because the academic Gurobi WLS licence
  allows 2 sessions;
- each failed run is attempted up to 5 times in all (`RETRIES` 5), 2 minutes apart;
- each task ends its WLS session with `gurobipy.disposeDefaultEnv()`
  (`run_case.py::release_licence`). An exited process's session otherwise stays counted for a few
  minutes, and 3 tasks of the first array failed on it;
- a run is skipped when its outputs exist;
- the conda environment and run outputs live in `/n/fs/yubi-iw`.

## 3. Rules

`collect.py` re-ranks every case over the λ grid and writes the tables
([Appendix A.3](#a3-outputs)). It applies the four rules below. Their constants are at the top of
`collect.py`:

| constant | value | meaning |
|---|---|---|
| `TOL_REGRESSION` ($\tau_r$) | 0.005 | the largest loss in mean composite a control may take against λ = 0; also the "gain" threshold |
| `TOL_PLATEAU` ($\tau_p$) | 0.001 | λ within this of the best overall mean are on the plateau |
| `MIN_GAIN_CASES` ($g$) | 2 | a λ's shortfall counts only if at least this many cases gain |
| `FAIL_COMPOSITE` | 0.95 | a bench case below this mean composite fails |
| `SCALE_DEPTH_TOL` | 0.3 | relative distance of depth per copy from the bench median that counts as a scale error |
| `SCALE_RATIO_MAX` | 1.5 | segment-to-junction scale ratio beyond which (either way) it is a scale error |
| `JUNCTION_MATCH` | 0.5 | junction recall at which an inferred species "carries" a truth species |

### 3.1 Per-case score

For case $c$ at penalty $\lambda$, let $T_c$ be the in-scope truth species. Let
$\mathrm{comp}_{c,t}(\lambda)$ be the composite of truth species $t$ against the rank-1
decomposition. The composite is defined in
[docs/pipeline/evaluation.md §3.7, composite score](../../docs/pipeline/evaluation.md#37-composite-score); a missed
species scores 0. The case score is

```math
m_c(\lambda) = \frac{1}{|T_c|} \sum_{t \in T_c} \mathrm{comp}_{c,t}(\lambda)
```

At λ = 0 the rank 1 is the `lambda_0` run's (`collect.py::run_row`). For λ > 0 it is the pool's
rank 1 after re-ranking at λ (`collect.py::pool_rows`, `rck.ecdna.rerank.rerank`).

### 3.2 The `merged` flag

> **Provenance:** (c) new, RCK yubi_dev@949f796. **Code:** `collect.py::merged_flag`.

Nested species share reference segments by construction, so segment overlap cannot tell a merged
walk from a correct one. In these simulations the species share no novel adjacency, so the flag
uses junctions. Let $S$ be the significant (not `low_weight`) inferred species of the evaluated
decomposition. Let $R(t, s)$ be the unphased junction recall of truth species $t$ in inferred
species $s$, at 10 kb tolerance (`rck.ecdna.evaluate::junction_scores`). Then

```math
\mathrm{carries}(s) = \{\, t \in T_c : R(t, s) \ge 0.5 \,\}
```

```math
\mathrm{merged} = \Big[\, \exists\, s \in S : |\mathrm{carries}(s)| \ge 2 \,\Big] \;\lor\; \Big[\, \exists\, t \in T_c^{\mathrm{missed}},\ \exists\, s \in S : t \in \mathrm{carries}(s) \,\Big]
```

where $T_c^{\mathrm{missed}}$ is the set of truth species with no one-to-one matched inferred
species in `evaluation.tsv`.

### 3.3 λ selection rule

> **Provenance:** (c) new, RCK yubi_dev@949f796; plateau rule changed @2fdf73d (2026-10-29).
> **Code:** `collect.py::choose`, `collect.py::gain_cases`, `collect.py::is_control`.

Let $C$ be the cases with a score at every $\lambda \in \Lambda$. The controls $K \subseteq C$ are
the regression cases and the `ctrl_*` scenarios. Define

```math
M(\lambda) = \frac{1}{|C|} \sum_{c \in C} m_c(\lambda), \qquad \lambda^{\mathrm{best}} = \operatorname*{arg\,max}_{\lambda \in \Lambda} M(\lambda)
```

```math
G(\lambda) = \{\, c \in C : m_c(\lambda^{\mathrm{best}}) > m_c(\lambda) + \tau_r \,\}
```

```math
P = \{\, \lambda \in \Lambda : M(\lambda) \ge M(\lambda^{\mathrm{best}}) - \tau_p \ \text{ or } \ |G(\lambda)| < g \,\}
```

```math
\mathrm{Reg}(\lambda) = \{\, c \in K : m_c(\lambda) < m_c(0) - \tau_r \,\}
```

```math
\lambda^{\star} = \min \{\, \lambda \in P : \mathrm{Reg}(\lambda) = \varnothing \,\}
```

where:

- $M(\lambda)$ is the mean case score over all complete cases;
- $\lambda^{\mathrm{best}}$ is the grid value with the highest mean (the first such value on a
  tie; pandas `idxmax`);
- $G(\lambda)$ is the set of cases that gain more than $\tau_r$ from $\lambda$ to
  $\lambda^{\mathrm{best}}$;
- $P$ is the plateau;
- $\mathrm{Reg}(\lambda)$ is the set of controls that lose more than $\tau_r$ against λ = 0;
- $\lambda^{\star}$ is the reported recommendation.

`collect.py` also lists the λ values whose rank-1 answers are identical to $\lambda^{\star}$'s in
every case: the same species count and composites string per case. `--lambda X` reports another
λ in `pipeline.tsv` and the scenario table. **RCK's default `--step-penalty` stays 0.** The rule
only makes a recommendation.

**Why the $|G(\lambda)| < g$ clause.** Without it, the 0.001 plateau test alone picked λ = 0.5 on
both the first benchmark and v2. In each case the reason was one case at the edge of the grid,
`nest_samehap_3to1_s301`, which edges up (0.693 at 0.2, 0.740 at 0.5) but stays a failure. The
rule now needs at least two cases to gain before a smaller λ counts as short of the best. With it,
v2 gives 0.2 (was 0.5) and the first benchmark gives 0.05 (was 0.5) (RCK CHANGELOG). The rule is
a heuristic for this benchmark. [citation needed: published practice for choosing a penalty
parameter from a benchmark plateau with a control-regression constraint]

### 3.4 Failure and failure stage

> **Provenance:** (c) new, RCK yubi_dev@949f796; `scale` stage added @2fdf73d.
> **Code:** `collect.py::failure_stage`, `collect.py::scale_check`, `diag/scale_ratio.py::scales`;
> subgraph checks from amplicon-allele-analysis `simulations/bench_nested.py::subgraph_check`.

A bench case **fails** at the chosen λ when any of these hold:

```math
m_c < 0.95 \quad\lor\quad \mathrm{merged} \quad\lor\quad n_c^{\mathrm{sig}} \ne |T_c|
```

where $n_c^{\mathrm{sig}}$ is the number of significant (not `low_weight`) species in rank 1. The
stage that explains a failure is the first that applies:

1. **`subgraph`**: the subgraph misses a true junction, or covers less than 80% (`COVERAGE`) of
   some true segment. This is `subgraph_fail` in `subgraph_check.json`, written by
   amplicon-allele-analysis `bench_nested.py`.
2. **`scale`**: either test below holds.

   ```math
   \left| \frac{d_c}{\tilde d} - 1 \right| > 0.3 \qquad\text{or}\qquad \rho_c = \frac{\sigma_c^{\mathrm{seg}}}{\sigma_c^{\mathrm{junc}}} \notin \left[\tfrac{1}{1.5},\ 1.5\right]
   ```

   with

   ```math
   \sigma_c^{\mathrm{seg}} = \frac{\sum_{i} \ell_i\, \hat x_i}{\sum_{i} \ell_i\, x_i^{\mathrm{true}}}, \qquad
   \sigma_c^{\mathrm{junc}} = \frac{\sum_{a} \hat y_a}{\sum_{a} y_a^{\mathrm{true}}}
   ```

   where:

   - $d_c$ is the subgraph's depth per copy (`scale.json`, through `subgraph_check.json`);
   - $\tilde d$ is its median over the bench cases;
   - $i$ runs over the segments of length $\ell_i \ge 20$ kb that carry ecDNA flow in
     `flow_segments.tsv` and have a truth copy number;
   - $\hat x_i$ is the observed total segment copy number (`obs_total`);
   - $x_i^{\mathrm{true}}$ is the summed copy number of the truth species covering the segment's
     midpoint;
   - $a$ runs over the novel adjacencies that carry ecDNA flow in `flow_adjacencies.tsv` and
     match a truth junction (both ends in the same 2 kb bin, `round(pos / 2000)`);
   - $\hat y_a$ is the observed adjacency copy number (`obs_cn`, reads ÷ junction reads per copy
     θ);
   - $y_a^{\mathrm{true}}$ is the truth copy number summed over the walk steps that traverse the
     junction.

   The ratio bound 1.5 lies between the noise seen on v2 (ratios 1.02–1.40) and the error of a
   wrong ploidy class (about 2, `nest_cross_3sp_s403` at class 6).
3. **`RCK`**: everything else. The note lists, as they apply:
   - the merged species;
   - the wrong species count;
   - ecDNA phase blocks wrong after the flow (`collect.py::flip_check`: blocks over ecDNA segments
     whose truth orientation has purity ≥ 0.9 and that the flow resolved);
   - the imputed edges used;
   - the composites.

**Why the `scale` stage.** On v2, `nest_cross_3sp_s403` failed because `rck_inputs.py` fitted
ploidy ≈ 6 and halved depth per copy (6.51 against about 13.8). The old rule called it an `RCK`
failure because it checked only coverage and missing junctions. The stage separates input-scale
errors from decomposition errors (RCK CHANGELOG, "a `scale` failure stage").

### 3.5 Method-agnostic metrics (opt-in)

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. **Code:** `collect.py::aa_metrics`,
> `diag/aa_regression.py`.

`collect.py --aa-metrics` scores every (case, λ) a second time through the AA-format path:

1. `rck.ecdna.aa_export` writes the rank-1 answer and its first 2 same-fit alternatives;
2. `rck.ecdna.aa_eval` scores them against the truth mixture graph of amplicon-allele-analysis
   `simulations/truth_mixture.py`.

It appends the main sample-level columns (`strict.*`, `best_of.*`, `aa_path_matches`) to
`bench_sweep.tsv`, after the existing ones. It also writes `aa_metrics.tsv` and the block between
the `aa-metrics` markers below. The λ rule and the failure rule read only the existing columns, so
the flag does not change them.

The metrics, the strict and best-of modes, and the CoRAL code and supplementary-methods variants
are defined in [docs/pipeline/evaluation.md §4, method-agnostic
scorer](../../docs/pipeline/evaluation.md#4-the-method-agnostic-scorer-aa_eval). (The text that `collect.py` writes into the block still points to
"docs/ecDNA.md, section 6", the location before this rewrite.)

## 4. Results (grown bench, 40 cases, 2026-11-01)

The block below is written by `python collect.py` over the 40 cases:

- the 22 re-baseline runs (RCK Slurm 32041913);
- the 18 new cases of step d (RCK Slurm 32050408).

All were run on subgraphs built with the `clip_path` fix and `--rck-ploidy 2`. Do not edit
between the markers: `collect.py` replaces that text.

<!-- results:start -->
**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.2**, the smallest λ on the plateau of best mean composite over all 40 cases (0.2, 0.5; a λ whose shortfall to the best comes from fewer than 2 cases is on it), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.2. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.7713 | 0.8155 | 0.8417 | 0.8417 | 0.8531 | 0.8542 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |
| cases gaining from λ to the best λ | demo7_totalonly, nest_cross_1to1_s507, nest_cross_2to1_s104, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301, nest_samehap_3to1_s305, nest_samehap_3to1_s308 | demo7_totalonly, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_samehap_3to1_s301 | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.2 | merged at 0 | merged at 0.2 |
|---|---|---|---|---|---|
| ctrl_dup | 6 | 0.856 | 0.856 | 0 | 0 |
| ctrl_separate | 2 | 0.979 | 0.979 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 6 | 0.813 | 0.889 | 3 | 3 |
| nest_cross_2to1 | 2 | 0.609 | 0.967 | 1 | 0 |
| nest_cross_3sp | 2 | 0.533 | 0.761 | 2 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 6 | 0.744 | 0.843 | 2 | 1 |
| nest_samehap_6to1 | 6 | 0.553 | 0.553 | 0 | 0 |

Remaining failures at λ = 0.2 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 17 of 32 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.964,ecdna2=0.423 |
| nest_samehap_3to1_s307 | RCK | 1 species for 2 true; composites ecdna1=0.985,ecdna2=0.000 |
| nest_samehap_3to1_s308 | RCK | composites ecdna1=0.973,ecdna2=0.924 |
| nest_cross_3sp_s403 | RCK | merged: ec1 carries ecdna1+ecdna2; ec2 carries ecdna1+ecdna3; composites ecdna1=0.519,ecdna2=0.593,ecdna3=0.523 |
| nest_cross_1to1_s502 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.967,ecdna2=0.943 |
| nest_cross_1to1_s505 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.944,ecdna2=0.586 |
| nest_cross_1to1_s507 | RCK | merged: ec3 carries ecdna1+ecdna2; 4 species for 2 true; composites ecdna1=0.969,ecdna2=0.952 |
| nest_cross_1to1_s508 | RCK | 3 species for 2 true; composites ecdna1=0.683,ecdna2=0.662 |
| ctrl_dup_s703 | RCK | 2 species for 1 true; composites ecdna1=0.712 |
| ctrl_dup_s708 | RCK | 2 species for 1 true; composites ecdna1=0.732 |
| ctrl_dup_s710 | RCK | 3 species for 1 true; imputed edges imputed_adj_intra_chr2_10; composites ecdna1=0.723 |
| nest_samehap_6to1_s801 | RCK | 1 species for 2 true; composites ecdna1=0.984,ecdna2=0.000 |
| nest_samehap_6to1_s802 | RCK | 1 species for 2 true; composites ecdna1=0.996,ecdna2=0.000 |
| nest_samehap_6to1_s803 | RCK | 1 species for 2 true; composites ecdna1=0.986,ecdna2=0.000 |
| nest_samehap_6to1_s804 | RCK | imputed edges imputed_adj_intra_chr3_10; composites ecdna1=0.987,ecdna2=0.695 |
| nest_samehap_6to1_s805 | RCK | 1 species for 2 true; composites ecdna1=0.991,ecdna2=0.000 |
| nest_samehap_6to1_s806 | RCK | 1 species for 2 true; composites ecdna1=0.995,ecdna2=0.000 |

Scale: 1 of 32 cases have depth per copy more than 30% from the median or segment / junction scales apart by more than 1.5x (nest_cross_2to1_s104).
Subgraph stage: 0 of 32 subgraphs miss a true junction or cover under 80% of a true segment.

Compute: 7030 s of RCK wall time over 40 cases (one λ = 0 run and one pool run each, 4 threads).
<!-- results:end -->

**Reading.** λ = 0.2 is the smallest grid value with every real gain. Each new failure of the
grown bench is at the `RCK` stage. All 18 new subgraphs pass the subgraph stage, and only s104 is
flagged at the `scale` stage, as before. Per grown scenario at λ = 0.2 (mean over 6 seeds):

| scenario | mean | failing seeds |
|---|---|---|
| nest_samehap_3to1 | 0.843 | 3: s301, s307, s308 |
| nest_cross_1to1 | 0.889 | 4: s502, s505, s507, s508 |
| ctrl_dup | 0.856 | 3: s703, s708, s710 |
| nest_samehap_6to1 | 0.553 | 6 (in 5 of them the nested species scores 0) |

### 4.1 Method-agnostic metrics (2026-11-03)

The block below is `collect.py --aa-metrics` on the 40 cases; no run was repeated. Definitions:
[docs/pipeline/evaluation.md §4, method-agnostic scorer](../../docs/pipeline/evaluation.md#4-the-method-agnostic-scorer-aa_eval). Do not
edit between the markers.

<!-- aa-metrics:start -->
**Method-agnostic metrics** (`python collect.py --aa-metrics`; `aa_metrics.tsv`; definitions in docs/ecDNA.md, section 6, "Method-agnostic scorer"). Every case and λ was scored again from AA-format files only (`rck.ecdna.aa_export` + `rck.ecdna.aa_eval`): 240 of 240 (case, λ) reproduce the existing composites, species counts and merged flags exactly.

| metric | bench, λ = 0 | all, λ = 0 | bench, λ = 0.2 | all, λ = 0.2 |
|---|---|---|---|---|
| composite, strict (rank 1, as today) | 0.735 | 0.771 | 0.820 | 0.853 |
| composite, best-of (rank 1 + 2 same-fit alternatives, every cycle) | 0.845 | 0.862 | 0.917 | 0.931 |
| composite-H0, strict | 0.699 | 0.738 | 0.799 | 0.835 |
| composite-H0, best-of | 0.806 | 0.827 | 0.897 | 0.914 |
| exact cycles (100 bp), strict | 32 / 60 | 42 / 72 | 43 / 60 | 55 / 72 |
| exact cycles (100 bp), best-of | 38 / 60 | 48 / 72 | 49 / 60 | 61 / 72 |
| CoRAL cycle-interval overlap (CoRAL code), strict | 0.737 | 0.772 | 0.839 | 0.871 |
| CoRAL cyclic LCS (CoRAL code), strict | 0.730 | 0.766 | 0.839 | 0.871 |
| CoRAL cyclic LCS (CoRAL code), best-of | 0.845 | 0.861 | 0.948 | 0.958 |
| secondary: cycle-interval overlap (CoRAL supplementary methods), strict | 0.821 | 0.846 | 0.872 | 0.898 |
| secondary: cyclic LCS (CoRAL supplementary methods), strict | 0.836 | 0.858 | 0.877 | 0.901 |
| secondary: cyclic LCS (CoRAL supplementary methods), best-of | 0.983 | 0.986 | 0.983 | 0.986 |
| species-count error (mean absolute, inferred − true) | 0.531 | 0.425 | 0.500 | 0.400 |
| missed species | 8 | 8 | 6 | 6 |
| false species | 9 | 9 | 10 | 10 |
| AA error-free reconstruction (share of cases, every amplicon, 100 bp) | 0.000 | 0.000 | 0.000 | 0.000 |
| AA error-free, benign copy-number splits merged (share of cases) | 0.906 | 0.897 | 0.906 | 0.897 |
| CoRAL breakpoint-graph accuracy | 1.000 | 1.000 | 1.000 | 1.000 |
| novel-adjacency FNR | 0.000 | 0.000 | 0.000 | 0.000 |
| novel-adjacency FDR | 0.008 | 0.007 | 0.008 | 0.007 |
| allele segment CN error (copies per base) | 1.452 | 1.318 | 1.452 | 1.318 |
| infinite-sites violations (sum) | 19 | 25 | 19 | 25 |
| CN explained, heaviest cycle | 0.709 | 0.712 | 0.705 | 0.703 |
| CN explained, 3 heaviest cycles | 0.892 | 0.881 | 0.885 | 0.876 |
| long-read subwalks satisfied, strict | 0.802 | 0.802 | 0.833 | 0.833 |
<!-- aa-metrics:end -->

The graph truth is the mixture graph of amplicon-allele-analysis `simulations/truth_mixture.py`,
written to `truth/aa/` of the 32 bench simulations (in `/n/fs/yubi-iw/runs/`). The tracked demo
runs have theirs in amplicon-allele-analysis (`simulations/runs/demo{5,6,7,8}_*/truth/aa/mixture_*`,
since 2026-11-03). demo4's (RCK's `tests/ecdna/data/demo4_truth`) is still built under
`/n/fs/yubi-iw/runs/bench_nested_aa/mixtures/`. Moving demo5–demo8 to the repository copies changed
only the paths in `aa_eval.json`. The copies are byte-identical to the ones beside the
simulations, and every graph and species table of those cases is unchanged.

- **Regression.** `diag/aa_regression.py` exports each answer to AA format, scores it with
  `aa_eval`, and compares it with today's scores. It covers λ = 0 (the λ = 0 run's rank 1) and
  λ = 0.2 (the pool's rank 1), 80 (case, λ) in all, compared at full precision. Identical:
  - every part of the composite: segment and haplotype Jaccard, junction F1, phased junction F1,
    phasing accuracy, order match, copy-number error;
  - every composite, the matched species, `n_species`, `n_truth` and `merged`.

  `collect.py --aa-metrics` repeats the check for all 240 (case, λ). It leaves unchanged the 19
  existing columns of `bench_sweep.tsv`, `pipeline.tsv`, `rerank_validation.tsv` and the results
  section above. Two things were needed:
  - A species' other orderings: `nest_samehap_3to1_s304` matches its order only through ordering
    2, which RCK's own `aa/` files leave out.
  - Reading the old tables with `float_precision="round_trip"`: pandas' default parser reads some
    stored values one ulp off.
- **Graph level (100 bp).** Only 1 of the 42 truth amplicons is an AA error-free reconstruction
  (demo6's chr17 amplicon). 38 of 42 are once benign splits are merged (`error_free_merged`,
  added 2026-11-03): neighbouring segments joined only by a concordant edge, with no discordant
  edge at the boundary, are merged in both graphs. Of the 4 left:
  - two use an imputed edge (below);
  - two (`demo4_amp0`, `nest_samehap_6to1_s803`) miss a copy-number-2 flank segment inside the
    truth interval that the subgraph shortens.

  In 40 of the 42 every true discordant edge is found and none is false. The other two,
  `ctrl_dup_s710` and `nest_samehap_6to1_s804`, use one imputed edge each (novel-adjacency FDR
  0.007 overall). Every failure therefore comes from segments. All 532 inferred segment boundaries
  that match no true boundary are imputed extremities of the subgraph (`subgraph_construct.py`):
  its 10 kb grid, and pieces about 3 kb from a junction.

  Breakpoint-graph accuracy is 1.000, because every subgraph holds every true junction (subgraph
  stage above). The allele-specific segment copy-number error is 1.3 copies per base: 1–11% of
  the truth's copy number on the bench cases (worst `ctrl_dup_s710`, 6.4 copies), and up to 17% on
  demo4. RCK's graphs have 25 infinite-sites violations, the same 25 as the truth, because the
  simulated duplications reuse extremities.
- **Species level.** Best-of gives 0.931 at λ = 0.2, against 0.853 strict. Best-of takes rank 1,
  its first 2 same-fit alternatives, and every cycle of them including `low_weight` ones. Its
  score is above the decomposition-level set metric of NEXT_STEPS (k = 2: 0.923), as it should be,
  because best-of picks per truth species. composite-H0 is 0.02 below the composite, so nothing in
  the composite rests on phasing alone. Strict, 55 of 72 truth species are exact cycles at 100 bp
  (61 best-of).

  The CoRAL cycle metrics are reported as CoRAL's code computes them (headline since 2026-11-03):
  cyclic LCS 0.871 strict and 0.958 best-of at λ = 0.2. The supplementary-methods version
  (secondary, `*_supp`) is lenient by construction (0.986 best-of). It normalises by the true
  cycle's length, so a merged walk that contains a nested species scores 1 for it. The code's
  normaliser (the bins' total with multiplicity) charges the merged walk's extra bases.
- **Sample level.** At λ = 0.2 the rank-1 answers miss 6 truth species and add 10 false ones (the
  6to1 and ctrl_dup failures of the table above). The heaviest cycle explains 0.70 of the graph's
  length-weighted copy number, and the three heaviest 0.88. 83% of the multi-junction long-read
  paths appear in order in some rank-1 cycle.

## 5. Known limitations and open failures

The full analysis is in [NEXT_STEPS.md](NEXT_STEPS.md), "Grown bench: failure analysis and
revised plan (2026-11-01)", with `diag/grow_failures.py` and `diag/grow_rules.py`. It uses the
saved pools only. In short:

- **nest_samehap_6to1.** The light nested species is split into 2–3 sub-threshold `low_weight`
  walks, so rank 1 has one significant species. The data cannot tell this from the right answer.
- **ctrl_dup s708 and s710.** The truth traverses a segment 4 times, above `rck-ecdna
  --max-multiplicity` 3, so it cannot be expressed as one species. amplicon-allele-analysis
  `check_seed` does not compare the two.
- **`low_weight` disagreement.** The output re-flags `low_weight` on refitted weights, while the
  ranking uses integer weights. They disagree in 8 of the 40 rank-1 answers, in both directions.
- **Decomposition of a noisy integer flow** (s301, s703, s505, s403; Appendix B.5). The
  decomposition must cover an integer flow exactly. A 1–3 copy dip from simulated depth bias costs
  as much as one species, so exact-fit alternatives win (extra leftovers or wrong merges). The
  opt-in noise-tolerant residual did not change any rank 1 (`diag/rck_noise/SUMMARY.md`).
- **Single seed.** Every number is one Gurobi run per case. The graph route's round trip found
  the decomposition changes in 29 of 36 cases at an unchanged flow objective
  ([../graph_route/README.md](../graph_route/README.md)), so case-level numbers carry
  tie-breaking noise.
- **The rule depends on one case.** The choice of 0.2 over 0.1 rests on `nest_cross_3sp_s402`,
  whose switch lies between λ = 0.10 and 0.11 (Appendix B.5).
- **The composite and the scale stage.** The composite scores s403's 2× scale error (0.592) above
  the right scale (0.545); read the `scale` stage, not that difference.

## Appendix A. How to run (Princeton ionic)

### A.1 Commands

```bash
# 1. simulations, seed checks, subgraphs and their checks (no Gurobi; 8 tasks, one per scenario).
#    The tracked subgraphs are built at ploidy class 2 (re-baseline, 2026-10-31); without the
#    variable, rck_inputs.py fits the class, as before.
#    The sbatch's default (--accept 2) rebuilds 2 seeds per scenario and rewrites the grown
#    scenarios' seed records to 2 slots, so build those four with the step d arguments:
cd ~/amplicon-allele-analysis/simulations && mkdir -p logs
BENCH_NESTED_RCK_PLOIDY=2 sbatch --array=0,1,3,5 bench_nested.sbatch
BENCH_NESTED_RCK_PLOIDY=2 BENCH_NESTED_ARGS="--accept 6 --candidates 30" sbatch --array=2,4,6,7 bench_nested.sbatch

# 2. first check the re-ranking on the step sweep's hardest cases (demo3, demo7_totalonly, demo8)
cd ~/RCK/results/bench_nested && mkdir -p logs && sbatch --array=0,6,7%2 bench_nested.sbatch
python collect.py --validate          # rerank_validation.tsv: every "pool" row should match

# 3. every case (finished runs are skipped). The sbatch's own --array is still 0-21 (the 22-case
#    bench), so pass the 40 cases explicitly; `python run_case.py list` gives the index of each.
sbatch --array=0-39%2 bench_nested.sbatch   # 2 at a time (academic WLS licence)
python collect.py                     # bench_sweep.tsv, pipeline.tsv, results above
python collect.py --aa-metrics        # optional: aa_metrics.tsv and the method-agnostic block
```

Both arrays skip every step whose outputs exist, so a failed or interrupted array can be
resubmitted as is. Estimated cost:

- simulation: about 1.5 min per seed;
- subgraph: about 1 min per seed;
- in all: ≈ 14–25 simulations with rejections, ≈ 0.5–1 CPU-hour;
- RCK: about 10–60 s per λ = 0 run and 30 s–3 min per pool run with 4 threads. For 22 cases that
  is ≈ 1–3 CPU-hours, about 30 min of wall time at 2 concurrent tasks.

The shell may carry stale `BENCH_NESTED_ROOT` / `BENCH_NESTED_OUT` from earlier work; unset them
before `collect.py` unless you mean to use them.

### A.2 Environment variables and locations

- `BENCH_NESTED_ROOT`, as for amplicon-allele-analysis `bench_nested.py`: take the bench cases from
  a benchmark built elsewhere (`$BENCH_NESTED_ROOT/configs/seeds.tsv`, `subgraphs/<run>/union` and
  `runs/<run>`).
- `BENCH_NESTED_OUT`: put the RCK runs elsewhere.
- `BENCH_NESTED_RESULTS`: write `collect.py`'s tables and a results-only `README.md` to another
  folder, so a new benchmark does not overwrite this one's.
- `BENCH_NESTED_AA_OUT`: where `--aa-metrics` writes its exported files (default
  `/n/fs/yubi-iw/runs/bench_nested_aa`).
- `AAA`: the amplicon-allele-analysis checkout (default `~/amplicon-allele-analysis`).

Unset, none of them changes anything.

Run outputs live in project space, `/n/fs/yubi-iw/runs/` (moved from scratch on 2026-11-04).

- The local `runs` link points at `/n/fs/yubi-iw/runs/bench_nested_v2/grow/rck_runs` (since
  2026-11-01). That folder holds the 18 new cases' runs, plus links to the re-baseline's runs in
  `bench_nested_v2/rebuild/rck_runs` for the other 22. The 8 regression cases there link on to
  `/n/fs/yubi-iw/runs/bench_nested/rck_runs`.
- amplicon-allele-analysis `simulations/runs/bench_nested` points at
  `/n/fs/yubi-iw/runs/bench_nested_v2/grow/runs`: the new simulations plus links to v2's.
- The diagnostics are in `/n/fs/yubi-iw/runs/bench_nested_v2/diag/`.

Everything this README or `NEXT_STEPS.md` needs to keep is in the repositories. The simulations
can be regenerated from the configs with `bench_nested.py`, and the RCK runs with
`bench_nested.sbatch`.

### A.3 Outputs

- **`rerank_validation.tsv`**: rank-1 composites and ambiguity fields, re-ranked (and the λ = 0
  run), against `results/step_sweep/sweep.tsv`, per case and λ.
- **`bench_sweep.tsv`**: one row per case and λ. Columns:
  - rank-1 significant species against the truth count, weights, steps, residual;
  - composite per truth species and their mean;
  - `merged` (§3.2);
  - `ambiguous`, `decided_by`, `gap`, `n_alternatives`;
  - runtime. `source` is `run` for λ = 0 and `pool` for re-ranked λ, whose runtime is the pool
    run's.

  With `--aa-metrics`, the main method-agnostic columns follow.
- **`pipeline.tsv`**: one row per bench case at the chosen λ. Columns:
  - subgraph scores: segments and junctions recovered; ecDNA phase blocks oriented right, wrong
    or not at all; and the junction reads per copy `rck_inputs.py` would have fitted against the
    6.5 used;
  - the scale columns (`depth_per_copy`, `depth_vs_median`, `segment_scale`, `junction_scale`,
    `scale_ratio`, `scale_fail`);
  - the flow objective and the imputed edges used;
  - block flips against the truth orientation;
  - the RCK scores;
  - `failure_stage` and `failure` (§3.4).
- **`aa_metrics.tsv`** (`--aa-metrics`): one row per case and λ, with every sample-level metric
  of `rck.ecdna.aa_eval` for the strict and best-of modes (`strict.*`, `best_of.*`), and
  `aa_path_matches`: whether the AA-format path reproduced the case's composites, species count
  and merged flag. Per-species and per-amplicon tables, and the exported AA files, go to
  `/n/fs/yubi-iw/runs/bench_nested_aa/` (`BENCH_NESTED_AA_OUT`). `diag/aa_regression.py` writes
  the same there, with the field-by-field comparison as `aa_regression.tsv`.
- **`comparison.md`, `comparison.tsv`**: the first benchmark against v2.
- **`runs/`**: every run's full output (not tracked).

### A.4 Diagnostic scripts (`diag/`)

| script | what it does |
|---|---|
| `aa_regression.py` | AA-format re-scoring and field-by-field comparison (§3.5) |
| `scale_ratio.py` | segment and junction scales against the truth; `scales(run)` is used by `collect.py` |
| `rck_diag.sbatch`, `summarize.sh` | diagnostic RCK runs; `DIAG_POOL_LAMBDA` sets the pool run's λ, `DIAG_LAMBDA0=0` skips the λ = 0 run |
| `pool_summary.py`, `flow_vs_truth.py`, `obs_vs_truth.py` | pool contents; flow against truth; segment-depth and junction-read errors against the truth copy numbers |
| `noise_summary.py` | per case and setting of the noise-tolerant residual: rank 1, the pool's best, inexact fits; `--archive` copies small outputs into `diag/rck_noise/runs/` |
| `leftover_recount.py` | re-flags saved pool species under candidate leftover rules and re-ranks with `rank_pool` |
| `grow_failures.py`, `grow_rules.py` | per-case anatomy of rank 1 against the pool's best; re-ranking of all pools under candidate rules (outputs in `grow_2026-11-01/`) |

## Appendix B. Earlier runs

The sections below are kept as they were written. Their numbers are those of the runs named. The
λ rule and the failure stages changed on 2026-10-29 (§3.3, §3.4); each section says which rule it
used.

### B.1 History

- **First benchmark (2026-10-25).** Per-step haplotype labels, 22 cases (14 bench).
- **bench_nested_v2 (2026-10-27).** Rebuilt with inherited haplotype labels, then adopted on
  2026-10-29.
- **Re-baseline (2026-10-31; step c of `NEXT_STEPS.md`).** The same 14 simulations, with subgraphs
  rebuilt with the long-read path fix and an explicit `--rck-ploidy 2`
  (`BENCH_NESTED_RCK_PLOIDY=2`), and RCK rerun on all 22 cases.
- **Grown bench (2026-11-01; step d).** 4 more accepted seeds each for nest_samehap_3to1,
  nest_cross_1to1 and ctrl_dup, and the new nest_samehap_6to1 scenario with 6 seeds. That adds 18
  simulations, for 32 bench cases and 40 cases in all; the 22 existing cases were not rerun.

Tables in git history:

- the first benchmark's: 9b74b87;
- v2 before the re-baseline: 7778991;
- the 22-case re-baseline: 143c24d.

### B.2 Grown bench (2026-11-01)

Step d of `NEXT_STEPS.md`. Runs:

- simulations and subgraphs: amplicon-allele-analysis Slurm 32049944
  (`BENCH_NESTED_ARGS="--accept 6 --candidates 30"`, `BENCH_NESTED_RCK_PLOIDY=2`, root
  `/n/fs/yubi-iw/runs/bench_nested_v2/grow`);
- RCK on the 18 new cases only: Slurm 32050408 (18 of 18 COMPLETED; 7030 s of RCK wall time over
  all 40 cases).

The run record is in `grow_2026-11-01/` (`README.md`, `collect.log`, the three tables,
`compare_existing.txt`).

- **The 22 existing cases are unchanged.** `bench_sweep.tsv` (132 rows) and `rerank_validation.tsv`
  (48 rows) are identical row for row. In `pipeline.tsv` the 14 existing rows differ in one column
  only, `depth_vs_median`, by 0.000–0.005. That column is relative to the median depth per copy
  over the bench cases, which moves from 13.762 (14 cases) to 13.7065 (32 cases). Each case's
  `depth_per_copy` and `scale_fail` are unchanged. (`compare_existing.txt` compared
  `bench_sweep.tsv` only.)
- **Adopted** into amplicon-allele-analysis (identical to `grow/configs` and `grow/subgraphs`):
  - the configs: the 21 new ones, including the 3 rejected ctrl_dup seeds, as for the existing
    scenarios;
  - `seeds.tsv`;
  - the 18 new subgraph folders and `subgraph_checks.tsv`.

  And here, the three tables. `python collect.py` reproduces all three tables and the results
  block byte for byte when run with no `BENCH_NESTED_*` variables, the tracked subgraphs and the
  repointed `runs` links.
- Mean composite at λ = 0.2: 0.9262 over 22 cases → **0.8531 over 40**. λ = 0.2 is still the
  recommendation. Failures: 17 of 32 bench cases (5 of 14 before).

### B.3 Re-baseline (2026-10-31)

Step c of `NEXT_STEPS.md`. The 14 subgraphs were rebuilt from the unchanged v2 simulations, from
amplicon-allele-analysis `simulations/`, with
`BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2/rebuild BENCH_NESTED_RCK_PLOIDY=2 python bench_nested.py rebuild`.
RCK was run on all 22 cases (Slurm 32041913, outputs in `rebuild/rck_runs`). The run record and
the comparison with v2 are in `rebuild_2026-10-31/` (`comparison_vs_v2.txt`).

- Graphs:
  - `nodes`, `edges`, `segment_blocks`, `phase_blocks` and `cn_profile` are byte-identical in 13
    of 14 cases;
  - `paths.tsv` grows in all 14 (0–10 → 6–23 paths);
  - s403 changes copy numbers only (ploidy class 6 → 2).
- Rank 1 at λ = 0.2 is unchanged in 21 of 22 cases. The path fix moves no composite.
- `nest_cross_3sp_s403`: 0.592 → 0.545.
  - Its scale is now right: depth per copy 13.92, 1% from the median; segment and junction CN
    1.03× and 1.04× truth. So `pipeline.tsv` puts its failure at the `RCK` stage, not `scale`.
  - Rank 1 is exactly the 2026-10-27 `--rck-ploidy 2` diagnostic's: the same flow and
    observations; species at 5.75 / 4.06 / 3.33; ec1 carries ecdna1 + ecdna2, ec2 carries
    ecdna1 + ecdna3.
  - The pool's best is again 0.729 (#15 of 83; #16 of 85 in the step 2 w1 pool), not picked. The
    failure is in the decomposition, not the scale.
  - The composite scores the 2× scale error (0.592) above the right scale (0.545); read the
    `scale` stage, not that drop.
- Mean composite at λ = 0.2: 0.9283 → 0.9262 (−0.0021, all from s403). λ = 0.2 is still the
  recommendation. `collect.py` with `BENCH_NESTED_OUT=.../rebuild/rck_runs` and the tracked
  subgraphs reproduces the three tables byte for byte.
- `subgraph_checks.tsv` gains `seg_calls_borderline`, and all 14 `verify_subgraph.py` verdicts
  are PASS (6 were CHECK: s103, s104, s402, s403, s502, s505). Segment calls are now scored
  against summed haplotype copies (amplicon-allele-analysis CHANGELOG, `verify_subgraph.py`).

#### Re-baseline results (22 cases, 2026-10-31; current rule)

From RCK Slurm 32041913 on the 14 rebuilt subgraphs and the 8 regression cases; the tables are in
git history (143c24d).

**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.2**, the smallest λ on the plateau of best mean composite over all 22 cases (0.2, 0.5; a λ whose shortfall to the best comes from fewer than 2 cases is on it), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.2. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.8253 | 0.8579 | 0.9055 | 0.9055 | 0.9262 | 0.9283 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |
| cases gaining from λ to the best λ | demo7_totalonly, nest_cross_2to1_s104, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | demo7_totalonly, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_samehap_3to1_s301 | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.2 | merged at 0 | merged at 0.2 |
|---|---|---|---|---|---|
| ctrl_dup | 2 | 0.855 | 0.855 | 0 | 0 |
| ctrl_separate | 2 | 0.979 | 0.979 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 2 | 0.860 | 0.860 | 2 | 2 |
| nest_cross_2to1 | 2 | 0.609 | 0.967 | 1 | 0 |
| nest_cross_3sp | 2 | 0.533 | 0.761 | 2 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 2 | 0.836 | 0.836 | 1 | 1 |

Remaining failures at λ = 0.2 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 5 of 14 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.964,ecdna2=0.423 |
| nest_cross_3sp_s403 | RCK | merged: ec1 carries ecdna1+ecdna2; ec2 carries ecdna1+ecdna3; composites ecdna1=0.519,ecdna2=0.593,ecdna3=0.523 |
| nest_cross_1to1_s502 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.967,ecdna2=0.943 |
| nest_cross_1to1_s505 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.944,ecdna2=0.586 |
| ctrl_dup_s703 | RCK | 2 species for 1 true; composites ecdna1=0.712 |

Scale: 1 of 14 cases have depth per copy more than 30% from the median or segment / junction scales apart by more than 1.5x (nest_cross_2to1_s104).
Subgraph stage: 0 of 14 subgraphs miss a true junction or cover under 80% of a true segment.

Compute: 2762 s of RCK wall time over 22 cases (one λ = 0 run and one pool run each, 4 threads).

### B.4 v2 results before the re-baseline (2026-10-29; current rule)

From Slurm 32000320 on the v2 subgraphs before the rebuild: no long-read path fix, and the ploidy
class fitted (6 for `nest_cross_3sp_s403`). The tables are in git history (7778991). Mean
composite at λ = 0.2: 0.9283, against 0.9262 after the re-baseline; the whole difference is s403.

**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.2**, the smallest λ on the plateau of best mean composite over all 22 cases (0.2, 0.5; a λ whose shortfall to the best comes from fewer than 2 cases is on it), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.2. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.8269 | 0.8600 | 0.9076 | 0.9076 | 0.9283 | 0.9304 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |
| cases gaining from λ to the best λ | demo7_totalonly, nest_cross_2to1_s104, nest_cross_3sp_s402, nest_cross_3sp_s403, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | demo7_totalonly, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_samehap_3to1_s301 | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.2 | merged at 0 | merged at 0.2 |
|---|---|---|---|---|---|
| ctrl_dup | 2 | 0.855 | 0.855 | 0 | 0 |
| ctrl_separate | 2 | 0.979 | 0.979 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 2 | 0.860 | 0.860 | 2 | 2 |
| nest_cross_2to1 | 2 | 0.609 | 0.967 | 1 | 0 |
| nest_cross_3sp | 2 | 0.550 | 0.785 | 2 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 2 | 0.836 | 0.836 | 1 | 1 |

Remaining failures at λ = 0.2 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 5 of 14 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.964,ecdna2=0.423 |
| nest_cross_3sp_s403 | scale | depth per copy 6.51 is -53% from the median 13.74; segment / junction scale 2.12 (segment CN 2.20x, junction CN 1.04x truth); composites ecdna1=0.635,ecdna2=0.483,ecdna3=0.657 |
| nest_cross_1to1_s502 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.967,ecdna2=0.943 |
| nest_cross_1to1_s505 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.944,ecdna2=0.586 |
| ctrl_dup_s703 | RCK | 2 species for 1 true; composites ecdna1=0.712 |

Scale: 2 of 14 cases have depth per copy more than 30% from the median or segment / junction scales apart by more than 1.5x (nest_cross_2to1_s104, nest_cross_3sp_s403).
Subgraph stage: 0 of 14 subgraphs miss a true junction or cover under 80% of a true segment.

Compute: 4349 s of RCK wall time over 22 cases (one λ = 0 run and one pool run each, 4 threads).

### B.5 bench_nested_v2 (2026-10-27; old λ rule and failure stages)

These are the results of the rebuilt benchmark (inherited haplotype labels), run on ionic (Slurm
32000320, the 14 bench cases; the 8 regression cases are links to the first benchmark's runs).
They went to `/n/fs/yubi-iw/runs/bench_nested_v2/rck_results/` and became this folder's tables,
with `comparison.md` and `comparison.tsv` against the first benchmark's. This section is the
2026-10-27 reading, under the old λ rule and failure stages.

**Why v2 was built.** The first benchmark was simulated with per-step haplotype labels. These give
some junctions several phasings and some fold-backs a haplotype crossing. Neither can arise from
one physical join, and they caused part of the `ctrl_dup_s702` and `nest_cross_3sp_s402`
failures. See amplicon-allele-analysis `CHANGELOG.md`, `hap_labels`.

v2 was first built and run outside the tracked folders. The regression cases are unchanged, so
their existing runs were linked in rather than rerun:

```bash
cd ~/RCK/results/bench_nested && mkdir -p logs
export BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2      # cases: its seeds.tsv and subgraphs
export BENCH_NESTED_OUT=$BENCH_NESTED_ROOT/rck_runs                # runs
export BENCH_NESTED_RESULTS=$BENCH_NESTED_ROOT/rck_results         # collect.py tables and README
mkdir -p $BENCH_NESTED_OUT
for c in demo3 demo4_amp0 demo4_amp1 demo5 demo6_union demo7_blockflip demo7_totalonly demo8; do
  ln -sfn /n/fs/yubi-iw/runs/bench_nested/rck_runs/$c $BENCH_NESTED_OUT/$c
done
sbatch --array=8-21%2 bench_nested.sbatch                          # the 14 bench cases (sbatch passes the environment)
python collect.py                                                  # afterwards, with the same three variables set
```

The two benchmarks share scenarios and seed blocks. 9 bench cases keep their name:

- the 4 samehap cases are identical;
- the 5 cross and control cases have the same walks and copy numbers but new haplotype labels.

5 are replacement seeds: s104, s505, s603, s703 and s403 (in place of s406). They replace s101,
s501, s602 and s701, which have no crossing under inherited labels; s403 now passes its read
check.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all 22 cases | 0.8269 | 0.8600 | 0.9076 | 0.9076 | 0.9283 | 0.9304 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |

- `ctrl_dup_s702` is fixed: 0.744 → 0.997 at every λ. It has no wrong junction label, no imputed
  edge and no ecDNA block wrong after the flow.
- `ctrl_separate_s601` is fixed by `--max-extend`: 0.759 → 0.982. Its replacement pair, s603,
  scores 0.976.
- `nest_cross_3sp_s402` has clean labels, flips and imputed edges. It is merged at λ ≤ 0.1
  (ecdna2 + ecdna3, 2 species for 3, 0.523), and recovered at λ ≥ 0.2 (3 species, 0.978).
- **The λ = 0.05–0.5 plateau no longer holds**, and 0.1 is no longer on it. λ = 0.2 is the
  smallest grid value with every real gain.
  - The old rule picks 0.5 again, for `nest_samehap_3to1_s301` alone (0.693 → 0.740, still
    merged).
  - The choice of 0.2 rests on one case. Re-ranking the s402 pool on a 0.01 grid puts its switch
    between 0.10 and 0.11. At 0.1, one extra species (penalty 1.0) exactly ties 10 fewer steps, so
    0.1 sits on the tie.
- `nest_cross_1to1_s502` finishes (13 s + 32 s). `nest_cross_3sp_s406` is replaced by s403.

Remaining failures at λ = 0.1, 6 of 14 bench cases (5 at λ = 0.2):

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | unchanged: ec2 carries ecdna1+ecdna2; 3 species for 2; 0.693 |
| nest_cross_3sp_s402 | RCK (λ) | ec2 carries ecdna2+ecdna3; 2 species for 3; 0.523; fixed at λ ≥ 0.2 |
| nest_cross_3sp_s403 | subgraph (scale), then RCK | `rck_inputs` ploidy fit 5.81 (depth per copy 6.51, about 13.8 elsewhere): segment CN 2× truth, flow objective 41, 181 steps; 0.592; `pipeline.tsv` says RCK. Rebuilt with `--rck-ploidy 2`, the scale is right (segment and junction CN 1.03 / 1.04× truth), but it still fails (0.545); see Follow-up checks |
| nest_cross_1to1_s502 | RCK | 3 species for 2: a 1.48-copy species carries both species' junctions; 0.955 |
| nest_cross_1to1_s505 | RCK | the same, with 2.84 copies; ecdna2 0.586 |
| ctrl_dup_s703 | RCK | 2 species for 1: the B insertion's observed CN (16.1) is below the truth (19.2) |

The old `nest_cross_1to1_s501`, `s502` and `nest_cross_3sp_s406` subgraphs had the same scale
error (ploidy ≈ 4, depth per copy 7.4–7.9). s501's old failure therefore comes from the subgraph
stage (its old segment CN was 1.81× the truth), not from RCK. `pipeline.tsv`'s stage rule at the
time checked coverage and missing junctions only.

RCK wall time: 4349 s over 22 cases, 2860 s of it in `nest_cross_3sp_s403`.

#### Follow-up checks (bench_nested_v2, 2026-10-27)

These are diagnostic runs on the existing outputs and subgraphs, with no new simulation:

- outputs in `/n/fs/yubi-iw/runs/bench_nested_v2/diag/`;
- Slurm 32005476 and 32005494;
- scripts `rck_diag.sbatch`, `pool_summary.py`, `flow_vs_truth.py`, `scale_ratio.py` and
  `summarize.sh`, kept in `diag/` here.

The plan that follows from them is in `NEXT_STEPS.md`. "Fitted θ" is the junction reads per copy
that `rck_inputs` fits; the bench fixes it at 6.5.

| case | check | result |
|---|---|---|
| nest_cross_3sp_s402 | re-rank the pool, λ = 0.10–0.20 in steps of 0.01 | 3 species, 0.978, from λ = 0.11 on |
| nest_cross_3sp_s403 | subgraph with `--rck-ploidy 2`, then RCK | scale fixed, still 0.545 at every λ; decomposition optimal; best pool solution 0.729 |
| nest_cross_1to1_s502 | RCK with fitted θ 5.11 | **resolved**: 2 species, 0.963 at λ ≥ 0.02 (junctions observed at 0.77× truth with θ 6.5) |
| nest_cross_1to1_s505 | RCK with fitted θ 4.93 | worse, 0.600 (θ 6.5 was right: junctions at 1.00× truth); no good pool solution with θ 6.5 (best = rank 1, 0.765) |
| nest_samehap_3to1_s301 | pool, and RCK with fitted θ 5.09 | the pool holds a 0.973 answer (4 species, 2 of them leftovers ≤ 2.1 copies), ranked below the 3-species rank 1; fitted θ: 0.691 |
| ctrl_dup_s703 | window depth; pool; RCK with fitted θ 4.53 | the B insertion's depth per copy is 12.0 against 12.9–16.4 elsewhere (the simulator's `sigma` 0.1 correlated depth bias); the pool holds a 0.971 answer (the full walk at 16.2 plus a 1.7 leftover), ranked below rank 1; fitted θ: 0.718 |

Reading:

- s402 is a λ choice: any λ ≥ 0.11 resolves it.
- s502 is a junction-scale input error. The fitted θ fixes it here, but it is wrong in s505 (4.93
  against a true ≈ 6.5) and s403 (2.12), so it cannot replace 6.5 as is.
- s301, s703, s505 and s403 (with the right scale) are decomposition failures.
  - The ecDNA flow is integer and fitted to noisy observations, and the decomposition must cover
    it exactly.
  - The residual weight is 1 per copy on a segment ≥ 50 kb (`adj_residual_weight` 1.0 on
    junctions), against `species_penalty` 1.
  - A 1–3 copy dip from the simulated depth bias (±10% per 200 kb, so ±2 copies at 20) therefore
    costs as much as an extra species.
  - Exact-fit alternatives then win: extra leftover species (s703, s301) or wrong merges (s505,
    s403).
  - In s301 and s703 the near-correct answer is in the pool but loses on the ranking. In s505 and
    s403 it is not in the pool at all.

### B.6 First benchmark (2026-10-25; replaced by bench_nested_v2; old λ rule)

Scenario settings of the first benchmark differed from §1.1: per-step haplotype labels, with
`p_cross` 0.2 in the cross and control scenarios. It had 7 scenarios × 2 seeds; nest_samehap_6to1
did not exist yet.

#### Status and reading

Run on ionic: Slurm 31981878/31981889 (validation) and 31981902/31981984 (RCK). 20 of 22 cases
are complete. `nest_cross_3sp_s406` and `nest_cross_1to1_s502` were cancelled at the session's
time limit:

- s406 was still at an 81% MIP gap after 800 s in its λ = 0 decomposition;
- s502 finished λ = 0 (488 s) but not its pool run.

Both need a longer Slurm limit (`--time=05:00:00`) and are excluded from the choice below until
then. The RCK runs live in project space (`runs` → `/n/fs/yubi-iw/runs/bench_nested/rck_runs`; a
5 GB home quota cannot hold them).

The (old) rule below picks 0.5 because of a single case, `nest_samehap_3to1_s301`, which stays
merged and wrong at every λ (0.693 → 0.740). 0.5 is also the edge of the grid. Every other case
gives identical answers from λ = 0.05 to 0.5. **λ = 0.05 is the smallest value that gives the
large gain**:

- mean composite 0.800 → 0.896;
- merges fixed on `nest_cross_2to1`, `nest_samehap_2to1` and demo7 total-only;
- no control or demo case worse than at λ = 0.

The earlier recommendation of 0.1 sits on the same plateau.

#### Results of the first benchmark

**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.5**, the smallest λ on the plateau of best mean composite over all 20 cases (0.5), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.5. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.8002 | 0.8026 | 0.8958 | 0.8958 | 0.8958 | 0.8981 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.5 | merged at 0 | merged at 0.5 |
|---|---|---|---|---|---|
| ctrl_dup | 2 | 0.870 | 0.870 | 0 | 0 |
| ctrl_separate | 2 | 0.670 | 0.670 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 1 | 0.716 | 0.729 | 0 | 0 |
| nest_cross_2to1 | 2 | 0.564 | 0.972 | 2 | 0 |
| nest_cross_3sp | 1 | 0.596 | 0.632 | 1 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 2 | 0.836 | 0.859 | 1 | 1 |

Remaining failures at λ = 0.5 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 6 of 12 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 4 species for 2 true; composites ecdna1=0.962,ecdna2=0.517 |
| nest_cross_3sp_s402 | RCK | merged: ec3 carries ecdna1+ecdna3; composites ecdna1=0.628,ecdna2=0.691,ecdna3=0.577 |
| nest_cross_1to1_s501 | RCK | composites ecdna1=0.893,ecdna2=0.564 |
| ctrl_separate_s601 | subgraph | missing junctions none; min segment coverage 0.0% |
| ctrl_separate_s602 | subgraph | missing junctions none; min segment coverage 0.0% |
| ctrl_dup_s702 | RCK | 1 ecDNA phase block(s) wrong after the flow; imputed edges imputed_adj_intra_chr16_6; composites ecdna1=0.744 |

Subgraph stage: 2 of 12 subgraphs miss a true junction or cover under 80% of a true segment (ctrl_separate_s601, ctrl_separate_s602).

Compute: 754 s of RCK wall time over 20 cases (one λ = 0 run and one pool run each, 4 threads).
