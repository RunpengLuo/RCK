# bench_IW_intpsi: results (2026-11-15)

Branch `yubi_CNP_ionic`. Arrays 32280369 (`rck_bulk`), 32280370 (`rck_mc`), 32280371 (`coral`),
32280372 (`ce`). Scored and collected 2026-11-15 with `bench_iw.py --set iw_intpsi score --jobs 8`
and `collect`. Tables in this folder:
- `bench_iw_intpsi_table.tsv`: one row per (case, method variant), 144 rows;
- `bench_iw_intpsi_pooled.tsv`: means and medians per method variant, per (variant, class) and per
  (variant, class, S); `n` / `n_scored` columns;
- `bench_iw_intpsi_proportions.tsv`: per truth/inferred species proportion pairs.

**Setup.** 2 classes (`simple`; `hard` = nested species) x S in {1, 2, 4} species x 4 replicates =
24 cases, with one shared set of clones and whole-number psi per clone (so the multi-clone input
equals the scoring truth). How it was built: `NEXT_STEPS.md` and both repos' CHANGELOGs.
- `rck_bulk`: the bulk route. `rck_mc`: the multi-clone route given the true clone profiles.
- CoRAL and CycleExtractor (CE) are **decomposition-only baselines**: they decompose `rck_bulk`'s
  exported graph (fitted copy numbers) plus the long-read path constraints. They do not build their
  own graph, so their scores measure decomposition only.
- Each CoRAL / CE answer is scored twice: linear walks kept (`:walks_kept`) and dropped
  (`:walks_dropped`).
- Strict cycles: `strict_rule` `status` for the RCK routes (species vs `low_weight`) and
  `cn-frac:0.1` for CoRAL / CE. Species counts (missed / false) are therefore not on exactly the
  same rule across the two groups.
- `rck_mc`'s missed / false are scored **within phase blocks** (one A/B swap per block, chosen to
  agree best with the truth); the absolute counts are given next to them in table 3. `rck_bulk` is
  scored absolutely. CoRAL / CE carry no alleles, so they have no phased columns.

## Runs and failures

- **All 96 runs finished** (`run_status` ok), all 144 method variants scored, no `score_error.txt`.
- No time-outs, no out-of-memory kills (peak RSS 796 MB of 2 G), no CE crashes (CE crashed on 2
  of 24 graphs on bench_IW).
- **Two solves stopped at their own time limit** and are scored on their incumbent:
  - `rck_bulk` decomposition, hard_S4_s16003: gap 79%;
  - CoRAL, simple_S4_s13002: gap 9%.
- Every `rck_mc` solve was optimal (longest 4.2 s; the 60 s limit was never near).
- CE ran with at most 4 Gurobi threads. Every CoRAL / CE run read `graph_source: bulk+paths` with
  `path_constraints_dropped: 0`.
- **`elapsed_s` is not solver time.** It includes waiting for the 2-session Gurobi licence: 42 of
  96 runs needed more than one attempt (up to 9, 2 min apart). Gurobi reported 5 active sessions
  against a baseline of 2, so the licence was also used outside these arrays. The longest
  `rck_mc` task took 731 s of its 15 min limit, nearly all of it waiting.

## 1. Headline (all 24 cases; means)

ecDNA-only, haplotype-blind, bp-weighted copy-number recall / precision; proportions over strict
cycles (total-variation distance; lower is better), with the `low_weight`-inclusive version; mean
missed / false species per case; best-match blind / phased Jaccard of segments (J_seg) and
adjacencies (J_adj) per truth species; exact cycles per case (truth: 2.33 species per case, 56
in all).

| method variant | ecDNA cn recall | ecDNA cn precision | prop. TV (strict) | prop. TV (+low_weight) | missed | false | blind J_seg | blind J_adj | blind exact cycles | phased J_seg |
|---|---|---|---|---|---|---|---|---|---|---|
| rck_bulk | 0.96 | 0.98 | 0.07 | 0.08 | 0.12 | 0.21 | 0.92 | 0.89 | 1.67 | 0.90 |
| rck_mc | 0.99 | 1.00 | 0.06 | 0.08 | 0.00 | 0.33 | 0.97 | 0.87 | 1.67 | 0.97 |
| coral:walks_kept | 0.85 | 0.48 | 0.23 | 0.25 | 0.46 | 0.29 | 0.66 | 0.80 | 1.00 | – |
| coral:walks_dropped | 0.69 | 0.91 | 0.26 | 0.27 | 0.88 | 0.04 | 0.75 | 0.73 | 1.14 | – |
| ce:walks_kept | 0.96 | 0.40 | 0.33 | 0.36 | 0.25 | 2.00 | 0.79 | 0.80 | 1.04 | – |
| ce:walks_dropped | 0.95 | 0.42 | 0.26 | 0.28 | 0.29 | 0.88 | 0.79 | 0.80 | 1.04 | – |

## 2. Per class and S (means over 4 replicates)

| method variant | class | S | ecDNA cn recall | ecDNA cn precision | prop. TV (strict) | prop. TV (+low_weight) | missed | false | blind J_seg | blind J_adj | blind exact cycles |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rck_bulk | simple | 1 | 0.97 | 0.95 | 0.17 | 0.20 | 0.00 | 0.75 | 0.97 | 0.93 | 0.75 |
| rck_bulk | simple | 2 | 0.95 | 0.99 | 0.03 | 0.05 | 0.00 | 0.25 | 1.00 | 1.00 | 2.00 |
| rck_bulk | simple | 4 | 0.96 | 1.00 | 0.01 | 0.01 | 0.25 | 0.25 | 1.00 | 1.00 | 4.00 |
| rck_bulk | hard | 1 | 1.00 | 0.97 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 1.00 | 1.00 |
| rck_bulk | hard | 2 | 0.96 | 0.99 | 0.02 | 0.04 | 0.00 | 0.00 | 0.92 | 0.86 | 1.50 |
| rck_bulk | hard | 4 | 0.94 | 0.94 | 0.21 | 0.20 | 0.50 | 0.00 | 0.61 | 0.56 | 0.75 |
| rck_mc | simple | 1 | 1.00 | 1.00 | 0.08 | 0.09 | 0.00 | 0.50 | 1.00 | 1.00 | 1.00 |
| rck_mc | simple | 2 | 0.99 | 1.00 | 0.11 | 0.11 | 0.00 | 0.25 | 0.96 | 0.76 | 1.25 |
| rck_mc | simple | 4 | 1.00 | 1.00 | 0.13 | 0.13 | 0.00 | 1.00 | 0.99 | 0.85 | 3.00 |
| rck_mc | hard | 1 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 1.00 | 0.92 | 0.75 |
| rck_mc | hard | 2 | 0.99 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.95 | 0.89 | 1.50 |
| rck_mc | hard | 4 | 0.93 | 1.00 | 0.03 | 0.18 | 0.00 | 0.25 | 0.94 | 0.81 | 2.50 |
| coral:walks_kept | simple | 1 | 0.87 | 0.10 | 0.21 | 0.21 | 0.00 | 0.75 | 0.73 | 1.00 | 0.50 |
| coral:walks_kept | simple | 2 | 0.92 | 0.77 | 0.06 | 0.18 | 0.00 | 0.00 | 0.88 | 1.00 | 1.75 |
| coral:walks_kept | simple | 4 | 0.84 | 0.78 | 0.18 | 0.17 | 0.50 | 0.00 | 0.69 | 0.88 | 2.50 |
| coral:walks_kept | hard | 1 | 0.90 | 0.34 | 0.16 | 0.16 | 0.00 | 0.75 | 0.75 | 1.00 | 0.50 |
| coral:walks_kept | hard | 2 | 0.92 | 0.48 | 0.27 | 0.28 | 0.00 | 0.25 | 0.57 | 0.64 | 0.50 |
| coral:walks_kept | hard | 4 | 0.68 | 0.39 | 0.49 | 0.46 | 2.25 | 0.00 | 0.31 | 0.28 | 0.25 |
| coral:walks_dropped | simple | 1 | 0.65 | 0.67 | 0.21 | 0.21 | 0.25 | 0.25 | 0.97 | 0.90 | 0.67 |
| coral:walks_dropped | simple | 2 | 0.88 | 1.00 | 0.10 | 0.21 | 0.25 | 0.00 | 0.87 | 0.88 | 1.75 |
| coral:walks_dropped | simple | 4 | 0.77 | 1.00 | 0.31 | 0.31 | 1.25 | 0.00 | 0.69 | 0.69 | 2.50 |
| coral:walks_dropped | hard | 1 | 0.62 | 0.99 | 0.12 | 0.12 | 0.25 | 0.00 | 1.00 | 1.00 | 0.67 |
| coral:walks_dropped | hard | 2 | 0.74 | 1.00 | 0.37 | 0.38 | 0.50 | 0.00 | 0.57 | 0.60 | 0.50 |
| coral:walks_dropped | hard | 4 | 0.49 | 0.74 | 0.43 | 0.41 | 2.75 | 0.00 | 0.41 | 0.29 | 0.33 |
| ce:walks_kept | simple | 1 | 1.00 | 0.11 | 0.47 | 0.47 | 0.00 | 2.75 | 0.97 | 0.93 | 0.75 |
| ce:walks_kept | simple | 2 | 0.99 | 0.50 | 0.10 | 0.11 | 0.00 | 1.00 | 0.84 | 1.00 | 1.50 |
| ce:walks_kept | simple | 4 | 0.97 | 0.37 | 0.15 | 0.19 | 0.25 | 3.00 | 0.74 | 0.94 | 2.25 |
| ce:walks_kept | hard | 1 | 1.00 | 0.34 | 0.36 | 0.36 | 0.00 | 2.00 | 1.00 | 1.00 | 1.00 |
| ce:walks_kept | hard | 2 | 0.97 | 0.48 | 0.47 | 0.50 | 0.00 | 2.50 | 0.75 | 0.65 | 0.50 |
| ce:walks_kept | hard | 4 | 0.81 | 0.62 | 0.45 | 0.55 | 1.25 | 0.75 | 0.45 | 0.29 | 0.25 |
| ce:walks_dropped | simple | 1 | 0.99 | 0.11 | 0.37 | 0.37 | 0.00 | 1.75 | 0.97 | 0.93 | 0.75 |
| ce:walks_dropped | simple | 2 | 0.98 | 0.50 | 0.05 | 0.06 | 0.00 | 0.25 | 0.84 | 1.00 | 1.50 |
| ce:walks_dropped | simple | 4 | 0.97 | 0.39 | 0.07 | 0.07 | 0.25 | 0.75 | 0.74 | 0.94 | 2.25 |
| ce:walks_dropped | hard | 1 | 0.99 | 0.41 | 0.23 | 0.23 | 0.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| ce:walks_dropped | hard | 2 | 0.95 | 0.48 | 0.43 | 0.45 | 0.00 | 1.25 | 0.75 | 0.65 | 0.50 |
| ce:walks_dropped | hard | 4 | 0.81 | 0.62 | 0.42 | 0.53 | 1.50 | 0.25 | 0.44 | 0.28 | 0.25 |

## 3. RCK routes: total copy number and species counts

Total copy-number recall / precision in the cell-averaged tumour frame, normal cells included,
over the truth ecDNA bases first and over the whole union subgraph second. "missed / false" for
`rck_mc` are within phase blocks; "(absolute)" is the same pairing without the per-block swap.

| method | class | S | total cn recall (ecDNA bases) | total cn precision (ecDNA bases) | total cn recall (subgraph) | total cn precision (subgraph) | missed | false | missed (absolute) | false (absolute) | phased exact cycles |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rck_bulk | simple | 1 | 0.98 | 0.95 | 0.97 | 0.99 | 0.00 | 0.75 | – | – | 0.75 |
| rck_bulk | simple | 2 | 0.97 | 0.99 | 0.98 | 0.99 | 0.00 | 0.25 | – | – | 2.00 |
| rck_bulk | simple | 4 | 0.97 | 0.99 | 0.96 | 1.00 | 0.25 | 0.25 | – | – | 3.75 |
| rck_bulk | hard | 1 | 1.00 | 0.96 | 1.00 | 0.99 | 0.00 | 0.00 | – | – | 1.00 |
| rck_bulk | hard | 2 | 0.99 | 0.95 | 1.00 | 0.94 | 0.00 | 0.00 | – | – | 1.50 |
| rck_bulk | hard | 4 | 0.97 | 0.94 | 0.99 | 0.92 | 0.50 | 0.00 | – | – | 0.75 |
| rck_mc | simple | 1 | 1.00 | 1.00 | 1.00 | 1.00 | 0.00 | 0.50 | 0.00 | 0.50 | 1.00 |
| rck_mc | simple | 2 | 1.00 | 1.00 | 1.00 | 1.00 | 0.00 | 0.25 | 0.50 | 0.75 | 1.25 |
| rck_mc | simple | 4 | 1.00 | 1.00 | 1.00 | 1.00 | 0.00 | 1.00 | 1.50 | 2.50 | 3.00 |
| rck_mc | hard | 1 | 1.00 | 1.00 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.75 |
| rck_mc | hard | 2 | 1.00 | 1.00 | 1.00 | 1.00 | 0.00 | 0.00 | 0.00 | 0.00 | 1.50 |
| rck_mc | hard | 4 | 0.94 | 1.00 | 0.99 | 1.00 | 0.00 | 0.25 | 0.75 | 1.00 | 2.50 |

## Findings

**ecDNA copy numbers.**
- `rck_mc` is near exact: recall 0.99, precision 1.00. Every case is at recall >= 0.97 except
  hard_S4_s16003 (0.71; below).
- `rck_bulk`: 0.96 / 0.98. It loses most on hard S=4 (precision 0.94, J_seg 0.61).
- Inferred ecDNA mass / truth mass, median per case: `rck_bulk` 0.99, `rck_mc` 1.00,
  CoRAL 3.30 (walks kept) and 0.77 (walks dropped), CE 3.2.

**Species and proportions.**
- `rck_mc` misses no species in any case (within phase blocks). Its 8 false species come from 5
  cases: simple_S1_s11004 (2), simple_S2_s12004 (1), simple_S4_s13002 (3), simple_S4_s13004 (1)
  and hard_S4_s16003 (1). In the four simple ones copy-number recall and precision are 1.00, so
  the extra cycles redistribute copies without changing copy numbers.
  - *Checked example:* simple_S1_s11004. Clone 2's 6 copies of the single species come back as
    1 copy of the right cycle plus two smaller cycles over the same segments. That gives 2 false
    species and proportion TV 0.32.
  - simple_S2_s12004 and simple_S4_s13002 have proportion TV 0.40 each; their cycles were not
    inspected.
- On hard (nested) cases `rck_mc` keeps J_seg 0.94 at S=4, against 0.61 for `rck_bulk` and 0.31–0.45
  for the baselines.
- Absolute vs within-block counts for `rck_mc`: 11 missed / 19 false absolute against 0 / 8 within
  blocks over the 24 cases. The difference is A/B labels that are identified only within a phase
  block.

**`rck_mc` on hard_S4_s16003.** Copy-number recall 0.71 and total copy-number recall over the ecDNA
bases 0.76, with blind J_seg 0.76 and 7 cycles for 4 species. The per-clone psi are within 1 of the
truth for every species. This is the case where, on bench_IW, the multi-clone model already moved
its given profiles. Not investigated further here; worth a look before the write-up of the
multi-clone route.

**Total copy number (RCK routes).** Over the ecDNA bases `rck_bulk` gets 0.98 / 0.96 and `rck_mc`
0.99 / 1.00. Over the whole subgraph both are 0.97–1.00. The purity effect expected in
NEXT_STEPS step 3 turned out small on this bench.

**CoRAL / CE (decomposition-only, on the bulk graph).**
- Even given `rck_bulk`'s graph and the read paths, both do worse than `rck_bulk`'s own
  decomposition on every headline metric, except CE's copy-number recall (equal, 0.96) and
  CoRAL's false species with walks dropped (0.04 against 0.21, at the cost of 0.88 missed).
- *CoRAL:* walks kept: recall 0.85, precision 0.48. Walks dropped: precision rises to 0.91, but
  recall falls to 0.69 and it misses 0.88 species per case (2.75 at hard S=4).
- *CE:* high recall (0.95–0.96) but precision about 0.4 either way. Dropping linear walks barely
  changes copy numbers because CE returns the background mostly as closed cycles.
- In 4 cases (simple_S1_s11002, simple_S1_s11004, simple_S2_s12004, hard_S1_s14004) the inferred
  ecDNA mass of CE, and of CoRAL with walks kept, is over 10x the truth: background returned as
  ecDNA. Neither walk variant removes this (the known limit in NEXT_STEPS decision 8).
- Both degrade most on hard S=4 (recall 0.49–0.81, J_seg 0.31–0.45).
