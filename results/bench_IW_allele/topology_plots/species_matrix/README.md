# Within-phase-block recovery matrices and J_seg / J_adj species matrices (make_block_figures.py)

Written 2026-11-21. The runs' own `recovery_matrix.tsv` (and `rck_mc_recovery_matrix.png`) pair truth and inferred
species and score haplotypes in the absolute frame. A phase block's A/B orientation is not identifiable from the data,
so the benchmark tables score within phase blocks: one A/B swap per block (`evaluate.block_swaps`), chosen from the
haplotype-blind pairing of the strict species, then pairing again with haplotypes (`aa_eval.score --phase-blocks`).
Everything here uses that frame and pairing. Run folders are only read.

- `<set>/<case>/rck_mc_recovery_matrix_block.{tsv,png}`: the recovery matrix within phase blocks. Its composites equal
  the benchmark table's per species and on average on all 36 cases (`recovery_check.tsv`), with the same pairing.
- `<set>/<case>/rck_mc_species_matrix_block.{png,pdf}`: J_seg | J_adj, RCK MC within phase blocks (`aa_eval_block/`).
- `<set>/<case>/rck_mc_species_matrix_blind.{png,pdf}`: the same, haplotype-blind.
- `<set>/<case>/coral_species_matrix_blind.{png,pdf}`: CoRAL, haplotype-blind (the 8 `../compare/` cases only).
- `jmatrix_summary.tsv`: one row per (set, case, method / frame, truth species): paired species, paired J_seg / J_adj,
  row maxima and whether they differ from the pair, `exact_cycle`, `is_1_1` (both >= 0.995, i.e. 1.00 as printed).
- `jmatrix_not_1_1.tsv`: every RCK MC within-block paired cell that is not (1, 1), with the truth novel adjacencies
  the paired walk misses / adds and a reason.
- `recovery_check.tsv`: new mean composite vs the table's `mean_composite`, per case.

Sets: `iw_allelic` (12 allelic-sibling cases), `iw_intpsi` (24 bench_IW_intpsi cases). Redo (local, ~25 s, no Gurobi):

    python make_block_figures.py run --jobs 8 && python make_block_figures.py reasons && python make_block_figures.py compare
