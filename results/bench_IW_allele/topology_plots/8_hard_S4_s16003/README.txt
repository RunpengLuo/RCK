hard_S4_s16003: 4 nested species, both alleles amplified (70%); the hardest case, multi-clone's weakest
class hard, S = 4; truth run /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_intpsi/hard_S4_s16003

--- topology_compare.png (make_compare_plots.py --topology-case) ---
topology_compare.png: one genome axis; truth species, RCK MC strict species (haplotypes after the scorer's one A/B swap per phase block), CoRAL strict cycles (coral:walks_kept, decomposition-only); -> = the table scorer's pairing.

scores (../bench_IW/bench_iw_intpsi_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  4.00  5.00  0.00  1.00  2.00  2.00  0.76  0.66  0.71  0.99  0.12
                   CoRAL  4.00  1.00  3.00  0.00  0.00  -  0.00  0.24  0.28  0.03  0.85
