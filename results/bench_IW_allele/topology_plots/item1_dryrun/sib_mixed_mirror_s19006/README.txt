sib_mixed_mirror_s19006 (iw_allelic); truth run /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_allelic/sib_mixed_mirror_s19006

--- topology_compare.png (make_compare_plots.py --topology-case) ---
topology_compare.png: one genome axis; truth species, RCK MC strict species (haplotypes after the scorer's one A/B swap per phase block), CoRAL strict cycles (coral:walks_kept, decomposition-only); -> = the table scorer's pairing.

scores (bench_iw_allelic_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  2.00  2.00  0.00  0.00  1.00  1.00  1.00  0.83  0.96  1.00  0.01
                   CoRAL  2.00  1.00  1.00  0.00  1.00  -  0.69  0.50  0.87  1.00  0.26
