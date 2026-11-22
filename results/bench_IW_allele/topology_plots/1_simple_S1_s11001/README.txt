simple_S1_s11001: 1 species, one allele
class simple, S = 1; truth run /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_intpsi/simple_S1_s11001

--- topology_compare.png (make_compare_plots.py --topology-case) ---
topology_compare.png: one genome axis; truth species, RCK MC strict species (haplotypes after the scorer's one A/B swap per phase block), CoRAL strict cycles (coral:walks_kept, decomposition-only); -> = the table scorer's pairing.

scores (../bench_IW/bench_iw_intpsi_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  1.00  1.00  0.00  0.00  1.00  1.00  1.00  1.00  1.00  1.00  0.00
                   CoRAL  1.00  2.00  0.00  1.00  1.00  -  1.00  1.00  1.00  0.19  0.17
