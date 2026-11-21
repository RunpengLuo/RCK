sib_opp_shared_s18004 (iw_allelic): outer + nested sibling on the opposite allele, equal clone weights (no clone-wise allele signal in the input), clean. RCK MC: both species exact. CoRAL: one cycle over both, sibling missed.
truth run: /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_allelic/sib_opp_shared_s18004

scores (bench_iw_allelic_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  2.00  2.00  0.00  0.00  2.00  2.00  1.00  1.00  1.00  1.00  0.00
                   CoRAL  2.00  1.00  1.00  0.00  0.00  -  0.42  0.50  0.84  1.00  0.72

RCK MC sibling-vs-outer relation (sibling_relation.tsv): opposite 413000 bp, same 0 bp in the truth; 0 bp disagreeing in the answer
