sib_mixed_mirror_s19012 (iw_allelic): sibling on the opposite allele on some shared intervals and the same allele on others, mirrored clones, clean. RCK MC: both exact. CoRAL: sibling missed.
truth run: /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_allelic/sib_mixed_mirror_s19012

scores (bench_iw_allelic_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  2.00  2.00  0.00  0.00  2.00  2.00  1.00  1.00  1.00  1.00  0.00
                   CoRAL  2.00  1.00  1.00  0.00  1.00  -  0.77  0.50  0.88  1.00  0.16

RCK MC sibling-vs-outer relation (sibling_relation.tsv): opposite 167000 bp, same 190000 bp in the truth; 0 bp disagreeing in the answer
