sib_opp_mirror_s17004 (iw_allelic): outer + nested sibling on the opposite allele, mirrored clones (outer in clone 1, sibling in clone 2), clean; the pilot case. RCK MC: both exact, psi per clone exact. CoRAL: 1 of 2 exact, a false species, CN precision 0.35.
truth run: /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_allelic/sib_opp_mirror_s17004

scores (bench_iw_allelic_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  2.00  2.00  0.00  0.00  2.00  2.00  1.00  1.00  1.00  1.00  0.00
                   CoRAL  2.00  3.00  0.00  1.00  1.00  -  0.97  0.75  0.90  0.35  0.29

RCK MC sibling-vs-outer relation (sibling_relation.tsv): opposite 680000 bp, same 0 bp in the truth; 0 bp disagreeing in the answer
