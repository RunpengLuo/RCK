hard_S2_s15001 (iw_intpsi): truth vs RCK MC vs CoRAL vs CE topology comparison (make_compare_plots_ce.py, 2026-11-26)

topology_compare_ce.png (150 dpi) / .pdf: one genome axis, one row per step in walk order; panels truth, RCK MC (strict species, haplotypes after the scorer's one A/B swap per phase block), CoRAL (strict cycles, coral:walks_kept, decomposition-only, grey = no alleles), CE (CycleExtractor 1.0.0, strict cycles, ce:walks_kept, decomposition-only, grey). -> = the table scorer's pairing (strict mode of each run's aa_species_metrics.tsv). Same drawing as make_compare_plots.topology_compare.

pairings drawn:
  RCK MC: ec1 -> ecdna1, ec2 -> ecdna2
  CoRAL: cycle1 -> ecdna2, cycle2 -> ecdna1
  CE: cycle2 -> ecdna2, cycle4 -> ecdna1
steps over 5 Mb clipped to the axis: CoRAL cycle3

inputs (read only):
  truth: /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_intpsi/hard_S2_s15001
  RCK MC answer: /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/rck_mc/score/answer/hard_S2_s15001_rck_mc_cycles.txt
  RCK MC pairing: /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/rck_mc/score/aa_eval/aa_species_metrics.tsv
  RCK MC psi: /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/rck_mc/decompositions.tsv
  phase blocks: /u/ym3292/amplicon-allele-analysis-cnp/simulations/subgraphs/bench_IW_intpsi_cnp/hard_S2_s15001/union/rck_input/rck.fragments.tsv
  clone proportions: /u/ym3292/amplicon-allele-analysis-cnp/simulations/subgraphs/bench_IW_intpsi_cnp/hard_S2_s15001/union/rck_input/clone_proportions.tsv
  CoRAL answer (coral:walks_kept): /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/coral/amplicon1/score/answer/coral_amplicon1_cycles.txt
  CoRAL pairing: /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/coral/amplicon1/score/aa_eval/aa_species_metrics.tsv
  CE answer (ce:walks_kept): /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/ce/amplicon1/score/answer/ce_amplicon1_cycles.txt
  CE pairing: /n/fs/yubi-iw/runs/bench_IW_intpsi/methods/hard_S2_s15001/ce/amplicon1/score/aa_eval/aa_species_metrics.tsv

scores (../bench_IW/bench_iw_intpsi_table.tsv; phased = within phase blocks, RCK MC only):
                  method  n_truth  n_strict  n_missed  n_false  blind_exact_cycles  phased_exact_cycles  blind_jseg_best_mean  blind_jadj_best_mean  cn_recall  cn_precision  proportion_tv
                  RCK MC  2.00  2.00  0.00  0.00  2.00  2.00  1.00  1.00  1.00  1.00  0.00
                   CoRAL  2.00  3.00  0.00  1.00  0.00  -  0.59  0.50  0.97  0.36  0.18
                      CE  2.00  4.00  0.00  2.00  0.00  -  0.74  0.58  1.00  0.35  0.72
