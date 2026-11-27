sib_opp_shared_s18004 (iw_allelic): two-column print-size topology comparison (make_compare_plots_ce_2col.py, 2026-11-26)

topology_compare_ce_2col.pdf (vector) / topology_compare_ce_2col.png (300 dpi): the walks of topology_compare_ce.png (same data, pairing and drawing; see README.txt) re-laid out for a printed page.

layout: two columns sharing one genome x-axis (same regions and scale) and one row height (0.125 in per step); left = Truth (top), RCK MC; right = CoRAL (top), CE; columns start at the same height. One row per step in walk order (reference-contiguous steps merged), numbered on the left; arrows = strand; connectors in traversal order, dotted = closing join. Label above each walk: name -> the truth species the table scorer pairs it with, bulk copies; psi per clone (Truth, RCK MC) or cycle / linear walk (CoRAL, CE). RCK MC haplotypes after the scorer's one A/B swap per phase block. One shared legend at the bottom; no figure title (goes in the caption).
size: 6.50 x 8.50 in (drawn at the final printed size, full 6.5 in text width; column width 2.97 in).
fonts (pt): panel titles 8 bold, walk labels 7, step numbers 6.5, axis ticks 6.5, axis label and chromosome names 7, legend 7; smallest 6.5. Fonts embedded as TrueType (pdf.fonttype 42).
walks drawn: Truth 2, RCK MC 2, CoRAL 1, CE 2.

pairings drawn:
  RCK MC: ec1 -> ecdna2, ec2 -> ecdna1
  CoRAL: cycle1 -> ecdna1
  CE: cycle1 -> ecdna1, cycle2 -> ecdna2

inputs (read only):
  truth: /u/ym3292/amplicon-allele-analysis-cnp/simulations/runs/bench_IW_allelic/sib_opp_shared_s18004
  RCK MC answer: /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/rck_mc/score/answer/sib_opp_shared_s18004_rck_mc_cycles.txt
  RCK MC pairing: /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/rck_mc/score/aa_eval/aa_species_metrics.tsv
  RCK MC psi: /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/rck_mc/decompositions.tsv
  phase blocks: /u/ym3292/amplicon-allele-analysis-cnp/simulations/subgraphs/bench_IW_allelic_cnp/sib_opp_shared_s18004/union/rck_input/rck.fragments.tsv
  clone proportions: /u/ym3292/amplicon-allele-analysis-cnp/simulations/subgraphs/bench_IW_allelic_cnp/sib_opp_shared_s18004/union/rck_input/clone_proportions.tsv
  CoRAL answer (coral:walks_kept): /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/coral/amplicon1/score/answer/coral_amplicon1_cycles.txt
  CoRAL pairing: /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/coral/amplicon1/score/aa_eval/aa_species_metrics.tsv
  CE answer (ce:walks_kept): /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/ce/amplicon1/score/answer/ce_amplicon1_cycles.txt
  CE pairing: /n/fs/yubi-iw/runs/bench_IW_allelic/methods/sib_opp_shared_s18004/ce/amplicon1/score/aa_eval/aa_species_metrics.tsv
