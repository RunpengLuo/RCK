<!-- results:start -->
**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.2**, the smallest λ on the plateau of best mean composite over all 22 cases (0.2, 0.5; a λ whose shortfall to the best comes from fewer than 2 cases is on it), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.2. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.8253 | 0.8579 | 0.9055 | 0.9055 | 0.9262 | 0.9283 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |
| cases gaining from λ to the best λ | demo7_totalonly, nest_cross_2to1_s104, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | demo7_totalonly, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_samehap_3to1_s301 | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.2 | merged at 0 | merged at 0.2 |
|---|---|---|---|---|---|
| ctrl_dup | 2 | 0.855 | 0.855 | 0 | 0 |
| ctrl_separate | 2 | 0.979 | 0.979 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 2 | 0.860 | 0.860 | 2 | 2 |
| nest_cross_2to1 | 2 | 0.609 | 0.967 | 1 | 0 |
| nest_cross_3sp | 2 | 0.533 | 0.761 | 2 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 2 | 0.836 | 0.836 | 1 | 1 |

Remaining failures at λ = 0.2 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 5 of 14 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.964,ecdna2=0.423 |
| nest_cross_3sp_s403 | RCK | merged: ec1 carries ecdna1+ecdna2; ec2 carries ecdna1+ecdna3; composites ecdna1=0.519,ecdna2=0.593,ecdna3=0.523 |
| nest_cross_1to1_s502 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.967,ecdna2=0.943 |
| nest_cross_1to1_s505 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.944,ecdna2=0.586 |
| ctrl_dup_s703 | RCK | 2 species for 1 true; composites ecdna1=0.712 |

Scale: 1 of 14 cases have depth per copy more than 30% from the median or segment / junction scales apart by more than 1.5x (nest_cross_2to1_s104).
Subgraph stage: 0 of 14 subgraphs miss a true junction or cover under 80% of a true segment.

Compute: 2762 s of RCK wall time over 22 cases (one λ = 0 run and one pool run each, 4 threads).
<!-- results:end -->
