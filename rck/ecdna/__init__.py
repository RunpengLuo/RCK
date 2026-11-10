"""Allele-specific ecDNA reconstruction on top of RCK.

The package turns an ecDNA breakpoint subgraph (the ``nodes.tsv`` / ``edges.tsv`` written by
amplicon-allele-analysis ``subgraph_construct.py``) into haplotype-specific circular
structures in four stages:

1. :mod:`rck.ecdna.subgraph_io` reads the subgraph and normalizes it into a diploid interval
   adjacency graph (DIAG): one segment edge per interval, reference / novel / imputed
   adjacencies between segment extremities.
2. :mod:`rck.ecdna.scale` converts read depth and junction read counts into copy numbers.
3. :mod:`rck.ecdna.ilp_circular` solves a two-layer (ecDNA + chromosomal background) flow MILP
   adapted from RCK's ``OptModelMultiClone``, with strict copy-number balance and no telomeres
   in the ecDNA layer.
4. :mod:`rck.ecdna.decompose` and :mod:`rck.ecdna.walk` decompose the ecDNA flow into the
   fewest circular species and extract ordered, phased walks, with alternative
   decompositions and alternative orderings.

:mod:`rck.ecdna.evaluate` scores the result against simulator truth. The command-line entry
point is ``rck-ecdna`` (:mod:`rck.ecdna_run`). See ``docs/ecDNA.md``.
"""
