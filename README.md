# README #

### What is this repository for? ###

* A family-based phasing algorithm for sequence data
* Main Developer: Mara Battagin mara.battagin@roslin.ed.ac.uk
* Most updated branch: maradevel

### Description ###

AlphaFamSeq is a family-based phasing algorithm for variable-coverage sequence data that prioritises the minimisation of phasing errors and, conditional on this, maximise the proportion of alleles that are phased. The algorithm uses sequence data and at least three generations pedigree of sequenced individuals. The algorithm processes individuals from the oldest to the youngest and iterates until convergence. In the first step, the algorithm calculate allele probabilities based on iterative peeling. In the subsequent steps, the genotypes are phased using heuristics that derive information from the sequence reads of grandparents, parents, and progeny.