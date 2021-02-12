# gvcf-call
Call variants from gVCF file.

Given a reference panel (normalized indels, multi-allelics decomposed) and a single-sample gVCF file
genotypes the positions in the reference panel given the evidence in the gVCF file. The gVCF file should
not be decomposed.

# Assumptions

Panel variants are annotated with genotype calls (GT) and likelihoods (PL) from sample gVCF file
if the variants 'match'.

1. Panel variant overlaps a reference block. Every panel variant within a reference block will have a homozygous reference genotype and the same PL based on the catch-all <NON_REF> ALT.
2. Panel variant matches a sample variant.

Two variants are defined as matching if and only if:
- they match on POS and REF
- The panel ALT is the most-likely ALT in the sample variant

# Requirements

- Python3
- cyvcf2
