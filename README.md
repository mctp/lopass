# lopass

Pipeline for low-coverage genotyping, imputation, phasing.

*Step 1*: gVCF Variant Calling (GATK4)

lopass-gatk.py A gVCF file is created calling variation within a BAM file. 

*Step 2*: Extraction of genotype likelihoods (PLs) for panel variants

lopass-genotype.py A VCF file is created with PLs for all variants in a reference panel.

Given a reference panel (normalized indels, multi-allelics decomposed) and a single-sample gVCF file
genotypes the positions in the reference panel given the evidence in the gVCF file. The gVCF file should
not be decomposed.

*Step 3*: Imputation and Phasing

lopass-glimpse.py Imputation and phasing is performed on the resulting VCF file.

## Installation

Lopass requires a Linux system with configured to compile `htslib` and `boost`, `Python3` and a number of third-party optionally managed by conda.

### Dependencies

- conda (system installation)
- boost (bundled)
- htslib (bundled)
- Python3 (managed via conda)
- cyvcf2 (managed via conda)
- gatk4 (managed via conda)

### Installation

1. Install conda
2. Execute `build.sh` script, this will compile `GLIMPSE`
3. Execute `setup-conda.sh` this will create a `lopass` environment. Alternatively,  install `Python3`, `cyvcf2` and `GATK4` manually.

Before running lopass setup paths and activate the conda environment:

```
source env.sh
```

## Details

Panel variants are annotated with genotype calls (GT) and likelihoods (PL) from sample gVCF file. A panel variant may overlap a gVCF reference block or variant site.

1. Panel variant overlaps a reference block. Every panel variant within a reference block will have a homozygous reference genotype and the same PL based on the catch-all <NON_REF> ALT.
2. Panel variant matches a sample variant. Two variants are defined as matching if and only if:
- they match on POS and REF
- The panel ALT is the most-likely ALT in the sample variant
- The sample genotype is homozygous reference, heterozygous, or homozygous alternative for the panel ALT

In all other cases, the site is returned as unobserved.

### Sex chromosomes

Chromosome X is fully supported with 'correct' handling of XX and XY individuals and PAR1 and PAR2 regions.

Chromosome Y is not supported.

