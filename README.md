1;5202;0c11;rgb:2323/2727/292911;rgb:2323/2727/2929# lopass

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

### Building

1. Install conda
2. Execute `build.sh` script, this will compile `GLIMPSE`
3. Execute `setup-conda.sh` this will create a `lopass` environment. Alternatively,  install `Python3`, `cyvcf2` and `GATK4` manually.

### Running

Before running lopass setup paths and activate the conda environment:

```
source env.sh
```

Lopass is currently only compatible with GRCh38-style chromosome names we need to re-code the genetic maps included in GLIMPSE

```
scripts/glimpse-recode-maps.sh GLIMPSE/maps/genetic_maps.b38 data/hg38/maps
```

Imputation panels need to be chunked by GLIMPSE prior to running:

```
scripts/glimpse-chunk-panel.sh <panel> data/hg38/chunks
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

### Example

Step 1:
```
./lopass-gatk.py -t1 -p4 -L chr22:1-18709565:2:2 -R grch38_full_analysis_set_plus_decoy_hla.fa --sex XX -i data/chr22-tiny/NA12878-10M-chr22-tiny.bam -o data/chr22-tiny/NA12878-10M-chr22-tiny
```

Step 2:
```
./lopass-genotype.py data/chr22-tiny/kg-hg38-nygc2020-chr22-tiny-sites.bcf data/chr22-tiny/NA12878-10M-chr22-tiny.vcf.gz -r chr22:1-18709565 | bcftools view -Oz -o data/chr22-tiny/NA12878-10M-chr22-tiny-call.vcf.gz
tabix data/chr22-tiny/NA12878-10M-chr22-tiny-call.vcf.gz
```

Step 3:
```
./lopass-glimpse.py -c data/chr22-tiny/chunks -m data/chr22-tiny/maps data/chr22-tiny/kg-hg38-nygc2020-chr22-tiny-genotypes-noNA12878.bcf  data/chr22-tiny/NA12878-10M-chr22-tiny-call.vcf.gz -o data/chr22-tiny/NA12878-10M-chr22-tiny-glimpse -r chr22 -d
```
