#!/usr/bin/env python3
import os
import sys
import math
import argparse
import multiprocessing
from collections import namedtuple

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ## I/O
    parser.add_argument("-i", '--input', help="Input BAM/CRAM file", required=True)
    parser.add_argument("-o", '--output', help="output VCF file", required=True)
    ## GATK options
    parser.add_argument("-R", '--ref-fasta', help="Reference FASTA file", required=True)
    parser.add_argument("-L", '--intervals', help="Set of call intervals", required=True)
    parser.add_argument("-J", '--java-opt', help="Java options")
    args=parser.parse_args()

####  final command
## gatk HaplotypeCaller -R /tpo/refs/grch38/assembly/grch38.d1.vd1.fa -ERC GVCF -L chr9:63252863-63542264 -I /mnt/download/glimpse_grch38/bam/NA12878.final-1x.bam -O test.vcf  --min-dangling-branch-length 1 --min-pruning 1 --indel-size-to-eliminate-in-ref-model 0 --heterozygosity 0.1 --indel-heterozygosity 0.1 --native-pair-hmm-threads 1 --base-quality-score-threshold 10
