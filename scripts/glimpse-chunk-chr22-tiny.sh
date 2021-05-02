#!/usr/bin/env bash
LOPASS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

## prepare genome chunks for imputation
PANEL=data/chr22-tiny/chr22-tiny-genotypes-noNA12878.bcf
OUTDIR=data/chr22-tiny/chunks
WINDOW_SIZE=2000000
BUFFER_SIZE=200000

mkdir -p $OUTDIR
chr=chr22
GLIMPSE_chunk --input $PANEL --window-size $WINDOW_SIZE --buffer-size $BUFFER_SIZE --region $chr \
              --output $OUTDIR/$chr.txt --log $OUTDIR/$chr.log &>> /dev/null
