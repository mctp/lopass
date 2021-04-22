#!/usr/bin/env bash

## re-codes GLIMPSE genetic maps for hg19/hg38 chromosome naming
MAPDIR=$1
OUTDIR=$2

mkdir -p $OUTDIR
for map in $MAPDIR/*; do
    zcat $map | awk -v OFS='\t' '{ if (NR!=1) { $2 = "chr"$2; print } else { print } }' | gzip > $OUTDIR/$(basename $map)
done
