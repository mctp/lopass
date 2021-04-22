#!/usr/bin/env bash

## prepare genome chunks for imputation
PANEL=$1
OUTDIR=$2
WINDOW_SIZE=2000000
BUFFER_SIZE=200000

mkdir -p $OUTDIR
for chr in chr{1..22}; do
    echo $chr
    GLIMPSE_chunk --input $PANEL --window-size $WINDOW_SIZE --buffer-size $BUFFER_SIZE --region $chr \
                  --output $OUTDIR/$chr.txt --log $OUTDIR/$chr.log &>> /dev/null
done

## chrX par1
GLIMPSE_chunk --input $PANEL --window-size $WINDOW_SIZE --buffer-size $BUFFER_SIZE --region chrX:1-2781479 \
              --output $OUTDIR/chrX_par1.txt --log $OUTDIR/chrX_par1.log &>> /dev/null
LAST_LINE=$(tail -n1 $OUTDIR/chrX_par1.txt  | cut -f1)

## chrX non-par
GLIMPSE_chunk --input $PANEL --window-size $WINDOW_SIZE --buffer-size $BUFFER_SIZE --region chrX:2781480-155701382 \
              --output $OUTDIR/chrX-tmp.txt --log $OUTDIR/chrX.log &>> /dev/null
awk -v "ll=$LAST_LINE" -v IFS='\t' -v OFS='\t' '{ $1 = $1 + ll + 1; print }' $OUTDIR/chrX-tmp.txt > $OUTDIR/chrX.txt
rm $OUTDIR/chrX-tmp.txt
LAST_LINE=$(tail -n1 $OUTDIR/chrX.txt | cut -f1)

## chrX par2
GLIMPSE_chunk --input $PANEL --window-size $WINDOW_SIZE --buffer-size $BUFFER_SIZE --region chrX:155701383-156030895 \
              --output $OUTDIR/chrX_par2-tmp.txt --log $OUTDIR/chrX_par2.log &>> /dev/null
awk -v "ll=$LAST_LINE" -v IFS='\t' -v OFS='\t' '{ $1 = $1 + ll + 1; print }' $OUTDIR/chrX_par2-tmp.txt > $OUTDIR/chrX_par2.txt
rm $OUTDIR/chrX_par2-tmp.txt
