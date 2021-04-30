#!/usr/bin/env python3
import os
import sys
import math
import cyvcf2
import argparse
from collections import namedtuple

GT2PL = {
    (0,):0,
    (1,):1,
    (2,):2,
    (3,):3,
    (4,):4,
    (5,):5,
    (6,):6,
    (0,0):0,
    (0,1):1,
    (1,1):2,
    (0,2):3,
    (1,2):4,
    (2,2):5,
    (0,3):6,
    (1,3):7,
    (2,3):8,
    (3,3):9,
    (0,4):10,
    (1,4):11,
    (2,4):12,
    (3,4):13,
    (4,4):14,
    (0,5):15,
    (1,5):16,
    (2,5):17,
    (3,5):18,
    (4,5):19,
    (5,5):20,
    (0,6):21,
    (1,6):22,
    (2,6):23,
    (3,6):24,
    (4,6):25,
    (5,6):26,
    (6,6):27,
}

GT2GT = {
    (0,):0,
    (1,):1,
    (2,):1,
    (3,):1,
    (4,):1,
    (5,):1,
    (6,):1,
    (0,0):(0,0),
    (0,1):(0,1),
    (1,1):(1,1),
    (0,2):(0,1),
    (2,2):(1,1),
    (0,3):(0,1),
    (3,3):(1,1),
    (0,4):(0,1),
    (4,4):(1,1),
    (0,5):(0,1),
    (5,5):(1,1),
    (0,6):(0,1),
    (6,6):(1,1),
}

VREC = namedtuple("Variant", ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
VCFREC = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n"

def make_header(vcf1, vcf2):
    head1_ln = vcf1.raw_header.strip().split("\n")
    head2_ln = vcf2.raw_header.strip().split("\n")
    head_ln = []
    head_ln.append('##fileformat=VCFv4.1')
    head_ln.extend([l for l in head1_ln if l.startswith("##contig")])
    head_ln.extend([
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GC,Number=1,Type=String,Description="Genotype call {ref,err,var,pra,npa,mpl}">',
        '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods rounded to integer">',
    ])
    head_ln.append(head2_ln[-1])
    head = "\n".join(head_ln) + "\n"
    return head
    
def process_query(vcf1, vcf2, query=None, gt_missing0=False, pl_missing0=False, pl_max=math.inf):

    if gt_missing0:
        GT_MISSING_HAPLOID = "0"
        GT_MISSING_DIPLOID = "0/0"
    else:
        GT_MISSING_HAPLOID = "."
        GT_MISSING_DIPLOID = "./."

    if pl_missing0:
        PL_MISSING_HAPLOID = "0,0"
        PL_MISSING_DIPLOID = "0,0,0"
    else:
        PL_MISSING_HAPLOID = "."
        PL_MISSING_DIPLOID = "."
    
    ## initialize parsers
    gen1 = vcf1(query)
    gen2 = vcf2(query)
    
    ## in a gvcf file every base is covered by at least one reference block or variant
    ## for each variant in vcf1 find the last block/variant in vcf2 which encloses it
    v2 = next(gen2)
    v2_next = next(gen2)
    gvcf_end = False
    for v1 in gen1:
        
        ## find v2 which matches v1
        while not gvcf_end and (v1.start >= v2_next.start):
            v2 = v2_next
            try:
                v2_next = next(gen2)
            except StopIteration:
                v2_next = None
                gvcf_end = True
        assert (v1.start>=v2.start & v1.start<=v2.end), "variant in panel did not overlap a variant/block in gVCF"

        ## missing
        missing = False
        ## get genotype likelihoods
        if not v2.format("PL") is None:
            # sample variant has PLs
            pl = v2.format("PL")[0] # genotype likelihoods per allele
            haploid = len(v2.genotypes[0]) == 2 # whether variant is haploid
            if v2.ALT[0]=="<NON_REF>":
                # panel variant overlaps a sample reference block
                if haploid:
                    GT = "0"
                    PL = "%s,%s" % (pl[GT2PL[(0,)]], pl[GT2PL[(1,)]])
                else:
                    GT = "0/0"
                    PL = "%s,%s,%s" % (pl[GT2PL[(0,0)]], pl[GT2PL[(0,1)]], pl[GT2PL[(1,1)]])
                match = "ref"
                qual = "."
            elif (v1.POS==v2.POS) and (v1.REF==v2.REF) and (v1.ALT[0] in v2.ALT):
                # panel variant overlaps a sample variant
                alt = v2.ALT.index(v1.ALT[0]) + 1 # index of the panel alternate allele
                gt = tuple(v2.genotype.alleles(0)) # genotype of sample
                ad = v2.format("AD")[0] # allele depths
                if (gt in ((0,0),(0,alt),(alt,alt),(0,),(alt,)) and gt in GT2GT):
                    # variant genotype has only ref and panel alt alleles
                    gt_ = GT2GT[gt]
                    if not ((0 in gt_ and ad[0]==0) or (1 in gt_ and ad[alt]==0)):
                        # variant genotype matches ADs
                        if haploid:
                            pl0 = pl[GT2PL[(0,)]]
                            pl1 = pl[GT2PL[(alt,)]]
                            GT = "%s" % gt_
                            PL = "%s,%s" % (pl0, pl1)
                        else:
                            pl00 = pl[GT2PL[(0,0)]]
                            pl01 = pl[GT2PL[(0,1)]]
                            pl11 = pl[GT2PL[(1,1)]]
                            GT = "%s/%s" % gt_
                            PL = "%s,%s,%s" % (pl00, pl01, pl11)
                        match = "var"
                        qual = "%.2f" % v2.QUAL
                    else:
                        missing = 'err'
                else:
                    # variant genotype includes non-panel alt alleles
                    # variant genotype has multiple non-reference alts (e.g. variant 1/2)
                    # variant genotype does not match panel alt (e.g. variant 0/1 panel 0/2)
                    # alt number > GT2PL keys (max 6)
                    missing = 'npa'
            else:
                # overlapping variants but pos, ref, alt mismatch
                missing = 'pra'
        else:
            # malformed variant with missing PLs
            missing = 'mpl'
                
        if missing:
            if haploid:
                GT = GT_MISSING_HAPLOID
                PL = PL_MISSING_HAPLOID
            else:
                GT = GT_MISSING_DIPLOID
                PL = PL_MISSING_DIPLOID
            match = missing
            qual = "."

        vrec = VREC(
            CHROM=v1.CHROM,
            POS=v1.POS,
            ID=v1.ID if v1.ID else ".",
            REF=v1.REF,
            ALT=v1.ALT[0],
            QUAL=qual,
            FILTER=".",
            INFO=".",
            FORMAT="GT:GC:PL",
            SAMPLE="%s:%s:%s" % (GT, match, PL)
        )
        
        yield vrec


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-o", '--out', help="output file")
    parser.add_argument("-H", '--no-header', help="suppress the header in VCF output")
    parser.add_argument("-r", '--regions', help="query regions in <chr>:<beg>-<end> format comma delimited")
    parser.add_argument('--gt_missing0', dest='gt_missing0', action='store_true', help="Use 0 or 0/0 for missing GTs")
    parser.add_argument('--pl_missing0', dest='pl_missing0', action='store_true', help="Use 0-array for missing PLs")
    parser.add_argument("panel", help="VCF file of sites in reference panel")
    parser.add_argument("gvcf", help="gVCF file of genotyped sample")
    args=parser.parse_args()

    vcf1 = cyvcf2.VCF(args.panel)
    vcf2 = cyvcf2.VCF(args.gvcf)

    if args.regions:
        queries = args.regions.split(",")
    else:
        queries = vcf1.seqnames

    if args.out:
        fh = open(args.out, 'wb')
    else:
        fh = sys.stdout

    if not args.no_header:
        fh.write(make_header(vcf1, vcf2))

    for query in queries:
        recs = process_query(vcf1, vcf2, query, args.gt_missing0, args.pl_missing0)
        try:
            for rec in recs:
                fh.write(VCFREC % rec)
        except BrokenPipeError:
            devnull = os.open(os.devnull, os.O_WRONLY)
            os.dup2(devnull, sys.stdout.fileno())
            sys.exit(1)
        fh.flush()
    fh.close()
    vcf1.close()
    vcf2.close()
