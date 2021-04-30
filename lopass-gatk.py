#!/usr/bin/env python3
import os
import sys
import math
import argparse
import subprocess
import multiprocessing
from collections import namedtuple

DEFAULT_INTERVALS = (
    ## chrom start end xx xy
    ("chr1", "1", "248956422", "2", "2"),
    ("chr2", "1", "242193529", "2", "2"),
    ("chr3", "1", "198295559", "2", "2"),
    ("chr4", "1", "190214555", "2", "2"),
    ("chr5", "1", "181538259", "2", "2"),
    ("chr6", "1", "170805979", "2", "2"),
    ("chr7", "1", "159345973", "2", "2"),
    ("chr8", "1", "145138636", "2", "2"),
    ("chr9", "1", "138394717", "2", "2"),
    ("chr10", "1", "133797422", "2", "2"),
    ("chr11", "1", "135086622", "2", "2"),
    ("chr12", "1", "133275309", "2", "2"),
    ("chr13", "1", "114364328", "2", "2"),
    ("chr14", "1", "107043718", "2", "2"),
    ("chr15", "1", "101991189", "2", "2"),
    ("chr16", "1", "90338345", "2", "2"),
    ("chr17", "1", "83257441", "2", "2"),
    ("chr18", "1", "80373285", "2", "2"),
    ("chr19", "1", "58617616", "2", "2"),
    ("chr20", "1", "64444167", "2", "2"),
    ("chr21", "1", "46709983", "2", "2"),
    ("chr22", "1", "50818468", "2", "2"),
    ("chrX", "1", "2781479", "2", "1"),
    ("chrX", "2781480", "155701382", "2", "2"),
    ("chrX", "155701383", "156030895", "2", "1"),
    ("chrY", "1", "57227415", "0", "1")
)

DEFAULT_CHROMS = ["chr%s" % c for c in range(1,23)] + ["chrX"]

GVCF_CMD = [
    "gatk", "HaplotypeCaller",
    "--emit-ref-confidence", "GVCF",
    "--min-dangling-branch-length", "1",
    "--min-pruning", "1",
    "--indel-size-to-eliminate-in-ref-model", "0",
    "--heterozygosity", "0.1",
    "--indel-heterozygosity", "0.1",
    "--base-quality-score-threshold", "10",
    "--dont-use-soft-clipped-bases", "true",
    "--create-output-variant-index", "false"
]

GATHER_CMD = [
    "gatk", "GatherVcfs",
    "--CREATE_INDEX", "false"
]

def run(cmd):
    out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    return(out)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ## I/O
    parser.add_argument("-i", '--input', help="Input BAM/CRAM file", required=True)
    parser.add_argument("-o", '--output', help="output file prefix", required=True)
    parser.add_argument("-s", '--sex', help="Genetic sex XX or XY", required=True)
    parser.add_argument("-y", '--ycall', help="Call on chrY", action="store_true")
    parser.add_argument("-c", '--chroms', help="Optionally restrict intervals to chromosomes")
    parser.add_argument("-d", '--delete', help="delete temporary files", action="store_true")
    ## GATK options
    parser.add_argument("-R", '--reference', help="Reference FASTA file", required=True)
    parser.add_argument("-L", '--intervals', help="File with call intervals")
    parser.add_argument("-J", '--javaopt', help="Java options")
    ## Parallel options
    parser.add_argument("-p", '--processes', help="Number of parallel GATK processes", default=4)
    parser.add_argument("-t", '--threads', help="Number of parallel GATK pair-hmm threads", default='1')
    args=parser.parse_args()

    #### prepare calling intervals
    intervals = []
    if args.intervals:
        if os.path.exists(args.intervals) and os.path.isfile(args.intervals):
            with open(args.intervals, "r") as fh:
                for line in fh:
                    interval = tuple(line.strip().split("\t"))
                    intervals.append(interval)                    
        else:
            lines = args.intervals.split(",")
            for line in lines:
                chr, begend, xx, xy = line.split(":")
                beg, end = begend.split("-")
                intervals.append((chr, beg, end, xx, xy))
    else:
        intervals = DEFAULT_INTERVALS

    ## prepare chromosomes
    if args.chroms:
        chroms = args.chroms.split(",")
    else:
        chroms = DEFAULT_CHROMS
        if args.ycall:
            chroms.append("chrY")

    ## prepare commands
    cmds = []
    vcfs = []
    for (i, (chr, beg, end, xx, xy)) in enumerate(intervals):
        ## call only selected chromosomes
        if not chr in chroms:
            continue
        ## I/O
        index = "%0.4d" % i
        interval = "%s:%s-%s" % (chr, beg, end)
        vcf = "-".join((args.output, index)) + ".vcf"
        ## ploidy
        if args.sex.lower() == "xy":
            ploidy = xy
        else:
            ploidy = xx

        ## GATK4 command
        cmd = GVCF_CMD + [
            "-R", args.reference,
            "-L", interval,
            "-I", args.input,
            "-O", vcf,
            "--sample-ploidy", ploidy,
            "--native-pair-hmm-threads", args.threads,
            "--verbosity", "ERROR"
        ]
        cmds.append(cmd)
        vcfs.append(vcf)

    ## Run GATK4 in parallel
    with multiprocessing.Pool(processes=args.processes) as pool:
        for _ in pool.imap_unordered(run, cmds, chunksize=1):
            pass

    ## GatherVcfs
    outvcf = args.output + ".vcf"
    cmd = GATHER_CMD + [
        "--OUTPUT", outvcf,
        "--VERBOSITY", "ERROR"
    ]
    for vcf in vcfs:
        cmd.extend(("--INPUT", vcf))
    subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

    ## compress and index
    cmd = ["bgzip", outvcf]
    subprocess.run(cmd, subprocess.DEVNULL)
    cmd = ["tabix", outvcf + ".gz"]
    subprocess.run(cmd, subprocess.DEVNULL)

    ## clean-up
    if args.delete:
        for vcf in vcfs:
            os.remove(vcf)
