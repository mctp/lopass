#!/usr/bin/env python3
import os
import sys
import tempfile
import argparse
import subprocess
import functools
import pathlib
import multiprocessing

ALL_CHR = ["chr%s" % i for i in range(1, 22+1)] + ["chrX_par1", "chrX", "chrX_par2"]

def write_ploidy(fh, args):
    out = subprocess.run(["bcftools", "query", "-l", args.vcf], stdout=subprocess.PIPE)
    samples = out.stdout.decode("ascii").strip().split("\n")
    out = subprocess.run(["bcftools", "view", "-H", args.vcf, args.xsnp], stdout=subprocess.PIPE)
    formats = out.stdout.decode("ascii").strip().split("\t")[9:]
    gts = [fmt.split(":")[0] for fmt in formats]
    ploids = [2 if gt in ("0/0", "0/1", "1/0", "1/1", "./.") else 1 for gt in gts]
    for sample, ploid in zip(samples, ploids):
        line = "%s\t%s\n" % (sample, ploid)
        fh.write(line)
        fh.flush()

def glimpse(chrom, maps, chunks, args):
    ## run GLIMPSE on each chromosome
    map_ = maps[chrom]
    chunks_ = chunks[chrom]
    with open(chunks_, "r") as fh:
        outs_ = []
        for chunk_ in fh:
            id_, chr_, irg_, org_, _, _ = chunk_.split("\t")
            out_ = args.output + "-%s-%0.2d.bcf" % (chr_, int(id_))
            ## run GLIMPSE phase
            cmd = ["GLIMPSE_phase",
                   "--input", args.vcf, "--output", out_,
                   "--reference", args.panel, "--map", map_,
                   "--input-region", irg_, "--output-region", org_,
                   "--burnin", str(args.burnin),
                   "--main", str(args.main),
                   "--thread", str(args.glimpse_thread),
                   ]
            if chrom == "chrX":
                with tempfile.NamedTemporaryFile(mode="w") as temp:
                    write_ploidy(temp, args)
                    cmd = cmd + ["--samples-file", temp.name]
                    subprocess.run(cmd, stdout=subprocess.DEVNULL)
            else:
                subprocess.run(cmd, stdout=subprocess.DEVNULL)
            # index bcf file
            cmd = ["bcftools", "index", "-f", out_]
            subprocess.run(cmd)
            outs_.append(out_)

        # ligate
        imp_ = args.output + "-" + chrom + "-" + "imputed.bcf"
        phs_ = args.output + "-" + chrom + "-" + "phased.bcf"

        with tempfile.NamedTemporaryFile(mode="w") as temp:
            temp.writelines("%s\n" % l for l in outs_)
            temp.flush()
            cmd = ["GLIMPSE_ligate",
                   "--input", temp.name,
                   "--output", imp_,
                   "--thread", str(args.glimpse_thread),
                   ]
            subprocess.run(cmd, stdout=subprocess.DEVNULL)

        cmd = ["GLIMPSE_sample",
               "--input", imp_,
               "--output", phs_,
               "--solve",
               "--thread", str(args.glimpse_thread),
               ]
        subprocess.run(cmd, stdout=subprocess.DEVNULL)

        if args.delete:
            for fn in outs_:
                os.remove(fn)
                os.remove(fn + ".csi")
                
        return (chrom, (imp_, phs_))

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    ## main
    parser.add_argument("panel", help="Reference panel")
    parser.add_argument("vcf", help="VCF file returned by lopass-genotype")
    ## required
    parser.add_argument("-o", '--output', help="Output prefix", required=True)
    parser.add_argument("-c", '--chunks', help="GLIMPSE chunks directory", required=True)
    parser.add_argument("-m", '--maps', help="GLIMPSE maps directory", required=True)
    ## wrapper settings
    parser.add_argument("-r", '--chromosomes', help="Select subset of chromosomes to impute")
    parser.add_argument("-x", '--xsnp', help="Panel SNP in chrX non-PAR region to infer sex", default="chrX:100000139")
    parser.add_argument("-t", '--thread', help="Number of threads", default=4)
    parser.add_argument("-d", '--delete', help="delete temporary files", action="store_true")
    ## GLIMPSE settings
    parser.add_argument('--glimpse_thread', help="Number of threads to pass to GLIMPSE", default=1)
    parser.add_argument('--burnin', help="GLIMPSE_phase burnin", default=15)
    parser.add_argument('--main', help="GLIMPSE_phase main", default=15)
    args=parser.parse_args()

    ## create output directory
    out_dir = os.path.dirname(args.output)
    pathlib.Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    ## chromosomes
    if args.chromosomes:
        chroms = args.chromosomes.split(",")
        chroms = list(set(chroms) & set(ALL_CHR))
    else:
        chroms=ALL_CHR

    ## maps and chunks
    map_fns = os.listdir(args.maps)
    maps = dict([(fn.split(".")[0], os.path.join(args.maps, fn)) for fn in map_fns if fn.endswith(".gmap.gz")])
    chunk_fns = os.listdir(args.chunks)
    chunks = dict([(fn.split(".")[0], os.path.join(args.chunks, fn)) for fn in chunk_fns if fn.endswith(".txt")])

    ## run GLIMPSE
    impd = {}
    phsd = {}
    with multiprocessing.Pool(int(args.thread)) as pool:
        chrom_lp = pool.map(functools.partial(glimpse, maps=maps, chunks=chunks, args=args), chroms, chunksize=1)
        for chrom, (l, p) in chrom_lp:
            impd[chrom] = l
            phsd[chrom] = p
    impl = [impd[chrom] for chrom in ALL_CHR if chrom in impd]
    phsl = [phsd[chrom] for chrom in ALL_CHR if chrom in phsd]
    cmd = ["bcftools", "concat", "-o", args.output + "-imputed.bcf", "-Ob"] + impl
    subprocess.run(cmd, subprocess.DEVNULL)
    cmd = ["bcftools", "index", "-f", args.output + "-imputed.bcf"]
    subprocess.run(cmd, subprocess.DEVNULL)
    cmd = ["bcftools", "concat", "-o", args.output + "-phased.bcf", "-Ob"] + phsl
    subprocess.run(cmd, subprocess.DEVNULL)
    cmd = ["bcftools", "index", "-f", args.output + "-phased.bcf"]
    subprocess.run(cmd, subprocess.DEVNULL)
    
    if args.delete:
        for fn in impl+phsl:
            os.remove(fn)
            os.remove(fn + ".csi")
