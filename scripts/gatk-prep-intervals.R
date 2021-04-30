library(data.table)
library(stringr)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)

tmp <- fread("data/wgs_calling_regions.hg38.interval_list", skip = 3368)
rng.call <- with(tmp, GRanges(V1, IRanges(V2, V3), V4))

gseqi <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
gseqi <- keepStandardChromosomes(gseqi)
gseqi <- dropSeqlevels(gseqi, c("chrM"))
seqinfo(rng.call) <- gseqi
rng.call$region <- "target"

rng.gaps <- gaps(rng.call)
rng.gaps <- rng.gaps[strand(rng.gaps)=="+"]
start(rng.gaps) <- start(rng.gaps) - 1
rng.gaps <- trim(rng.gaps)
rng.gaps$region <- "gap"

rng <- reduce(sort(c(rng.call, rng.gaps)), min.gapwidth=0)
rng$xx.ploidy <- 2
rng$xy.ploidy <- 2
rng[seqnames(rng)=="chrY"]$xx.ploidy <- 0
rng[seqnames(rng)=="chrY"]$xy.ploidy <- 1

## PARs
par1.end <- GRanges("chrX", IRanges(2781479, 2781479), "+")
par2.beg <- GRanges("chrX", IRanges(155701383, 155701383), "+")

par1.split <- which(rng %over% par1.end)
end(rng[par1.split-1]) <- 2781479
start(rng[par1.split]) <- 2781479+1
rng[1:(par1.split-1)]$xy.ploidy <- 1

par2.split <- which(rng %over% par2.beg)
end(rng[par2.split-1]) <- 155701383-1
start(rng[par2.split]) <- 155701383
rng[par2.split:last(which(seqnames(rng)=="chrX"))]$xy.ploidy <- 1

rng.tbl <- as.data.table(rng)[,.(chr=seqnames, start, end, xx.ploidy, xy.ploidy)]
fwrite(rng.tbl, "data/wgs_calling_full_regions.hg38.intervals", sep="\t", col.names=FALSE)
