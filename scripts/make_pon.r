#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

if (!(length(args) %in% c(3,4)))
    stop("usage: make_pon.r input.tab metadata.csv output.tab [n.cores]")

inf <- args[1]
metaf <- args[2]
outf <- args[3]
n.cores <- 1
if (length(args) == 4)
    n.cores <- as.integer(args[4])
outf.gz <- paste0(outf, '.gz')


for (outfile in paste0(outf, c('', '.gz', '.gz.tbi'))) {
    if (file.exists(outfile))
        stop(sprintf("output file %s already exists, please delete it first", outf))
}

library(scan2)
library(Rsamtools)
library(future)
library(future.apply)
library(progressr)

if (n.cores > 1)
    plan(multicore, workers=n.cores)
else
    plan(sequential)


# These panels can be extremely large, so need to read in only the
# columns that really matter.
# Read the first 10 rows, just to get the header.
# Format is: chrom, pos, dbsnp, refnt, altnt, nalleles, N*(gt, ref, alt)
# where N* means repeat the 3 column pattern N times (N=sample number).
# Currently, this just skips the genotype and ref count columns for each
# sample.
read.alts.for.samples <- function(path, region, meta) {
    tf <- Rsamtools::TabixFile(path)
    open(tf)
    header <- read.tabix.header(tf)

    cols.to.read <- rep("NULL", tot.cols)
    cols.to.read[1:6] <- c('character', 'integer', 'character', 'character',
        'character', 'integer')
    df.samplenames <- c()
    for (i in seq(7, length(header), 3)) {
        if (header[i] %in% meta$sample) {
            cols.to.read[i+2] <- 'integer'
            sample.ids <- c(sample.ids, header[i])
        }
    }

    # Only read the alt counts for samples in the metadata table.
    tb <- read.tabix.data(path=path, region=region, header=header,
        colClasses=cols.to.read)

    # Since the GT column isn't read, sample names are lost
    colnames(tb)[-(1:6)] <- sample.ids
    tb
}
    
# Counts the number of unique cells and unique donors that support each mutation
# dmap - a named vector that maps sample ID -> donor ID
make.panel <- function(df, dmap, bulks) {
    bulk.idxs <- which(colnames(df) %in% bulks)
    cat("Bulks: ")
    print(bulks)
    cat("Bulk indexes: ")
    print(bulk.idxs)

    unique.cells <- rowSums(df[,-c(1:6, bulk.idxs)] > 0)
    unique.donors <- rowSums(do.call(cbind, lapply(unique(dmap), function(dn) {
        sns <- names(dmap)[dmap==dn]
        idxs <- which(colnames(df) %in% sns)
        rowSums(df[,idxs]) > 0
    })))
    unique.bulks <- rowSums(df[,bulk.idxs] > 0)
    # remove the entry with the maximum support
    outs <- apply(df[,-(1:6)], 1, function(row)
        row[-which.max(row)]
    )
    cat("Getting maximum out-group..\n")
    max.out <- apply(outs, 2, max)
    cat("Getting sum out-group..\n")
    sum.out <- apply(outs, 2, sum)
    cat("Getting sum bulk..\n")
    sum.bulk <- rowSums(df[,bulk.idxs])
    cat("Building final data.frame\n")
    data.table::data.table(df[,1:6], unique.donors=unique.donors,
        unique.cells=unique.cells,
        unique.bulks=unique.bulks,
        max.out=max.out, sum.out=sum.out,sum.bulk=sum.bulk,
        stringsAsFactors=FALSE)
}


# First load metadata. We will not read count data for any samples not
# in the metadata table.
cat(sprintf("Loading sample metadata %s..\n", metaf))
meta <- data.table::fread(metaf)

# GRanges intervals for chunked pipeline
grs <- tileGenome(seqlengths=seqinfo(genome.object)[as.character(1:22)],
                  tilewidth=10e6, cut.last.tile.in.chrom=TRUE)

# map donor <-> sample IDs
dmap <- meta$donor   
names(dmap) <- meta$sample

cat('Starting integrated table pipeline on', length(grs), 'chunks.\n')
cat('Parallelizing with', future::nbrOfWorkers(), 'cores.\n')

progressr::with_progress({
    p <- progressr::progressor(along=1:length(grs))
    xs <- future.apply::future_lapply(1:length(grs), function(i) {
        gr <- grs[i,]

        tb <- read.alts.for.samples(path=inf, region=gr, meta=meta)
        p(class='sticky', amount=0, paste0('read.alts.for.samples', i))

        ret <- make.panel(tb, dmap, meta$sample[meta$amp=='bulk'])
        p(class='sticky', amount=0, paste0('make.panel', i))

        p()
        ret
    })
})

pon <- rbindall(xs)

cat(sprintf("Saving PON to %s..\n", outf))
data.table::fwrite(pon, file=outf, sep='\t', quote=FALSE)
Rsamtools::bgzip(outf, outf.gz)
Rsamtools::indexTabix(file=outf.gz, format='vcf', comment='#')
