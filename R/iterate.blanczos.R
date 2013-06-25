#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

i <- as.numeric(args[2])
k <- as.numeric(args[3])
l <- as.numeric(args[4])
r <- as.numeric(args[5])

source('lib.R')
A <- read.geno(filename)
A.svd <- svd(A)

write.table(
  replicate(r, { SVDs2FOMs(A.svd, fastPCA.blanczos.iterate(A, i, k, l))}),
            row.names=F, col.names=F, sep="\t")
