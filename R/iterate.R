#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

k <- args[2]
l <- args[3]
i <- args[4]
r <- args[5]

source('lib.R')
A <- read.pop(filename)
A.svd <- svd(A)

A.T.svds <- fastPCA.iterate(A, i, k, l)
write.table(
  replicate(r, { SVDs2FOMs(A.svd, fastPCA.iterate(A, i, k, l) }), 
            row.names=F, col.names=F, sep="\t")