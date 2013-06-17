#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

k <- args[2]
l <- args[3]
i.max <- args[4]
rep <- args[5]

source('lib.R')
A <- read.pop(args[1])
A.svd <- svd(A)

T.svds <- fastPCA.iterate(A, l, i.max, k)
SVDs2FOMs(A.svd, T.svds)
