#!/usr/bin/env Rscript

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)

source(paste(script.basename, 'lib.R', sep='/'))
  
args <- commandArgs(trailingOnly = TRUE)

X <- read.geno(args[1])
Y <- scale.geno(X)
Y.svd <- svd(Y)

e <- Y.svd$d**2
e <- e * length(e) / sum(e)
e[1] <- paste('#', e[1], sep='')

write.table(Y.svd$v[,1:10],
            sep="\t",
            row.names = F, col.names = e[1:10], quote=F)
