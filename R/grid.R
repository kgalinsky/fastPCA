#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]

k <- as.numeric(args[2])
L <- as.numeric(args[3])
I <- as.numeric(args[4])

source('lib.R')
A <- read.geno(filename)
A.svd <- svd(A)

FOM.std <- matrix(nrow=I+1, ncol=L-k)
FOM.mod <- matrix(nrow=I+1, ncol=L-2*k+1)
FOM.blc <- matrix(nrow=I+1, ncol=L-k)

G <- rG(A, L)
for (l in seq(k+1,L)) {
  write(l, file=stderr())
  R <- calc.R.blanczos(G[1:l,], A, I)

  FOM.std[,l-k] <- sapply(
    lapply( R, function(r) { calc.T.svd(r, A, k) } ),
    FOM.KG, A.svd=A.svd )

  if (l >= 2*k) {
    FOM.mod[,l-2*k+1] <- sapply(
      lapply( R, function(r) { calc.T.svd.mod(r, A, k) } ),
      FOM.KG, A.svd=A.svd )
  }
  
  FOM.blc[,l-k] <- sapply(
    lapply( seq(1, I+1), function(i) { calc.T.svd.mod(do.call(rbind, R[1:i]), A, k) } ),
    FOM.KG, A.svd=A.svd )
}

print("G")
write.table(G, row.names=F, col.names=F, sep="\t")

print("FOM.std")
write.table(FOM.std, row.names=F, col.names=F, sep="\t")

print("FOM.mod")
write.table(FOM.mod, row.names=F, col.names=F, sep="\t")

print("FOM.blc")
write.table(FOM.blc, row.names=F, col.names=F, sep="\t")
