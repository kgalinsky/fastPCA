#!/usr/bin/env Rscript

args <- as.numeric(commandArgs(trailingOnly = TRUE))
print(args)

m <- args[1]
n <- args[2]
k <- args[3]
l <- args[4]
i.max <- args[5]
rep <- args[6]

gauss_matrix <- function (m, n) {
  matrix(rnorm(m*n), m, n)
}

fastPCA.iterate <- function (A, k=10, i.max=50, l=20, G) {
  if (missing(G)) {
    G <- gauss_matrix(l, dim(A)[1])
  }
  
  R <- G %*% A
  R <- R / mean(R)

  T.svds <- list();
  
  for (i in 1:i.max) {
    R <- tcrossprod(R, A)
    R <- R %*% A
    R <- R / sqrt(mean(R^2))

    Q <- svd(R, nv=k)$v
    T <- A %*% Q

    T.svd <- svd(T)
    T.svd$v <- Q %*% T.svd$v
    
    T.svds[[i]] <- T.svd
  }
  
  return(T.svds)
}

FOM.KG <- function (A.svd, T.svd) {
  k <- length(T.svd$d)
  prod(sqrt(rowSums(crossprod(A.svd$u[,1:k], T.svd$u)^2)))
}

SVDs2FOMs <- function (A.svd, T.svds) {
  sapply(T.svds, FOM.KG, A.svd=A.svd)
}

A <- gauss_matrix(m, n)
A.svd <- svd(A)

A.T.svds <- replicate(rep, fastPCA.iterate(A, k, i.max, l), simplify=F)
A.T.FOMs <- lapply(A.T.svds, SVDs2FOMs, A.svd=A.svd)

write.table(do.call(cbind, A.T.FOMs))
