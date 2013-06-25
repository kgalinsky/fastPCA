read.geno <- function (filename, ...) {
  sapply(
    strsplit(
      scan(filename, what='character', ...),
      '' ),
    as.numeric, simplify=T )
}

scale.geno <- function (x) {
  apply( x, 2, scale.geno.col )
}

scale.geno.col <- function (x) {
  m <- mean(x)
  p = m/2
  (x-m) / sqrt(2*p*(1-p))
}

rnorm.matrix <- function (m, n, ...) {
  matrix(rnorm(m*n, ...), m, n)  
}

rG <- function (A, l) { rnorm.matrix(l, dim(A)[1]) }

rgeno <- function(n, p) {
  rbinom(n, 2, p)
}

calc.R <- function (G, A, i) {
  R <- G %*% A
  for (j in 1:i) {
    R <- calc.R.it(R, A)
  }
  return(R)
}

calc.R.blanczos <- function (G, A, i) {
  R <- list()
  R[[1]] <- G %*% A
  for (j in 1:i) {
    R[[j+1]] <- calc.R.it(R[[j]], A)
  }
  return(R)
}

calc.R.it <- function (R, A) {
  return( tcrossprod(R, A) %*% A * prod(dim(R)) / sqrt(sum(R^2)) )
}

calc.T.svd <- function (R, A, k) {
  Q <- svd(R, nv=k)$v
  T <- A %*% Q
  
  T.svd <- svd(T)
  T.svd$v <- Q %*% T.svd$v
  
  return(T.svd)
}

calc.T.svd.mod <- function (R, A, k) {
  Q <- svd(R)$v
  T <- A %*% Q
  
  Tp.svd <- svd(T)
  
  T.svd <- list(u=Tp.svd$u[,1:k], d=Tp.svd$d[1:k])
  T.svd$v <- (Q %*% Tp.svd$v)[,1:k]
  
  return (T.svd)
}

calc.T.svd.smart <- function (R, A, k) {
  if (2*k <= dim(R)[1]) {
    return(calc.T.svd.mod(R, A, k))
  }
  else {
    return(calc.T.svd(R, A, k))
  }
}

fastPCA <- function (A, i=50, k=10, l=20, G) {
  if (missing(G)) { G <- rG(A, l) }
  return(calc.T.svd(calc.R(G, A, i), A, k))
}

fastPCA.mod <- function (A, i=50, k=10, l=20, G) {
  if (missing(G)) { G <- rG(A, l) }
  return(calc.T.svd.mod(calc.R(G, A, i), A, k))
}

fastPCA.iterate <- function (A, i=50, k=10, l=20, G) {
  if (missing(G)) { G <- rG(A, l) }
  return( lapply( calc.R.blanczos(G, A, i)[1:i+1],
                  function(r) { calc.T.svd(r, A, k) } ) )
}

fastPCA.mod.iterate <- function (A, i=50, k=10, l=20, G) {
  if (missing(G)) { G <- rG(A, l) }
  return( lapply( calc.R.blanczos(G, A, i)[1:i+1],
                  function(r) { calc.T.svd.mod(r, A, k) } ) )
}

fastPCA.blanczos <- function (A, i=50, k=10, l=20, G) {
  if (missing(G)) { G <- rG(A, l) }
  return( calc.T.svd.mod( do.call(rbind, calc.R.blanczos(G, A, i)), A, k ) )
}

fastPCA.blanczos.iterate <- function (A, i=50, k=10, l=20, G) {
  if (missing(G)) { G <- rG(A, l) }
  R <- calc.R.blanczos(G, A, i)
  return( lapply( seq(1, i+1),
                  function(j) { calc.T.svd.mod(do.call(rbind, R[1:j]), A, k) } ) )
}

FOM.KG <- function (A.svd, T.svd) {
  k <- length(T.svd$d)
  prod(sqrt(rowSums(crossprod(A.svd$u[,1:k], T.svd$u)^2)))
}

FOM.NP <- function (A.svd, T.svd) {
  k <- length(T.svd$d)
  diff <- A.svd$d[1:k]^2 - T.svd$d^2
  sqrt(sum(diff^2))
}

SVDs2FOMs <- function (A.svd, T.svds, FOM=FOM.KG) {
  sapply(T.svds, FOM, A.svd=A.svd)
}
