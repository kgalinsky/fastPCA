read.pop <- function (filename) {
  apply(
    sapply(
      strsplit(
        scan(filename, what='character'),
        '' ),
      as.numeric, simplify=T )
    , 2, function (x) {
      m <- mean(x)
      p <- m / 2
      ( x - m ) / 2 / p / (1-p)
    })
}

gauss_matrix <- function (m, n) {
  matrix(rnorm(m*n), m, n)
}

fastPCA.iterate <- function (A, l=20, i.max=50, k, G) {
  if (missing(G)) {
    G <- gauss_matrix(l, dim(A)[1])
  }
  if (missing(k)) {
    k <- l - 1
  }
  
  R <- G %*% A
  R <- R / mean(R)
  
  T.svds <- list();
  
  for (i in 1:i.max) {
    R <- tcrossprod(R, A)
    R <- R %*% A
    R <- R / mean(R)
    
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

FOM.NP <- function (A.svd, T.svd) {
  k <- length(T.svd$d)
  diff <- A.svd$d[1:k]^2 - T.svd$d^2
  sqrt(sum(diff^2))
}

SVDs2FOMs <- function (A.svd, T.svds, FOM=FOM.KG) {
  sapply(T.svds, FOM, A.svd=A.svd)
}