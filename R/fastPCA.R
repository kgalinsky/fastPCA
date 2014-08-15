read.geno <- function (filename) {
  X <- sapply( strsplit( readLines(filename), '' ), as.integer, simplify=T )
  X[X==9] <- NA
  return(t(X))
}

scale.geno <- function (X) {
   <- rowMeans(X, na.rm = TRUE)
  P <-  / 2
  S <- sqrt(2*P*(1-P))
  return((X-)/S)
}

rnorm.matrix <- function (m, n, ...) {
  matrix(rnorm(m*n, ...), m, n)  
}

rG <- function (A, l) { rnorm.matrix(ncol(A), l) }

# NORAL ALGORITHM

STEP1 <- function (A, l, i) {             # generate H
  G <- rG(A, l)
  H <- matrix(nrow=nrow(A), ncol=(i+1)*l)
  
  H[,1:l] <- A %*% G
  if(i > 0)
    for (j in 1:i)
      H[,j*l+1:l] <- A %*% crossprod(A, H[,(j-1)*l+1:l])
  
  return(H)
}

STEP2 <- function (H) qr.Q(qr(H))         # generate Q
STEP2.svd <- function (H) svd(H)$u

STEP3 <- function (A, Q) crossprod(A, Q)  # generate T

STEP45 <- function (T, Q, k) {
  # STEP 4
  T.svd <- svd(t(T), k, k)
  
  # STEP 5
  T.svd$u <- Q %*% T.svd$u
  T.svd$d <- T.svd$d[1:k]

  return(T.svd)
}

# ODIFICATIONS

STEP1.qr <- function (A, l, i) {             # generate H
  G <- rG(A, l)
  H <- matrix(nrow=nrow(A), ncol=(i+1)*l)
  
  H[,1:l] <- A %*% G
  if(i > 0)
    for (j in 1:i)
      H[,j*l+1:l] <- A %*% qr.Q(qr(crossprod(A, H[,(j-1)*l+1:l])))
  
  return(H)
}

# ALL TOGETHER

halko2011 <- function (A, k, l, i) {
  if (missing(l)) { l = k + 2 }
  if (missing(i)) { i = 2 }
  
  H <- STEP1(A, l, i)
  Q <- STEP2(H)
  T <- STEP3(A, Q)
  T.svd <- STEP45(T, Q, k)
  
  return(T.svd)
}

halko2011.qr <- function (A, k, l, i) {
  if (missing(l)) { l = k + 2 }
  if (missing(i)) { i = 2 }
  
  H <- STEP1.qr(A, l, i)
  Q <- STEP2(H)
  T <- STEP3(A, Q)
  T.svd <- STEP45(T, Q, k)
  
  return(T.svd)
}

halko2011.qr.svd <- function (A, k, l, i) {
  if (missing(l)) { l = k + 2 }
  if (missing(i)) { i = 2 }
  
  H <- STEP1.qr(A, l, i)
  Q <- STEP2.svd(H)
  T <- STEP3(A, Q)
  T.svd <- STEP45(T, Q, k)
  
  return(T.svd)
}


# TESTS

Fo <- function(A.svd, B.svd) {
  k <- min(ncol(A.svd$v), ncol(B.svd$v))
  f <- list(
    u = rowSums(crossprod(A.svd$u[,1:k], B.svd$u[,1:k])^2),
    v = rowSums(crossprod(A.svd$v[,1:k], B.svd$v[,1:k])^2),
    d = cor(A.svd$d[1:k], B.svd$d[1:k])^2
    )
  f['up'] = prod(f$u)
  f['vp'] = prod(f$v)
  return(f)
}
