source('lib.R')

M <- 1000
N <- 100
K <- 1
L <- 10
I <- 10

#X <- rpop(M, N)
X <- rnorm.matrix(M, N)
X.svd <- svd(t(X))

Ls <- c( 2, 3, 5, 10 )
legend <- sapply(Ls, function (L) { sprintf("L = %d", L) })
cols <- c( 'red', 'green', 'blue', 'purple' )

plot( c(), xlim=c(0,10), ylim=c(0,1), xlab='i', ylab='FoM' )
legend('bottomright', legend=legend, lw=1, col=colors)

for (i in 1:4) {
  L <- Ls[i]
  col <- cols[i]
  for (j in 1:5) {
    G <- rnorm.matrix(N, L)
    
    H <- list()
    H[[1]] <- X %*% G
    for (i in 1:I) {
      H[[i+1]] <- X %*% crossprod(X, H[[i]])
    }
    
    H.svd <- lapply(seq(1, I+1), function (i) {svd(do.call(cbind, H[1:i]))})
    T     <- lapply(H.svd, function (h.svd) { crossprod(X, h.svd$u) })
    T.svd <- lapply(T, function (x) { svd(x, k, nu=K) })
    
    lines( 
      seq(0, 10),
      sapply(T.svd, function(x) { FOM.KG(X.svd, x) }),
      col=col,
    )
  }
}