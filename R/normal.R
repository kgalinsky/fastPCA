library('RColorBrewer')
source('lib.R')

plot.traces <- function (X, X.svd, K, Ls, I, R=5) {
  if (missing(X.svd)) { X.svd <- svd(X) }

  pal=brewer.pal(n=8, 'Dark2')
  
  plot(c(),
       xlim=c(0, I), xlab='I',
       ylim=c(0, 1), ylab='FO'
  )
  abline(h=.99, lty=2)
  abline(h=.95, lty=3)
  
  leg <- sapply(Ls, function (L) { sprintf("L=%d", L) })
  legend('bottomright', legend=leg, lt=1, col=pal)
  
  Is <- seq(0, I)
  
  for (l in 1:length(Ls)) {
    L <- Ls[l]
    for (r in 1:R) {
      Ts.svd <- fastPCA.blanczos.iterate(X, i=I, l=L, k=K)
      Ts.FO <- sapply(Ts.svd, function (T.svd) { FOM.KG(X.svd, T.svd) })
      lines(is, Ts.FO, col=pal[l])      
    }
  }
}

lines.trace <- function(A, col='red', ...) {
  Ts.svd <- fastPCA.blanczos.iterate(X, ...)
  Ts.FO <- sapply(Ts.svd, function (T.svd) { FOM.KG(X.svd, T.svd) })
  lines(is, Ts.FO, col=col) 
}

#X <- rpop(, N)
X.10k.1k <- rnorm.matrix(1000, 10000)
X.10k.1k.svd <- svd(X)


plot.traces(X.10k.1k, X.10k.1k.svd, Ls=c(2,3,5,10), K=1, I=10, R=3)

X.5k.1k <- rnorm.matrix(1000, 5000)
X.5k.1k.svd <- svd(X.5k.1k)

plot.traces(X.5k.1k, X.5k.1k.svd, Ls=c(2,3,5,10), K=1, I=10, R=3)