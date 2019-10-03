sdmd <- function(X, Y, d_X) {
  svd <- svd(X)
  U <- svd$u
  S <- svd$d
  V <- svd$v
  
  diag_S = sort(S, decreasing = TRUE, index.return = TRUE)$x
  idx = sort(S, decreasing = TRUE, index.return = TRUE)$ix
  
  U = U[,idx[1:d_X]]
  S = diag(diag_S[1:d_X])
  V = V[,idx[1:d_X]]
  
  
  M <- Y%*%V%*%solve(S)
  A_til <- t(Conj(U))%*%M
  w <- eigen(A_til)$vectors
  D <- eigen(A_til)$values
  z <- eigen(t(A_til))$vectors
  
  #normalization
  N = Conj(t(Conj(z))%*%w)
  z = z%*%solve(N)
  
  D = sort(D, decreasing = TRUE, index.return = TRUE)$x
  idx = sort(D, decreasing = TRUE, index.return = TRUE)$ix
  phi = (M%*%w[,idx])/D
  kap = U%*%z[,idx]
  
  return(list(
    evalue = D,
    mode = Conj(kap),
    efun = t(Conj(phi))%*%X,
    levec = kap
  ))
}
