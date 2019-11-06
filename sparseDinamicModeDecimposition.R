library(cmvnorm)
library(invgamma)
library(jpeg)
img <-readJPEG("lady.jpg")[,,1] #画像ファイルの座標は左上が原点なので
image_upright <- function(x){    #列にrevを作用させて正立させる
  image(t(apply(x,2,rev)),col = grey(0:11/12))
}

image_upright(img)

X = dat2[1:49,]
Y = dat2[2:50,]

gibbsiter <- 1000
burnin <- 20
interval <- 4

mcemiter <- 5



N = dim(X)[1]
D = dim(X)[2]
K = 5





rgig <- function(P, a, b, sampleSize) {
  
  lambda = P
  omega = sqrt(a*b)
  alpha = sqrt(omega^2 + lambda^2) - lambda
  
  x = as.numeric(-psi(1, alpha, lambda))
  
  if((x >= 0.5) && (x <= 2)){
    t = 1
  }else if(x > 2){
    t = sqrt(2 / (alpha + lambda))
  }else if(x < 0.5){
    t = log(4/(alpha + 2*lambda))
  }
  
  x = as.numeric(-psi(1, alpha, lambda))
  if((x >= 0.5) && (x <= 2)){
    s = 1
  }else if(x>2){
    s = sqrt(4/(alpha*cosh(1) + lambda))
  }else if(x<0.5){
    s = min(1/lambda, as.numeric(log(1 + 1/alpha + sqrt(1/alpha^2+2/alpha))))
  }
  
  eta = -psi(t, alpha, lambda)
  zeta = -dpsi(t, alpha, lambda)
  theta = -psi(-s, alpha, lambda)
  xi = dpsi(-s, alpha, lambda)
  p = 1/xi
  r = 1/zeta
  td = t - r*eta
  sd = s - p*theta
  q = td + sd
  
  
  
  X = matrix(0, sampleSize, 1)
  for (sample in 1:sampleSize) {
    done = FALSE
    while (!done) {
      U = runif(1)
      V = runif(1)
      W = runif(1)
      
      if(U < as.numeric(q / (p + q + r))){
        X[sample] = -sd + q*V
      }else if(U < as.numeric((q + r) / (p + q + r))){
        X[sample] = (td - r*log(V))
      }else {
        X[sample] = (-sd + p*log(V))
      }
      f1 = exp(-eta - zeta*(X[sample]-t))
      f2 = exp(-theta + xi*(X[sample]+s))
      if(as.numeric(W*g(X[sample], sd, td, f1, f2)) <=
         as.numeric(exp(psi(X[sample], alpha, lambda)))){
        done = TRUE
      }
    }
  }
  X = exp(X) * (lambda / omega + sqrt(1 + (lambda/omega)^2))
  X = X / sqrt(a/b)
  
  # cat("p =",p,"\n")
  # cat("q =",q,"\n")
  # cat("r =",r,"\n")
  # cat("U =",U,"\n")
  # browser()
  return(X)
}


psi <- function(x, alpha, lambda) {
  -alpha*(cosh(x) - 1) - lambda*(exp(x) - x - 1)
}

dpsi <- function(x, alpha, lambda) {
  -alpha*sinh(x) - lambda*(exp(x) - 1)
}

g <- function(x, sd, td, f1, f2) {
  
  x = as.numeric(x)
  a = 0
  b = 0
  c = 0
  
  if((x >= as.numeric(-sd)) && (x <= as.numeric(td))){
    a = 1
  }else if(x > as.numeric(td)){
    b = f1
  }else if(x < as.numeric(-sd)){
    c = f2
  }
  return(a + b + c)
}


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


get_Xi <- function(W, Phi, k, tau) {
  K = nrow(W)
  Xi = Phi[,-k,drop=FALSE]%*%diag(tau[-k])%*%W[-k,,drop=FALSE]
  return(Xi)
}

get_Eta <- function(Lambda, W, Phi, k, tau) {
  if(k < 2)
    K = nrow(W)
  Eta = Phi[,-k,drop=FALSE]%*%diag(c(Lambda[,-k]))%*%diag(tau[-k])%*%W[-k,,drop=FALSE]
  return(Eta)
}



##hyperparameters
alpha = 1e-3
beta = 1e-3

##initialization K=<D only
fit = sdmd(t(Conj(X)), t(Conj(Y)), K)
Lambda = t(fit$evalue)
W = t(fit$mode)
Phi = t(fit$efun)

Lambda[is.nan(Lambda)] = 1.0
W[is.nan(W)] = 1.0
Phi[is.nan(Phi)] = 1.0

A = rep(1,K)
F1 = X - Phi%*%diag(A)%*%W
F2 = Y - Phi%*%diag(c(Lambda))%*%diag(A)%*%W
Sigma2 = sum(sum(F1*Conj(F1) + F2*Conj(F2)))/(2*N*D-1)
Nu2 = rep(1,K)
V2 = matrix(100,K,D)
Tau2 = rep(1,K)


Gamma_v = t(Conj(D*sqrt(Sigma2)/rowSums(abs(W),2)))

##memory allocation
num_sample  <- floor((gibbsiter - burnin)/interval)
sample.Lambda <- matrix(0, gibbsiter, K)
sample.W <- array(0, c(gibbsiter, K, D))
sample.Phi <- array(0, c(gibbsiter, N, K))
sample.Sigma2 <- matrix(0, gibbsiter, 1)
sample.Nu2 <- matrix(0, gibbsiter, K)
sample.V2 <- array(0, c(gibbsiter, K, D))


##main iteration of MC-EM
# for (mcem in 1:mcemiter) {

count_sample <- 1

for (gibbs in 1:gibbsiter) {
  
  ###update of Lambda
  
  Lambda_new <- Lambda
  
  for (k in 1:K) {
    Eta <- get_Eta(Lambda, W, Phi, k, A)
    p_l = A[k]*A[k]*Conj(W[k,])%*%(W[k,]) * 
      (t(Conj(Phi[,k]))%*%Phi[,k])/Sigma2 + 1/Nu2[k]
    m_l = A[k] *Conj(W[k,])%*%
      colSums(sweep(Y - Eta, MARGIN = 1,Conj(Phi[,k]), "*"))/
      p_l/Sigma2
    
    Lambda_new[k] = rnorm(1,as.numeric(m_l), as.numeric(1/p_l))
    # Lambda_new[k] = cnormrnd(m_l, 1/p_l)
  }
  
  Lambda = Lambda_new
  
  ###update of W
  W_new = W
  
  for (k in 1:K) {
    
    tmp = V2[k,]
    
    Xi = get_Xi(W, Phi, k, A)
    Eta <- get_Eta(Lambda, W, Phi, k, A)
    P_W = A[k]*A[k]* (1+Conj(Lambda[k])*Lambda[k])/Sigma2 *
      as.complex((t(Phi[,k])%*%Phi[,k]))*diag(rep(1,D))  + 
      diag(1.0/tmp)/Sigma2  ########
    P_inv_W = diag(1/diag(P_W))
    P_inv_W = 0.5*(P_inv_W + t(Conj(P_inv_W))) ####
    m_W =A[k]* ((colSums(sweep(X-Xi, MARGIN = 1, Conj(Phi[,k]),"*"))+
                   colSums(sweep(Y-Eta, MARGIN = 1, Conj(Phi[,k]),"*")))/
                  Sigma2)%*%Conj(t(P_inv_W))
    
    W_new[k,] = rmvnorm(1, mean = c(m_W), sigma = P_inv_W)
    # W_new[k,] = cnormrnd(c(m_W),  P_inv_W)
  }
  
  W = W_new
  
  ###update of Phi
  P = diag(A)%*%Conj(W)%*%t(W)%*%diag(A)/Sigma2 +
    diag(c(Conj(Lambda))) %*% diag(A) %*% Conj(W) %*% t(W)%*%
    diag(A) %*% diag(c(Lambda))/Sigma2 +
    diag(rep(1,K))
  
  P = 0.5*(P + t(Conj(P)))
  UP = eigen(P)$vectors
  DP = eigen(P)$values
  
  P_inv = UP%*%diag(1/DP)%*%t(Conj(UP))
  m = ((X %*% t(Conj(W)) %*% diag(A) + 
          Y %*% t(Conj(W))%*%diag(A)%*%
          diag(Conj(c(Lambda))))/Sigma2) %*%t(P_inv)
  for (n in 1:N) {
    Phi[n,] = rmvnorm(1, mean = c(m[n,]), sigma = P_inv) 
    # Phi[n,] = cnormrnd(c(m[n,]), P_inv) 
  }
  
  
  ###update of Sigma2
  F1 = X - Phi%*%diag(A)%*%W
  F2 = Y - Phi%*%diag(c(Lambda))%*%diag(A)%*%W
  a_sigma2 = 2*N*D + alpha + 0.5*K
  b_sigma2 = sum(sum(Conj(F1)*F1,2),1) + sum(sum(Conj(F2)*F2,2),1) + 
    beta +
    0.5*sum(sum(Conj(W)*W/V2))
  
  Sigma2 = invgamma::rinvgamma(1,a_sigma2,as.numeric(1/b_sigma2))
  
  
  
  ###update of V2
  a_V2 = Gamma_v^2
  b_V2 = matrix(as.numeric(Conj(W)*W/Sigma2), K , D) ###as.numeric
  for (k in 1:K) {
    for (d in 1:D) {
      V2[k,d] = rgig(0.5, a_V2[k], b_V2[k,d], 1)
    }
  }
  
  
  
  sample.Lambda[gibbs,] = Lambda
  sample.W[gibbs,,] = W
  sample.Phi[gibbs,,] = Phi
  sample.Sigma2[gibbs,] = Sigma2
  sample.Nu2[gibbs,] = Nu2
  sample.V2[gibbs,,]  = V2
  
}


##########################









plot(Re(Phi[,1]), type = "l")
plot(Re(Phi[,2]), type = "l")
plot(Re(Phi[,3]), type = "l")
plot(Re(Phi[,4]), type = "l")
plot(Re(Phi[,5]), type = "l")
plot(Re(Phi[,6]), type = "l")
plot(Re(Phi[,]), type = "l")

matplot(Re(W))


matplot(Phi%*%W,type = "l")
matplot(dat2,type = "l")
