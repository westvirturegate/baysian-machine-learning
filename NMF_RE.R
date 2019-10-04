library(tcltk)
library(tensor)
library(jpeg)
img <-readJPEG("lady.jpeg")[,,1] #画像ファイルの座標は左上が原点なので
img_upright <- function(x){    #列にrevを作用させて正立させる
  t(apply(x,2,rev))
}


NMF_RE <- function(X,M,seed=1234){
  set.seed(seed)
  a_W=1
  b_W=1
  a_H=1
  b_H=1
  D <-nrow(X)
  N <-ncol(X)
  W<-matrix(rgamma(D*M,a_W,b_W),D,M)
  H<-matrix(rgamma(M*N,a_H,b_H),M,N)
  
  logW <-log(W)
  logH <-log(H)
  pb <- txtProgressBar(min = 1, max = 1000, style = 3)
  
  for(i in 1:1000){
    pi=array(1,dim=c(D,M,N))
    for (m in 1:M) pi[,m,]=exp(outer(logW[,m], logH[m,], "+"))
    den=apply(pi, c(1,3), sum)
    pihat=sweep(pi, c(1,3), den, "/")
    S = sweep(pihat, c(1,3), X, "*")
    
    ahat_W =apply(S, c(1,2), sum)+a_W
    bhat_W =apply(H, 1, sum)+b_W
    ahat_H =apply(S, c(2,3), sum)+a_H
    bhat_H =apply(W, 2, sum)+b_H
    
    W = sweep(ahat_W, 2, bhat_W,"/")
    H = sweep(ahat_H, 1, bhat_H, "/")
    
    logW = sweep(digamma(ahat_W), 2,log(bhat_W))
    logH = sweep(digamma(ahat_H), 1, log(bhat_H))
    
    
    
    setTxtProgressBar(pb, i)
  }
  return(list(W=W,H=H))
}


outputNMF_RE <-function(img,n){
  Mn<-NMF_RE(round(img_upright(img)*1000),M=n)
  image((Mn$W%*%Mn$H),
        main=paste('M=',n), col = grey(0:11/12))
  return(Mn=Mn)
}

Mn=outputNMF_RE(img, 10)
