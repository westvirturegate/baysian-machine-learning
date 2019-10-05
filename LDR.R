library(png)
library(jpeg)
library(tcltk)
img <-readPNG("lady.png")[,,1] #画像ファイルの座標は左上が原点なので
img_upright <- function(x){    #列にrevを作用させて正立させる
  t(apply(x,2,rev))
}

image(img_upright(img), col = grey(0:11/12))#グレーで表示


LDR <-function(Y,M,seed=1234){  #データをYとして受け取り、潜在変数次元をMとして受け取る
  set.seed(seed)#シード値のセット
  N <-ncol(Y)#Yの列数をN
  D <-nrow(Y)#Yの行数をD
  X <- matrix(rnorm(M*N),M,N)#要素を正規分布からのサンプルで初期化
  W <- matrix(rnorm(M*D),D,M)#要素を正規分布からのサンプルで初期化
  I_D <- diag(1,D)#対角成分が１の対角行列を用意
  S_muinv <- diag(1,D)
  S_Winv <- diag(1,M)
  I_M <- diag(1,M)
  pb <- txtProgressBar(min = 1, max = 1000, style = 3)
  for(i in 1:1000){
    S_muinv <-N*I_D+S_muinv#S_muinvを更新
    mu <-drop(rowSums(Y - W%*%X)%*%solve(S_muinv))#solveでS_muにする
    #rowSums(Y - W%*%X)%*%S_muでmu(μの平均値)を更新
    #計算結果がnum [1, 1:D]の配列になっているので
    #drop()でラベルを消してベクトルにする
    
    W <-((Y-mu)%*%t(X)) %*% solve(X%*%t(X)+S_Winv)#solveでS_Wにする
    #Wを更新
    
    
    
    X <-t(t(Y-mu)%*%W %*% (solve(t(W)%*%W+I_M)))#Xを更新
    
    setTxtProgressBar(pb, i)
  }
  list(W=W,X=X,mu=mu)
}


M10<-LDR(img_upright(img),10)#潜在変数次元10で試す
image((M10$W%*%M10$X+M10$mu),main="M=10", col = grey(0:11/12))


output <-function(img,n){#潜在変数次元数を引数に取り
  Mn<-LDR(img_upright(img),n)#潜在変数を抽出し
  image((Mn$W%*%Mn$X+Mn$mu),#画像を復元して表示する
        main=paste('M=',n), col = grey(0:11/12))
  
}

output(img,50)
