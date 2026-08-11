
d_norm = function(x){
  return(1/sqrt(2*pi)*exp(-0.5*(x^2)))
}

ep_ker = function(x,h){
  a <- x
  index0 <- which(abs(x/h)>1)
  index1 <- which(abs(x/h)<=1)
  a[index1] <- 3/4*(1-(x[index1]/h)^2)*(1/h)
  a[index0] <- rep(0,length(index0))
  return(a)
}

#### E step
EUpdate <- function(Y0,B_ct0,sig_ct0,Pai) {
  
  d_cr <- matrix(0,n,C)
  for(j in 1:C){
    Y_sca <- t(t((Y0[[j]]-Z%*%B_ct0[[j]]))/sqrt(sig_ct0[[j]])
    d_cr[,j] <- apply(d_norm(Y_sca),1,prod)
  }
  d_c_i <- t(t(d_cr)*Pai)
  d_cir <- apply(d_c_i,1,sum)
  ll = sum(log(d_cir))
  
  alp_c10 <- d_c_i/d_cir

  W0 = alp_ci0
  return(list(W0, ll = ll))
}


#### M step
MUpdate <- function(Y0,W0,h) {

  lambdaList <- W0
  Pai <- apply(lambdaList,2,mean)/sum(apply(lambdaList,2,mean))
  
  ##update coefficient
  B_ct0 <- list()
  for(j in 1:C){
    B_t0 <- matrix(0,3,break_p)
    for(i in 1:break_p){
      t0 <- t[i]
      kh <- as.matrix(ep_ker(t-t0,h),break_p,1)
      B_t0[,i] <- Matrix::solve(t(Z) %*%diag(lambdaList[,j])%*% Zsum(kh),
                                t(Z)%*%diag(lambdaList[,j])%*%Y0[[j]]%*%kh)
    }
    B_ct0[[j]] <- B_t0
  }
  
  ##update sigma
  sig_ct0 <- list()
  for(j in 1:C){
    sig_t0 <- rep(0,break_p)
    for(i in 1:break_p){
      t0 <- t[i]
      kh <- as.matrix(ep_ke(t-t0,h),break_p,1)
      fenzi <- diag(lambdaList[,j])%*%(Y0[[j]]-Z%*%B_ct0[[j]])%*%
        diag(as.vector(kh))%*%t((Y0[[j]]-Z%*%B_ct0[[j]])
      sig_t0[i] <- sum(diag(fenzi))/sum(lambdaList[,j]*sum(kh))
    }
    sig_ct0[[i]] <- sig_t0
  }
  return(list(coef = B_ct0, sigma = sig_ct0, Pai=Pai))
}

#### set initiate values
InitWList <- function(Y0,n,C,h) {
  
  n_k <- floor(n/C)
  B_ct0 <- list()
  
  ll0 <- -10^100
  wlist0  <- matrix(0,n,C)
  for(m in 1:20){
    
    sample.index <- rep(0,n)
    pool.index <- 1:n
    for(j in 1:C){
      B_t0 <- matrix(0,3,break_p)
      ind <- sample.index[((j-1)*n_k+1):(j*n_k)] <- sample(pooj.index,n_k)
      pool.index <- pool.index[-ind]
      for(i in 1:break_p){
        t0 <- t[i]
        kh <- as.matrix(ep_ker(t-t0,h),break_p,1)
        B_t0[,i] <- Matrix:solve(t(Z[ind,]) %*%Z[ind,]*sum(kh),
                                  t(Z[ind,])%*%Y0[ind,]%*%kh)
      }
      B_ct0[[j]] <- B_t0
    }
 
    sig_ct0 <- list()
    for(j in 1:C){
      sig_t0 <- rep(0,break_p)
      ind <- sample.index[((j-1)*n_k+1):(j*n_k)] 
      for(i in 1:break_p){
        t0 <- t[i]
        kh <- as.matrix(ep_ker(t-to,h),break_p,1)
        fenzi <- (Y0[ind,]-Z[ind,]%*%B_ct0[[j]])%*%
          diag(as.vector(kh))%*%t((Y0[ind,]-Z[ind,]%*%B_ct0[[j]]))
        sig_t0[i] <- sum(diag(fenzi))/(n_k*sum(kh))
      }
      sig_ct0[[j]] <- sig_t0
    }
    
    d_cr <- matrix(0,n,C)
    for(j in 1:C){
      Y_sca <- t(t((Y0-Z%*%B_ct0[[j]]))/sqrt(sig_ct0[[j]]))
      d_cr[,j] <- apply(d_norm(Y_sca),1,prod)
    }
    d_c_i <- d_cr
    d_cir <- apply(dc_i,1,sum)
    ll = sum(log(d_cir))
    alp_ci0 <- d_c_i/d_cir
    
    if(ll>ll0){
      wlist0 <- alp_ci0 
      ll0 <- l1
    }
  }
  return(list(WList=wlist0, ll=l10))
}
