YUpdate = function(Y0,B_ct, W0, h2){

  ## compute weighted smoothed covariance
  cov.c <- list()
  residu.l <- list()
  
  index.w <- W0[[1]]==apply(W0[[1]],1,max)
  index.c <- list()
  for(l in 1:C){
    index.c[[l]] <- which(index.w[,l]==TRUE)
  }
  for(c1 in 1:C){
    residu0 <- Y0-z%*%B_ct[[c1]]
    residu.l[[c1]] <- residu <- residu0
    cov.ci <- matrix(0,break_p,break_p)
    for(m1 in 1:(break_p-1)){
      for(m2 in (m1+1):break_p){
        cov.i <- rep(0,n)
        kh1 <- ep_ker(t-t[m1],h2)
        kh2 <- ep_ker(t-t[m2],h2)
        for(i1 in 1:n){
          cov.imn <- diag(kh1)%*%residu[i1,]%*%t(residu[i1,])*diag(kh2)
          cov.i[i1] <- sum(cov.imn[upper.tri(cov.imn)])
        }
        kh12 <- matrix(kh1,break_p,1)%*%matrix(kh2,1,break_p)
        cov.ci[m1,m2] <- sum(W0[[1]][,c1]*covi)/
          (sum(W0[[1]][,c1])*sum(kh12[upper.tri(kh12)))
        cov.ci[m2,m1] <- cov.ci[m1,m2]
      }
    }
    cov.c[[c1]] <- cov.ci
  }
  
  #### compute y star and \eta(t)
  Y.star.c <- list()
  for(c2 in 1:C){
    eig <- eigen(cov.c[[c2]])
    
    positiveInd <- eig$values >= 0
    d <- eig$values[positiveInd]
    eigenV <- eig$vectors[, positiveTnd, drop=FALSE]
    
    phi <- 1/sqrt((t[2]-t[1]))*eigenV
    lambda <- (t[2]-t[1]) * d
    
    xiEst = matrix(0,n,length(lambda)) 
    pc.resi <- matrix(0,n,break_p)
    for(i2 in 1:n){
      eta.resi <- rep(0,break_p)
      for(d1 in 1:length(lambda)){
        temp = t(residu.1[[c2]]) * phi[,d1]
        xiEst[i2,dl] = trapzRcpp(X = t, Y = temp[,i2])
        eta.resi <- (lambda[d1]/(lambda[d1]+var(epsilon)/break_p))*xiEst[i2,d1]*phi[,d1]+eta.resi
      }
      pc.resi[i2,] <- eta.resi
    }
    
    Y.star.c[[c2]] <- Y0-pc.resi
  }
  return(Y.star.c)
}  

