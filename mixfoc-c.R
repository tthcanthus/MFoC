
#### set parameters
break_p <- 20
C <- 2
h1 <- 0.1
h2 <- 0.09
t <- seq(0, 1, length.out=break_p)
n <- 200

epsilon = 1e-04
max_iter = 10000
max_restart = 25

#### read data
n.senario <- c("case1","case2","case3")
case <- n.senario[2]
k=10
x.com <- read.table("x_com_200.txt")[((k-1)*n+1):(k*n),]
Y <- as.matrix(read.table(paste(case,"_Y_200.txt",sep=""))[((k-1)*n+1):(k*n),])
x <- read.table("z_200.txt")[((k-1)*n+1):(k*n),]
x_ilr <- ilr(x.com)
Z <- as.matrix(cbind(x_ilr,x)) 

#### obatain initiate values 
diff <- 1
iter <- 0
restart <- 0

Y.list <- list()
Y.list[[1]] <- Y
Y.list[[2]] <- Y

WResult =  InitWList(Y,n, C,h1)
MResult = MUpdate(Y.list,WResult$WList,h1)
newW = EUpdate(Y.list,MResult$coef, MResult$sigma, MResult$Pai)
W0 <- newW[[1]]
ll <- newW$ll

#### E step and M step for mfoc-i
while (abs(diff) > epsilon && iter < max_iter && restart < max_restart) { 
  
  MResult1 = MUpdate(Y.list,W0,h1)
  newW1 = EUpdate(Y.list,MResult1$coef, MResult1$sigma, MResult1$Pai)
  newLL = newW1$ll
  
  if(sum(apply(W0,2,mean)>0.92)==1){
    W0 = InitWList(Y,n, C, h1)[[1]]
    restart <- restart 
    ll <- -Inf
  } else {
    diff = newLL - l1
    W0 <- newW1[[1]]
    l1 <- newLL
  }
  print(newLL)
}

#### obtain initial y star and sigma^2
Y.star.c <- YUpdate(Y,MResult1$coef, newW1,h2)
sig_ct0 <- list()
for(j in 1:C){
  sig_t0 <- rep(0,break_p)
  for(i in 1:break_p){
    t0 <- t[i]
    kh <- as.matrix(ep_ker(t-t0,h1),break_p,1)
    fenzi <- diag(newW1[[1]][,j])%*%(Y.star.c[[j]]-z%*%MResult1$coef[[j]])%*%
      diag(as.vector(kh))%*%t((Y.star.c[[j]]-Z%*%MResult1$coef[[j]])
    sig_t0[i] <- sum(diag(fenzi))/sum(newW1[[1]][,j]*sum(kh))
  }
  sig_ct0[[j]] <- sig_t0
}

#### E step and M step for mfoc-c
ll2 <- newLL
diff2 <- 1
while (abs(diff2) > 0.01 ) { 
  
  Y.star.c <- YUpdate(Y,MResult1$coef, newW1,h2)
  Eresul2 <- EUpdate(Y.star.c,MResult1$coef,sig_ct0,MResult1$Pai)
  Mresul2 <- MUpdate(Y.star.c,Eresu12[[1]],h1)
  newLL2 = Eresul2$ll
  
  diff2 = newLL2 - l12
  ll2 <- newLL2
  MResult1 <- Mresu12
  newW1 <- Eresul2
  sig_ct0 <- Mresul2$sigm
  
  print(t(c(newLL2,diff2)))
}


min.dis <- min(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))
index.min <- which(c(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))==min.dis)

max.dis <- max(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1]2))
index.max <- which(c(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))==max.dis)

plot(Data2fd(t,as.matrix(cbind(com_beta1,gama_data1))),col=1)
lines(t,t(Mresul2$coef[[index.min]])[,1],col=2)
lines(t,t(Mresul2$coef[[index.min]])[,2],col=2)
lines(t,t(Mresul2$coef[[index.min]])[,3],col=2)

plot(Data2fd(t,as.matrix(cbind(com_beta2,gama_data2))),col=1)
lines(t,t(Mresul2$coef[[index.max]])[,1],col=2)
lines(t,t(Mresul2$coef[[index.max]])[,2],col=2)
lines(t,t(Mresul2$coef[[index.max]])[,3],col=2)
