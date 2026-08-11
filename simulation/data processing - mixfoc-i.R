
library(RcppEigen)
library(compositions)

#### set parameters
break_p <- 20
C <- 2
h1 <- 0.1
t <- seq(0, 1, length.out=break_p)
n <- 200

epsilon = 1e-04
max_iter = 10000
max_restart = 25

#### read data
n.senario <- c("case1","case2","case3")
case <- n.senario[1]
k=5
x.com <- read.table("x_com_200.txt")[(k-1)*n+1):(k*n),]
Y <- as.matrix(read.table(paste(case,"_Y_200.txt",sep=""))[(k-1)*n+1):(k*n),])
x <- read.table("z_200.txt")[(k-1)*n+1):(k*n),]
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
ll <- newW$l1

#### E step and M step for loop
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

#### plot the estimated coefficients
min.dis <- min(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))
index.min <- which(c(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))==min.dis)

max.dis <- max(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))
index.max <- which(c(sum((MResult1$coef[[1]][1,]-com_beta1[,1])),sum((MResult1$coef[[2]][1,]-com_beta1[,1])))==max.dis)

plot(Data2fd(t,as.matrix(cbind(com_beta1,gama_data1))),
     col=1,ylim=c(-1,1))
lines(t,MResult1$coef[[index.min]][,1],col=2)
lines(t,MResult1$coef[[index.min]][,2],col=2)
lines(t,MResult1$coef[[index.min]][,3],col=2)

plot(Data2fd(t,as.matrix(cbind(com_beta2,gama_data2))),
     col=1,ylim=c(-1,1))
lines(t,MResult1$coef[[index.max]][,1],col=2)
lines(t,MResult1$coef[[index.max]][,2],col=2)
lines(t,MResult1$coef[[index.max]][,3],col=2)
