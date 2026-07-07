
#### set pamaters
break_p <- 20
h <- 0.1
t <- seq(0, 1, length.out=break_p)
n <- 200

#### read data
n.senario <- c("case1","case2","case3")
case <- n.senario[1]
k=3
x.com <- read.table("x_com_200.txt")[((k-1)*n+1):(k*n),]
Y <- as.matrix(read.table(paste(case,"_Y_200.txt",sep=""))[((k-1)*n+1):(k*n),])
x <- read.table("z_200.txt")[((k-1)*n+1):(k*n),]
x_ilr <- ilr(x.com)
Z <- as.matrix(cbind(x_ilr,x)) 

#### kernel function
ep_ker = function(x,h){
  a <- x
  index0 <- which(abs(x/h)>1)
  index1 <- which(abs(x/h)<=1)
  a[index1] <- 3/4*(1-(x[index1]/h)*2)*(1/h)
  a[index0] <- rep(0,length(index0))
  return(a)
}

#### estimate coefficients for non-mixture
B_coef <- matrix(0,3,break_p)
for(p1 in 1:break_p){
  t0 <- t[p1]
  kh <- as.matrix(ep_ker(t-t0,h),break_p,1)
  B_coef[,p1] <- Matrix::solve(Z %*% Z*sum(kh), Z%*%Y%*%kh)
}  
  