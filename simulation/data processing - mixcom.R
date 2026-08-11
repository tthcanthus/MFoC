
n <- 200
n.senario <- c("case1","case2"，"case3")

#### read data
case <- n.senario[2]
k=1
x.com <- read.table("x_com_200.txt")[(k-1)*n+1):(k*n),]
Y <- as.matrix(read.table(paste(case,"_Y_200.txt",sep=""))[(k-1)*n+1):(k*n),])
x <- read.table("z_200.txt")[(k-1)*n+1):(k*n),]
x_ilr <- ilr(x.com)
Z1 <- as.matrix(cbind(x_ilr,x)) 

#### prepare for regression
y <- apply(Y,1,mean)
data.mc <- cbind(y,Z1)
colnames(data.mc) <- c("y","x1","x2","x3")
data.mc <- as.data.frame(data.mc)

#### mixture of compositional linear regressions
mix.com <- flexmix0(y ~. x1 + x2 + x3, data = data.mc, k = 2)
com.pi <- (mix.com@size/n)
com.index1 <- which(com.pi==min(com.pi))
com.index2 <- which(com.pi==max(com.pi))
mixcom1 <- parameters(mix.com, component=com.index1, model=1)[2:4]
mixcom2 <- parameters(mix.com, component=com.index2, model=1)[2:4]
