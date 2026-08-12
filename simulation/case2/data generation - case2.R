

break_p <- 20
##generate coefficients
x <- seq(0, 1, length.out=break_p)
n <- 200
pi0 <- 3.1415926

data11 <- 2*exp(-x)
data12 <- cos(0.5*pi0 * x) + 1
data13 <- -0.5*x+1.5

total1 <- data11 + data12 + data13
data11 <- data11 / total1
data12 <- data12 / total1
data13 <- data13 / total1
data_beta1 <- cbind(data11, data12, data13)
data_beta_ilr1 <- as.matrix(ilr(data_beta1))
gama_data1 <- 0.2*cos(1.5*pi0*x)

data21 <- exp(x)
data22 <- sin(pi0 * x) + 1.5
data23 <- 0.1*x^(1/3)+1.2
total2 <- data21 + data22 + data23
data21 <- data21 / total2
data22 <- data22 / total2
data23 <- data23 / total2

data_beta2 <- cbind(data21, data22, data23)
data_beta_ilr2 <- as.matrix(ilr(data_beta2))
gama_data2 <- 0.2*sin(2*pi0*x)

##generate compositional predictor
x.mean <- acomp(c(3, 2, 1))
x.var <- ilrvar2clr(matrix(c(1, -0.5, -0.5, 1), ncol = 2))
x_com <- rnorm.acomp(n, mean = x.mean, var = x.var)
x_ilr <- as.matrix(ilr(x_com))

z <- rnorm(n, 0.5, 0.5)

e <- rnorm(n, 0, 0.1)

index <- sample(c(0,1),n,replace = TRUE,prob=c(0.45,0.55))
index1 <- which(index==0)
index2 <- which(index==1)

n1 <- length(index1)
n2 <- length(index2)


####### generate functional response in Case2
error12  <- kronecker(rnorm(n1,0,sqrt(0.02)),sin(4*pi0*x)+
  kronacker(rnorm(n1,0,sqrt(0.01)),cos(4*pi0*x))
error12 <- t(matrix(error12,break_p,n1))+e[index1]

error22  <- kronecker(rnorm(n2,0,sqrt(0.02)),sin(3*pi0*x)+
  kronecker(rnorm(n2,0,sqrt(0.01)),cos(3*pi0*x))
error22 <- t(matrix(error22,break_p,n2))+e[index2]

Y2 <- matrix(0,n,break_p)
sig21 <- x_ilr[index1,] %*% t(data_beta_ilr1) + matrix(z[index1] %X%gama_data1,break_p,n1)
sig22 <- x_ilr[index2,] %*% t(data_beta_ilr2) + matrix(z[index2] %x%gama_data2,break_p,n2)
Y2[index1,] <- sig21 + error12
Y2[index2,] <- sig22 + error22

