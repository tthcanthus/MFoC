library(fda)
library(compositions)
library(fdapace)
library(robCompositions)
library(maps)
library(mapdata)
library(viridis)
library(RColorBrewer)
#############################read data############################
pm25w <- read.table("pm25.txt")
covaria <- read.table("covariate.txt",header = T)
com.x <- covaria[,5:7]/apply(covaria[,5:7],1,sum)
x.ilr <- pivotCoord(com.x, pivotvar = 1)
covaria[,8] <- log(covaria[,8])
scovari <- apply(covaria[,8:11],2,scale)
Y <- apply(as.matrix(pm25w),2,scale)
Z <- as.matrix(cbind(x.ilr,scovari))
colnames(Z) <- c("x1","x2","x3","x4","x5","x6")
break_p <- 26
C <- 2
h <- 0.09
h1 <- 0.09
h2 <- 0.16
#############################set initial parameters######################
t <- seq(0, 1, length.out=break_p)
n <- dim(Y)[1]
epsilon = 1e-04
max_iter = 10000
max_restart = 25
ep_ker = function(x,h){
  a <- x
  index0 <- which(abs(x/h)>1)
  index1 <- which(abs(x/h)<=1)
  a[index1] <- 3/4*(1-(x[index1]/h)^2)*(1/h)
  a[index0] <- rep(0,length(index0))
  return(a)
}
B_coef <- matrix(0,6,break_p)
for(p1 in 1:break_p){
  t0 <- t[p1]
  kh <- as.matrix(ep_ker(t-t0,h),break_p,1)
  B_coef[,p1] <- Matrix::solve(t(Z) %*% Z*sum(kh), t(Z)%*%Y%*%kh)
}
{
  diff <- 1
  iter <- 0
  restart <- 0
  Y.list <- list()
  for(i in 1:C){
    Y.list[[i]] <- Y
  }
  # initialize (E step)
  WResult =  InitWList(Y,n, C,h)
  # first M step
  MResult = MUpdate(Y.list,WResult$WList,h)
  # second E step
  newW = EUpdate(Y.list,MResult$coef, MResult$sigma, MResult$Pai)
  W0 <- newW[[1]]
  llc <- newW$llc
}
#############################fit MFoC-i model######################
while (abs(diff) > epsilon && iter < max_iter && restart < max_restart) {
  MResult1 = MUpdate(Y.list,W0,h)
  newW1 = EUpdate(Y.list,MResult1$coef, MResult1$sigma, MResult1$Pai)
  newLL = newW1$llc
  if(sum(apply(W0,2,mean)>0.92)==1){
    W0 = InitWList(Y,n, C, h)[[1]]
    restart <- restart + 1
    ll <- -Inf
  } else {
    diff = newLL - llc
    W0 <- newW1[[1]]
    llc <- newLL
  }
  print(newLL)
}
#############################fit MFoC-c model######################
Y.star.c <- YUpdate(Y,MResult1$coef, newW1)
sig_ct0 <- list()
for(j in 1:C){
  sig_t0 <- rep(0,break_p)
  for(i in 1:break_p){
    t0 <- t[i]
    kh <- as.matrix(ep_ker(t-t0,h),break_p,1)
    fenzi <- diag(newW1[[1]][,j])%*%(Y.star.c[[j]]-Z%*%MResult1$coef[[j]])%*%
      diag(as.vector(kh))%*%t((Y.star.c[[j]]-Z%*%MResult1$coef[[j]]))
    sig_t0[i] <- sum(diag(fenzi))/sum(newW1[[1]][,j]*sum(kh))
  }
  sig_ct0[[j]] <- sig_t0
}
ll2 <- newLL
diff2 <- 10
while (abs(diff2) > 0.9 ) { 
  Y.star.c <- YUpdate(Y,MResult1$coef, newW1)
  ## E update
  Eresul2 <- EUpdate(Y.star.c,MResult1$coef,sig_ct0,MResult1$Pai)
  ## M Update
  Mresul2 <- MUpdate(Y.star.c,Eresul2[[1]],h)
  newLL2 = Eresul2$ll
  diff2 = newLL2 - ll2
  ll2 <- newLL2
  MResult1 <- Mresul2
  newW1 <- Eresul2
  sig_ct0 <- Mresul2$sigma
  
  print(t(c(newLL2,diff2)))
}
###############################display results###########################
w.est <- Eresul2[[1]]
apply(w.est,2,mean)
w.index <- w.est==apply(w.est,1,max)
png("figure/china clust.png",width=5000,height=5000)
par(mfrow=c(1,1))
map("china",col="darkgray",ylim=c(18,54),lwd=3)
points(covaria[which(w.index[,1]),3],covaria[which(w.index[,1]),4],pch=15,col= 2,cex=10)
points(covaria[which(w.index[,2]),3],covaria[which(w.index[,2]),4],pch=16,col= "dodgerblue",cex=10)
text(covaria[,3],covaria[,4],covaria[,2],cex=3.5)
dev.off()
#########################show the coefficient functions########################
colors5 <- c("#FF0000", "#00FF00", "#0000FF", "#FFA500", "#800080")
betad_1 <- pivotCoordInv(t(Mresul2$coef[[1]])[,1:2])
betafd_1 <- Data2fd(t,betad_1)
png("figure/betad1.png",width=500,height=400)
par(mar=c(4,6,1,1))
plot(betafd_1[1,],col=colors5[1],ylim=c(0.21,0.67),lwd=3,lty=1,
     ylab=expression(paste(symbol(b)[1]^D,"(t)")),
     cex.axis=2,cex.lab=2,xlab="t")
lines(betafd_1[2,],col=colors5[2],lwd=3,lty=2)
lines(betafd_1[3,],col=colors5[3],lwd=3,lty=3)
legend(x=0.45,y=0.67,legend=c("primary industry","secondary industry","tertiary industry"), 
       lty=c(1,2,4), col=c(colors5[1],colors5[2],colors5[3]),
       lwd=3,cex=1.5,horiz=FALSE)
dev.off()

min1 <- min(MResult1$coef[[2]])-0.1
max1 <- max(MResult1$coef[[2]])+0.1
betax1_1 <- as.matrix(pivotCoord(betad_1, pivotvar = 1))
betax2_1 <- as.matrix(pivotCoord(betad_1, pivotvar = 2))
betax3_1 <- as.matrix(pivotCoord(betad_1, pivotvar = 3))
png("figure/betax1.png",width=500,height=400)
par(mar=c(4,6,1,1))
plot(Data2fd(t,betax1_1)[1,],col=colors5[1],ylim=c(-0.5,1.2),lwd=4,lty=1,
     ylab=expression(paste(symbol(b)[1],"*","(t)")),
     cex.axis=2.2,cex.lab=2.2,xlab="t")
lines(Data2fd(t,betax2_1)[1,],col=colors5[2],lwd=4,lty=2)
lines(Data2fd(t,betax3_1)[1,],col=colors5[3],lwd=4,lty=3)
legend(x=0.45,y=1.2,legend=c("primary industry","secondary industry","tertiary industry"), 
       lty=c(1,2,4), col=c(colors5[1],colors5[2],colors5[3]),
       lwd=3,cex=1.5,horiz=FALSE)
lines(t,rep(0,length(t)),lty=4,col="gray",lwd=2)
t1 <- seq(0,0.096,by=0.001)
polygon(c(t1, rev(t1)), c(rep(max1,length(t1)), rep(min1, length(t1))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
t2 <- seq(0.096,0.346,by=0.001)
polygon(c(t2, rev(t2)), c(rep(max1,length(t2)), rep(min1, length(t2))),
        col =rgb(1, 0.498, 0.314,0.2), border = NA)
t3<- seq(0.346,0.604,by=0.001)
polygon(c(t3, rev(t3)), c(rep(max1,length(t3)), rep(min1, length(t3))),
        col =rgb(0.855, 0.647, 0.125,0.15), border = NA)
t4<- seq(0.604,0.857,by=0.001)
polygon(c(t4, rev(t4)), c(rep(max1,length(t4)), rep(min1, length(t4))),
        col =rgb(0.5294, 0.8078, 0.9216,0.2), border = NA)
t5 <- seq(0.857,1,by=0.001)
polygon(c(t5, rev(t5)), c(rep(max1,length(t5)), rep(min1, length(t5))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
dev.off()

betad_2 <- pivotCoordInv(t(Mresul2$coef[[2]])[,1:2])
betafd_2 <- Data2fd(t,betad_2)
png("figure/betad2.png",width=500,height=400)
par(mar=c(4,6,1,1))
plot(betafd_2[1,],col=colors5[1],ylim=c(0.09,0.78),lwd=3,lty=1,
     ylab=expression(paste(symbol(b)[2]^D,"(t)")),
     cex.axis=2,cex.lab=2,xlab="t")
lines(betafd_2[2,],col=colors5[2],lwd=3,lty=2)
lines(betafd_2[3,],col=colors5[3],lwd=3,lty=3)
legend(x=0.45,y=0.78,legend=c("primary industry","secondary industry","tertiary industry"), 
       lty=c(1,2,4), col=c(colors5[1],colors5[2],colors5[3]),
       lwd=3,cex=1.5,horiz=FALSE)
dev.off()

betax1_2 <- as.matrix(pivotCoord(betad_2, pivotvar = 1))
betax2_2 <- as.matrix(pivotCoord(betad_2, pivotvar = 2))
betax3_2 <- as.matrix(pivotCoord(betad_2, pivotvar = 3))
png("figure/betax2.png",width=500,height=400)
par(mar=c(4,6,1,1))
plot(Data2fd(t,betax1_2)[1,],col=colors5[1],ylim=c(-1.2,1.75),lwd=4,lty=1,
     ylab=expression(paste(symbol(b)[2],"*","(t)")),
     cex.axis=2.2,cex.lab=2.2,xlab="t")
lines(Data2fd(t,betax2_2)[1,],col=colors5[2],lwd=4,lty=2)
lines(Data2fd(t,betax3_2)[1,],col=colors5[3],lwd=4,lty=3)
legend(x=0.45,y=1.75,legend=c("primary industry","secondary industry","tertiary industry"), 
       lty=c(1,2,4), col=c(colors5[1],colors5[2],colors5[3]),
       lwd=3,cex=1.5,horiz=FALSE)
lines(t,rep(0,length(t)),lty=4,col="gray",lwd=2)
polygon(c(t1, rev(t1)), c(rep(1.88,length(t1)), rep(-1.3, length(t1))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
polygon(c(t2, rev(t2)), c(rep(1.88,length(t2)), rep(-1.3, length(t2))),
        col =rgb(1, 0.498, 0.314,0.2), border = NA)
polygon(c(t3, rev(t3)), c(rep(1.88,length(t3)), rep(-1.3, length(t3))),
        col =rgb(0.5294, 0.8078, 0.9216,0.15), border = NA)
polygon(c(t4, rev(t4)), c(rep(1.88,length(t4)), rep(-1.3, length(t4))),
        col =rgb(0.855, 0.647, 0.125,0.2), border = NA)
polygon(c(t5, rev(t5)), c(rep(1.88,length(t5)), rep(-1.3, length(t5))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
dev.off()

png("figure/gamma2.png",width=500,height=400)
par(mar=c(4,5,1,1))
plot(Data2fd(t,t(Mresul2$coef[[2]]))[3,],ylim=c(-0.35,0.55),lwd=4,col=colors5[2],lty=1,
     ylab="coefficients",
     cex.axis=2.2,cex.lab=2.2,xlab="t")
lines(Data2fd(t,t(Mresul2$coef[[2]]))[4,],lwd=4,col=colors5[3],lty=2)
lines(Data2fd(t,t(Mresul2$coef[[2]]))[5,],lwd=4,col=colors5[1],lty=3)
lines(Data2fd(t,t(Mresul2$coef[[2]]))[6,],lwd=4,col=colors5[5],lty=4)
legend(x=0.57,y=0.55,legend=c("per GDP","tenperature","humidity","wind"), 
       lty=c(1,2,3,4), col=c(colors5[2],colors5[3],colors5[1],colors5[5]),
       lwd=3,cex=1.5,horiz=FALSE)
lines(t,rep(0,length(t)),lty=4,col="gray",lwd=2)
polygon(c(t1, rev(t1)), c(rep(1.88,length(t1)), rep(-1.3, length(t1))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
polygon(c(t2, rev(t2)), c(rep(1.88,length(t2)), rep(-1.3, length(t2))),
        col =rgb(1, 0.498, 0.314,0.2), border = NA)
polygon(c(t3, rev(t3)), c(rep(1.88,length(t3)), rep(-1.3, length(t3))),
        col =rgb(0.5294, 0.8078, 0.9216,0.15), border = NA)
polygon(c(t4, rev(t4)), c(rep(1.88,length(t4)), rep(-1.3, length(t4))),
        col =rgb(0.855, 0.647, 0.125,0.2), border = NA)
polygon(c(t5, rev(t5)), c(rep(1.88,length(t5)), rep(-1.3, length(t5))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
dev.off()

png("figure/gamma1.png",width=500,height=400)
par(mar=c(4,5,1,1))
plot(Data2fd(t,t(Mresul2$coef[[1]]))[3,],ylim=c(-0.68,0.98),lwd=4,col=colors5[2],lty=1,
     ylab="coefficients",
     cex.axis=2.2,cex.lab=2.2,xlab="t")
lines(Data2fd(t,t(Mresul2$coef[[1]]))[4,],lwd=4,col=colors5[3],lty=2)
lines(Data2fd(t,t(Mresul2$coef[[1]]))[5,],lwd=4,col=colors5[1],lty=3)
lines(Data2fd(t,t(Mresul2$coef[[1]]))[6,],lwd=4,col=colors5[5],lty=4)
legend(x=0.57,y=0.98,legend=c("per GDP","tenperature","humidity","wind"), 
       lty=c(1,2,3,4), col=c(colors5[2],colors5[3],colors5[1],colors5[5]),
       lwd=3,cex=1.5,horiz=FALSE)
lines(t,rep(0,length(t)),lty=4,col="gray",lwd=2)
polygon(c(t1, rev(t1)), c(rep(1.88,length(t1)), rep(-1.3, length(t1))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
polygon(c(t2, rev(t2)), c(rep(1.88,length(t2)), rep(-1.3, length(t2))),
        col =rgb(1, 0.498, 0.314,0.2), border = NA)
polygon(c(t3, rev(t3)), c(rep(1.88,length(t3)), rep(-1.3, length(t3))),
        col =rgb(0.5294, 0.8078, 0.9216,0.15), border = NA)
polygon(c(t4, rev(t4)), c(rep(1.88,length(t4)), rep(-1.3, length(t4))),
        col =rgb(0.855, 0.647, 0.125,0.2), border = NA)
polygon(c(t5, rev(t5)), c(rep(1.88,length(t5)), rep(-1.3, length(t5))),
        col =rgb(0.596, 0.851, 0.392,0.2), border = NA)
dev.off()

png("figure/pm25clust2.png",height =400,width = 500)
par(mar = c(4, 5, 1, 1),mfrow=c(1,1))
plot(Data2fd(1:26,t(as.matrix(pm25w[w.index[,2],]))),lwd=2,ylim=c(0,180),
     ylab="PM2.5",xlab="t",col=4,cex.axis=2,cex.lab=2)
dev.off()

png("figure/pm25clust1.png",height =400,width = 500)
par(mar = c(4, 5, 1, 1),mfrow=c(1,1))
plot(Data2fd(1:26,t(as.matrix(pm25w[w.index[,1],]))),lwd=2,ylim=c(0,180),
     ylab="PM2.5",xlab="t",col=2,cex.axis=2,cex.lab=2)
dev.off()
