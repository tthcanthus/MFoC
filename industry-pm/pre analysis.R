library(fda)
library(compositions)
library(cluster)
library(robCompositions)
library(viridis)
library(RColorBrewer)
library(flexmix)
library(nlme)
##############visualization of original data###########
pm25d <- read.table("pm25.txt")
png("figure/pm25-raw.png",height =400,width = 1000)
par(mar = c(4, 5, 1, 1))
set1_colors <- brewer.pal(12, "Set3")
my_colors <- rep(set1_colors, ceiling(217/12))[1:217]
plot(1:26,as.matrix(pm25w[1,]),lwd=0.6,ylim=c(5,190),
     ylab="PM2.5",xlab="biweek",cex.axis=1.6,cex.lab=1.6,type="b",col=my_colors[1],pch=19)
for(j in 2:220){
  lines(1:26,as.matrix(pm25w[j,]),lwd=0.6,
        ylab="PM2.5",cex.axis=1.6,cex.lab=1.6,type="b",col=my_colors[j],pch=19)
}
dev.off()

png("figure/pm25o.png",height = 700,width = 1000)
par(mar = c(4, 5, 1, 1))
plot(Data2fd(1:26,t(as.matrix(pm25w))),lwd=1.7,
     ylab="PM2.5",xlab="week",cex.axis=2,cex.lab=2)
dev.off()

covaria <- read.table("industry and others.txt",header = T)
com.x <- covaria[,5:7]/apply(covaria[,5:7],1,sum)
png("figure/hengjiemian.png",height =400,width = 500,type = "cairo",res=120)
par(mar = c(1, 2.9, 0.0001, 4))
plot.acomp(com.x,col=4,cex=0.9,pch=16,labels = c("primary","secondary","tertiary"),
           cex.lab=2)
dev.off()
#################mixture of scalar-on-composition to fit the data##############
com.ilr <- pivotCoord(com.x,pivotvar = 1)
data0 <- data.frame(com.ilr ,pm25w[,18])
colnames(data0) <- c("x1","x2","Y")
K <- 1:10
fits <- lapply(K, function(k) {
  flexmix(Y~x1+x2,data=data0, k = k)
})
bic_vals <- sapply(fits, BIC)
data.frame(K, BIC = bic_vals)

het1 = gls(Y~x1+x2,data = data0,
           weights = varExp(form =  ~ x1+x2)
)
BIC(het1)
####################triangular prim########################################
mix1 <- flexmix(Y~x1+x2,data=data0, k = 2)
clus1 <- which(mix1@cluster==1)
clus2 <- which(mix1@cluster==2)
l=18
open3d(windowRect = c(100, 100, 480, 400))
llay <- ceiling(min(pm25w[,l]))-2
ulay <- floor(max(pm25w[,l]))+2
# define the base
v1 <- c(0, 0, llay)
v2 <- c(55, 0, llay)
v3 <- c(55*0.5, 55*(sqrt(3)/2), llay)
# define the top
v4 <- c(0, 0, ulay)
v5 <- c(55, 0, ulay)
v6 <- c(55*0.5, 55*(sqrt(3)/2), ulay)
text3d((v1-0.2), texts = "primary",adj = c(0.7, 0.3),
       cex = 0.8)
text3d((v2+0.4), texts = "secondary",adj = c(-0.02, 0.1),
       cex = 0.8)
text3d(v3-0.3, texts = "tertiary",adj = c(0.7, 0.7),
       cex = 0.8)
v7 <- 0.5*(v3+v6)
text3d(v7[1]-0.4, v7[2]-0.4,v7[3],texts = "PM2.5",adj = c(1.5, 0),
       cex = 0.8)
triangles3d(rbind(v1, v2, v3), color="blue", alpha=0.05)
triangles3d(rbind(v4, v5, v6), color="blue", alpha=0.05)
triangles3d(rbind(v1, v2, v5), color="yellow", alpha=0.05)
triangles3d(rbind(v1, v5, v4), color="yellow", alpha=0.05)
triangles3d(rbind(v2, v3, v6), color="yellow", alpha=0.05)
triangles3d(rbind(v2, v6, v5), color="yellow", alpha=0.05)
triangles3d(rbind(v3, v1, v4), color="yellow", alpha=0.05)
triangles3d(rbind(v3, v4, v6), color="yellow", alpha=0.05)

segments3d(rbind(v3,v6,v1,v3,v4,v6,v2,v3,v5,v6), col="gray", lwd=1.1)
dir <- v6-v3 
dir <- dir / sqrt(sum(dir^2))  
#find perpendicular direction
perp <- c(-dir[2], dir[1], 0)
if(sum(perp^2) < 1e-6) perp <- c(0, -dir[3], dir[2])
perp <- perp / sqrt(sum(perp^2))
for(i in 0:5) {
  t <- 0.2*i
  p <- v3 + t * (v6 - v3)
  segments3d(rbind(p, p + 4*perp), lwd = 1.2)
  
  text3d(p[1], p[2], p[3],
         texts = sprintf("%.0f", (llay+(ulay-llay)*t)),
         adj = c(1.4, 0), cex = 0.75)
}

cordi.comx <- matrix(0,dim(com.x)[1],3)
cordi.comx[,1] <- 55*com.x[,2]+55*com.x[,3]*0.5
cordi.comx[,2] <- 55*com.x[,3]*0.5*sqrt(3)
cordi.comx[,3] <- pm25w[,l]
spheres3d(cordi.comx[clus1,1], cordi.comx[clus1,2], cordi.comx[clus1,3],col="red",radius = 0.8)
points3d(cordi.comx[clus2,1], cordi.comx[clus2,2], cordi.comx[clus2,3],col="blue",size=6.5)
rgl.snapshot("figure/triangular prim.png")

