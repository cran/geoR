##
## An example of a geostatistical analysis using the package geoR
## ==============================================================
##
## Data: Hydraulic conductivity measurements
##
## files: (they must be copied to your working directory)
##    - Cruciani.dat     (coordinates of the data locations and K measurements)
##    - Cruciani.borders (coordinates defining the limits of the area)
##
require(geoR)
##
## 1. Reading data
##
cru <- read.geodata("Cruciani.dat", head=T, coords.col=2:3, data.col=4)
cru.borda <- read.table("Cruciani.border", head=T)[,2:3]
cru.borda <- rbind(cru.borda, cru.borda[1,])
##
## 2. Exploratory plots
##
## Plotting original data 
plot(cru, bord=cru.borda)
points(cru, bord=cru.borda, cex.min=1, cex.max=1, col=gray(seq(0.9,0,l=length(cru$data))), xla="Coord X", ylab="Coord Y")
##
## now transforming the data (log)
plot(cru, bord=cru.borda, lambda = 0)
points(cru, bord=cru.borda, cex.min=1, cex.max=1, lambda = 0, col=gray(seq(0.9,0,l=length(cru$data))), xlab="Coord X", ylab="Coord Y")
##
## histograms for original and transformed data
par.ori <- par(no.readonly = TRUE)

par(mfrow=c(1,2), mar=c(3.5,3.5,1,1))
hist(cru$data, main="", xlab="K data")
hist(log(cru$data), main="", xlab="log(K) data")
par(par.ori)
##
## 3. Parameter estimation via maximum likelihood
##
## defining a grid of initial values for the numerical maximisation of the likelihood
ini.fit <- expand.grid((0:12)/2, (0:12)/2)
##
## Estimating parameters, including the transformation parameter lambda, under different correlation models
##
cru.ml <- list()
##
cru.ml$mat1 <- likfit(cru, ini=ini.fit, fix.nugget=F, fix.lambda=F, cov.model="mat", kappa=1)
summary(cru.ml$mat1)
cru.ml$gau <- likfit(cru, ini=ini.fit, fix.nugget=F, fix.lambda=F, cov.model="gau")
summary(cru.ml$gau)
cru.ml$sph <- likfit(cru, ini=ini.fit, fix.nugget=F, fix.lambda=F, cov.model="sph")
summary(cru.ml$sph)
cru.ml$exp <- likfit(cru, ini=ini.fit, fix.nugget=F, fix.lambda=F, cov.model="exp")
summary(cru.ml$exp)
##
## Now fixing lambda = 0
##
cru.ml0 <- list()
##
cru.ml0$mat1 <- likfit(cru, ini=ini.fit, fix.nugget=F, lambda=0, cov.model="mat", kappa=1)
summary(cru.ml0$mat1)
cru.ml0$gau <- likfit(cru, ini=ini.fit, fix.nugget=F, lambda=0, cov.model="gau")
summary(cru.ml0$gau)
cru.ml0$sph <- likfit(cru, ini=ini.fit, fix.nugget=F, lambda=0, cov.model="sph")
summary(cru.ml0$sph)
cru.ml0$exp <- likfit(cru, ini=ini.fit, fix.nugget=F, lambda=0, cov.model="exp")
summary(cru.ml0$exp)
##
## 4. Profile likelihoods
##
#prof.lambda <- proflik(cru.ml$gau,cru,sill.val=NULL, range.val=NULL, nugget.val=NULL, lambda.val=seq(-.4,.4,l=21))
#plot(prof.lambda)
##
#cru.prof <- list()
#cru.prof$ml.gau <- proflik(cru.ml$gau,cru,sill.val=seq(.5, 4, l=21), range.val=seq(.8, 5, l=21), nugget.val=seq(0,2,l=21), lambda.val=seq(-.4,.4,l=21))
#par(mfrow=c(2,2))
#plot(cru.prof$ml.gau)
#par(mfrow=c(1,1))
##
#cru.prof$ml0.gau <- proflik(cru.ml0$gau,cru,sill.val=seq(.5, 4, l=21), range.val=seq(.8, 5, l=21), nugget.val=seq(0,1.5,l=21))
#par(mfrow=c(2,2))
#plot(cru.prof$ml0.gau)
#par(mfrow=c(1,1))
##
## 5. Empirical variogram
##
cru.vario <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
plot(cru.vario)
lines(cru.ml0$mat1)
lines(cru.ml0$gau, col="blue")
lines(cru.ml0$sph, col="red")
lines(cru.ml0$exp, col="green")
##
## 6. Kriging
##
apply(cru.borda,2,range)
cru.grid <- expand.grid(seq(0,22.5,by=0.125), seq(0,13,by=0.125))

cru.sk <- krige.conv(cru, loc=cru.grid, krige=krige.control(cov.model="gau", obj.model=cru.ml0$gau, lambda=0, type="sk"), border=cru.borda)

image(cru.sk, border=cru.borda, loc=expand.grid(seq(0,22.5,by=0.125), seq(0,13,by=0.125)), col=gray(seq(.9,0,l=21)), x.leg=c(3,18), y.leg=c(-3, -2), coords.data=cru$coords, xlab="Coord X", ylab="Coord Y")

contour(cru.sk, coords.data=cru$coords, xlab="Coord X", ylab="Coord Y")

contour(cru.sk, coords.data=cru$coords, xlab="Coord X", ylab="Coord Y", filled=TRUE)

image(cru.sk, col=gray(seq(.9,0,l=21)), x.leg=c(3,18), y.leg=c(-3, -2), coords.data=cru$coords, xlab="Coord X", ylab="Coord Y")

image(cru.sk, values=cru.sk$krige.var, col=gray(seq(.9,0,l=21)), x.leg=c(3,18), y.leg=c(-3, -2), coords.data=cru$coords, xlab="Coord X", ylab="Coord Y")
##
## 7. Difficulties when using variograms for inference on model parameters
##     (just for illustration)
##
cru.v1 <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
cru.v1.exp <- variofit(cru.v1, ini=ini.fit, cov.model="exp")
cru.v1.gau <- variofit(cru.v1, ini=ini.fit, cov.model="gau")
cru.v2 <- variog(cru, uvec=seq(0,12, l=8), lambda=0)
cru.v2.exp <- variofit(cru.v2, ini=ini.fit, cov.model="exp")
cru.v2.gau <- variofit(cru.v2, ini=ini.fit, cov.model="gau")
cru.v3 <- variog(cru, uvec=seq(1,15, l=8), lambda=0)
cru.v3.exp <- variofit(cru.v3, ini=ini.fit, cov.model="exp")
cru.v3.gau <- variofit(cru.v3, ini=ini.fit, cov.model="gau")
cru.v4 <- variog(cru, uvec=seq(0,9, l=12), lambda=0)
cru.v4.exp <- variofit(cru.v4, ini=ini.fit, cov.model="exp")
cru.v4.gau <- variofit(cru.v4, ini=ini.fit, cov.model="gau")
cru.v5 <- variog(cru, uvec=seq(0,12, l=12), lambda=0)
cru.v5.exp <- variofit(cru.v5, ini=ini.fit, cov.model="exp")
cru.v5.gau <- variofit(cru.v5, ini=ini.fit, cov.model="gau")
cru.v6 <- variog(cru, uvec=seq(1,15, l=12), lambda=0)
cru.v6.exp <- variofit(cru.v6, ini=ini.fit, cov.model="exp")
cru.v6.gau <- variofit(cru.v6, ini=ini.fit, cov.model="gau")
##
par(mfrow=c(3,2), mar=c(3,3,.5,.5), mgp=c(1.8,.8,0))
plot(cru.v1, ylim=c(0, 3.6))
lines(cru.v1.exp)
lines(cru.v1.gau)
plot(cru.v2, ylim=c(0, 3.6))
lines(cru.v2.exp)
lines(cru.v2.gau)
plot(cru.v3, ylim=c(0, 3.6))
lines(cru.v3.exp)
lines(cru.v3.gau)
plot(cru.v4, ylim=c(0, 3.6))
lines(cru.v4.exp)
lines(cru.v4.gau)
plot(cru.v5, ylim=c(0, 3.6))
lines(cru.v5.exp)
lines(cru.v5.gau)
plot(cru.v6, ylim=c(0, 3.6))
lines(cru.v6.exp)
lines(cru.v6.gau)
par(par.ori)
