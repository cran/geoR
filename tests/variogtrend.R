require(geoR)
options(digits = 3, width = 80)
set.seed(44)
##
## Simulating data
##
ap <- grf(100, cov.pars=c(10, .25))
##
## replicating data and
## adding a "trend" to the simulated data
##
apt <- ap
apt$data <- 50 + apt$data - 15*ap$coords[,1] + 20*ap$coords[,2]
##
## creating a covariate
##
apt$covar <- 15*ap$coords[,1] + 20*ap$coords[,2]
##
## Comparing plots for data with and without trend
##
plot.geodata(ap)
plot.geodata(apt)
##
## Computing empirical variograms
##
## for data without trend
ap.v <- variog(ap, max.dist=1)
ap.v[1:4]
##
## for data with trend, ignoring the trend
apt.v0 <- variog(apt, max.dist=1)
apt.v0[1:4]
##
## for data with trend, with estimated 1st degree trend
apt.v1 <- variog(apt, trend="1st", max.dist=1)
apt.v1[1:4]
## Notice here that the following commands would produce exactly the same results:
## res1 <- lm(apt$data ~ apt$coords)$resid
## apt.v1.res1 <- variog(apt, data=res1, max.dist=1)
##
## for data with trend, with trend given by the covariate
apt.vc <- variog(apt, trend=~apt$covar, max.dist=1)
apt.vc[1:4]
## Notice here that the following commands would produce exactly the same results:
## res2 <- lm(apt$data ~ apt$covar)$resid
## apt.vc.res2 <- variog(apt, data=res2, max.dist=1)
##
## Ploting and coparing variograms
##
par(mfrow=c(2,2))
plot(ap.v)
plot(apt.v0)
plot(apt.v1)
plot(apt.vc)
##
