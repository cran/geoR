require(geoR)
options(digits = 3, width = 80)
set.seed(18)
##
## Examples and test for the function likfit
## using the options for realisations
##
##
## Simulating 3 data-set at the same locations
##
ap1 <- grf(40, cov.pars=c(3, .25))
ap2 <- grf(grid=ap1$coords, cov.pars=c(3, .25))
ap3 <- grf(grid=ap1$coords, cov.pars=c(3, .25))
##
## Pooling together the 3 simulations
##
sim01 <- list(coords=rbind(ap1$coords,ap2$coords,ap3$coords), data=c(ap1$data,ap2$data,ap3$data),realisation=rep(1:3,rep(40,3)))
##
## estimating parameters of for each relisation
##
ap1.ml <- likfit(ap1, ini=c(3,1))
ap2.ml <- likfit(ap2, ini=c(3,1))
ap3.ml <- likfit(ap3, ini=c(3,1))
ap1.ml
ap2.ml
ap3.ml
##
## estimating parameters for the pooled data
##
sim01.ml <- likfit(sim01, ini=c(3,1), reali = sim01$real)
##
## Changing the order of the coordinates
##
ind1 <- sample(1:40)
ind2 <- sample(1:40)
ind3 <- sample(1:40)
sim02 <- list(coords=rbind(ap1$coords[ind1,],ap2$coords[ind2,],ap3$coords[ind3,]), data=c(ap1$data[ind1],ap2$data[ind2],ap3$data[ind3]),realisation=rep(1:3,rep(40,3)))
sim02.ml <- likfit(sim02, ini=c(3,1), reali = sim02$real)
sim01.ml
sim02.ml
##
## Now a simulation where the different realisations have
## partially different locations 
##
ind1 <- sample(1:40, 30)
ind2 <- sample(1:40, 20)
ind3 <- sample(1:40, 25)
sim03 <- list(coords=rbind(ap1$coords[ind1,],ap2$coords[ind2,],ap3$coords[ind3,]), data=c(ap1$data[ind1],ap2$data[ind2],ap3$data[ind3]),realisation=rep(1:3,c(30, 20, 25)))
sim03.ml <- likfit(sim03, ini=c(3,1), reali = sim03$real)
sim03.ml
##
