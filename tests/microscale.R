require(geoR)
options(digits = 3, width = 80)
set.seed(50)

op <- par(no.readonly=T)

dloc <- cbind((0:10)/10, 0)
set.seed(274)
dat <- grf(grid=dloc, cov.pars=c(1, .25), nug=.25)
ploc <- cbind((0:1000)/1000, 0)

par(mfrow=c(3,1), mar=c(3,3,0,0), mgp=c(1.5,.8,0))

##
## Predctions with option signal = F
##
## Here predictions "honor the data"
kp1 <- krige.conv(dat, loc=ploc, krige=krige.control(cov.pars=c(1, .25), nug=.25, micro=0), out=output.control(signal=F))
image(dat)
lines(ploc[,1], kp1$pred)
## idem here
kp2 <- krige.conv(dat, loc=ploc, krige=krige.control(cov.pars=c(1, .25), nug=.25, micro=0.25), out=output.control(signal=F))
image(dat)
lines(ploc[,1], kp2$pred)
## and here
kp3 <- krige.conv(dat, loc=ploc, krige=krige.control(cov.pars=c(1, .25), nug=.25, micro=0.125), out=output.control(signal=F))
image(dat)
lines(ploc[,1], kp3$pred)

##
## Now predictions with signal = T
##
## Here we predict the signal
kp4 <- krige.conv(dat, loc=ploc, krige=krige.control(cov.pars=c(1, .25), nug=.25, micro=0), out=output.control(signal=T))
image(dat)
lines(ploc[,1], kp4$pred)
## here we "honor the data" because microscale = nugget
kp5 <- krige.conv(dat, loc=ploc, krige=krige.control(cov.pars=c(1, .25), nug=.25, micro=0.25), out=output.control(signal=T))
image(dat)
lines(ploc[,1], kp5$pred)
## and here we are between the prevuious two
kp6 <- krige.conv(dat, loc=ploc, krige=krige.control(cov.pars=c(1, .25), nug=.25, micro=0.125), out=output.control(signal=T))
image(dat)
lines(ploc[,1], kp6$pred)

##
##
##
data(s100)
loci <-  rbind(s100$coords[1,],c(0.5,0.5))

kc <- krige.control(cov.pars=c(1,1), nugget = 1, micro.scale = 0)

oc <- output.control(signal=F, mess=F)
lapply(krige.conv(s100, locations =loci, krige= kc, output=oc)[1:2], round, dig=4)
oc <- output.control(signal=T, mess=F)
lapply(krige.conv(s100, locations =loci, krige= kc, output=oc)[1:2], round, dig=4)
##
##
##
kc <- krige.control(cov.pars=c(1,1), nugget = 1, micro.scale = 1)

oc <- output.control(signal=F, mess=F)
lapply(krige.conv(s100, locations =loci, krige= kc, output=oc)[1:2], round, dig=4)
oc <- output.control(signal=T, mess=F)
lapply(krige.conv(s100, locations =loci, krige= kc, output=oc)[1:2], round, dig=4)
##
##
kc <- krige.control(cov.pars=c(1,1), nugget = 1, micro.scale = 0.5)

oc <- output.control(signal=F, mess=F)
lapply(krige.conv(s100, locations =loci, krige= kc, output=oc)[1:2], round, dig=4)
oc <- output.control(signal=T, mess=F)
lapply(krige.conv(s100, locations =loci, krige= kc, output=oc)[1:2], round, dig=4)


