##
## Rvario.R: a tutorial on variograms
##
library(geoR)
cru <- read.geodata("Cruciani.dat", head=T, coords.col=2:3, data.col=4)
cru.borda <- read.table("Cruciani.border", head=T)[,2:3]
cru.borda <- rbind(cru.borda, cru.borda[1,])

cru.cloud <- variog(cru, option="cloud")
plot(cru.cloud)

args(variog)
#options(helphtml = TRUE)
#help.start()
#help(variog)

cru0.cloud <- variog(cru, option="cloud", lambda=0)
plot(cru0.cloud)

cru.cloud.m <- variog(cru, option="cloud", est="modulus")
plot(cru.cloud.m)

cru.cloud0.m <- variog(cru, option="cloud", lam=0, est="modulus")
plot(cru.cloud0.m)
##
##
##
par.ori <- par(no.readonly=TRUE)
par(mfrow=c(3,2), mar=c(3,3,0,0), mgp=c(1.5,.7,0))

cru.v1 <- variog(cru)
plot(cru.v1)
cru.v1

cru0.v1 <- variog(cru, lambda=0)
plot(cru0.v1)
cru0.v1

cru0.v2 <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
plot(cru0.v2)
cru0.v3 <- variog(cru, uvec=seq(1,15, l=8), lambda=0)
plot(cru0.v3)
cru0.v4 <- variog(cru, uvec=seq(0,9, l=12), lambda=0)
plot(cru0.v4)
cru0.v5 <- variog(cru, uvec=seq(1,15, l=12), lambda=0)
plot(cru0.v5)

par(par.ori)
##
##
##
par.ori <- par(no.readonly=TRUE)
par(mfrow=c(3,2), mar=c(3,3,0,0), mgp=c(1.5,.7,0))

cru.v1 <- variog(cru)
plot(cru.v1, ylim=c(0,3.5))
cru.v1

cru0.v1 <- variog(cru, lambda=0, max.dist=12)
plot(cru0.v1, ylim=c(0,3.5))
cru0.v1

cru0.v2 <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
plot(cru0.v2, ylim=c(0,3.5))
cru0.v3 <- variog(cru, uvec=seq(1,15, l=8), lambda=0)
plot(cru0.v3, ylim=c(0,3.5))
cru0.v4 <- variog(cru, uvec=seq(0,9, l=12), lambda=0)
plot(cru0.v4, ylim=c(0,3.5))
cru0.v5 <- variog(cru, uvec=seq(1,15, l=12), lambda=0)
plot(cru0.v5, ylim=c(0,3.5))

par(par.ori)
##
##
##
cru.v1 <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
plot(cru.v1, ylim=c(0,2.5))
lines.variomodel(list(nugget=0, cov.pars=c(2, 2), cov.model="exp"), max.dist=9)

lines.variomodel(list(nugget=0.5, cov.pars=c(2, 2), cov.model="exp"), max.dist=9, col="red")
lines.variomodel(list(nugget=0.5, cov.pars=c(2, 3), cov.model="exp"), max.dist=9, col="blue")
lines.variomodel(list(nugget=0.5, cov.pars=c(2, 2), cov.model="gau"), max.dist=9, col="green")

cru.v1 <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
plot(cru.v1)
cru.v1.exp <- variofit(cru.v1, ini=c(2,2), cov.model="exp")
lines(cru.v1.exp)

cru.v1.gau <- variofit(cru.v1, ini=c(2,2), cov.model="gau")
lines(cru.v1.gau, lty=2)

args(variofit)
#help(variofit)

cru.v1.exp1 <- variofit(cru.v1, ini=c(2,2), cov.model="exp", wei="eq")
lines(cru.v1.exp, col="blue")
cru.v1.gau1 <- variofit(cru.v1, ini=c(2,2), cov.model="gau", wei="eq")
lines(cru.v1.gau, lty=2, col="blue")

cru.v1 <- variog(cru, uvec=seq(0,9, l=8), lambda=0)
cru.env1 <- variog.mc.env(cru, obj=cru.v1)
plot(cru.v1, env=cru.env1)
##
##
##
data(wolfcamp)
#help(wolfcamp)

