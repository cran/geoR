require(geoR)
set.seed(123)
options(digits = 3, width = 80)
data(s100)
ml <- likfit(s100, ini=c(0.5, 0.5), fix.nug = TRUE)
ml
summary(ml)
reml <- likfit(s100, ini=c(0.5, 0.5), fix.nug = TRUE, met = "REML")
summary(reml)
#plot(variog(s100))
#lines(ml)
#lines(reml, lty = 2)

ap <- grf(50, cov.pars=c(1, .3), nug=.3)
ml <- likfit(ap, ini=c(0.5, 0.5), nug=0.2)
ml
ml <- likfit(ap, ini=c(0.5, 0.5), fix.nug=TRUE, nug=0.2)
ml
ml <- likfit(ap, data=exp(ap$data), ini=c(0.5, 0.5), nug=0.2, fix.lambda=FALSE)
ml
ml <- likfit(ap, ini=c(0.5, 0.5), nug=0.2, fix.psiR = FALSE)
ml
ml <- likfit(ap, ini=c(0.5, 0.5), nug=0.2, fix.psiA = FALSE, fix.psiR = FALSE)
ml
ml <- likfit(ap, ini=c(0.5, 0.5), cov.model="matern", fix.kappa=TRUE, kappa=2)
ml
ml <- likfit(ap, ini=c(0.5, 0.5), cov.model="matern", fix.kappa=FALSE)
ml
ml <- likfit(ap, ini=expand.grid(c(.5,1,1.5),c(.1,.2,.3)))
ml
ml <- likfit(ap, ini=expand.grid(c(.5,1),c(.1,.2)), nug=c(.2,.3))
ml

## Multiple realizations
ap1 <- grf(40, cov.pars=c(1, .3))
ap1
ap2 <- grf(grid=ap1$coords[sample(1:40, 20),], cov.pars=c(1, .3))
ap2
ap3 <- grf(grid=ap1$coords[sample(1:40, 30),], cov.pars=c(1, .3))
ap3
ap <- list(coords = rbind(ap1$coords, ap2$coords, ap3$coords), data = c(ap1$data, ap2$data, ap3$data))
ap$realisations <- c(rep(1,40), rep(2,20), rep(3,30))
ap.fit <- likfit(ap, ini=c(.5, .5), reali=ap$real)
ap.fit


