require(geoR)
set.seed(123)
options(digits = 3, width = 80)

ex.data <- grf(50, cov.pars=c(10, .25))
ex.post <- krige.bayes(ex.data)
ex.post
ex1 <- krige.bayes(ex.data, prior = list(phi.prior = "fixed", phi = 0.3))
ex1 <- krige.bayes(ex.data, model = list(cov.model="spherical"))
##ex1 <- krige.bayes(ex.data, output = list(n.posterior = 100))
ex1 <- krige.bayes(ex.data, output = list(n.posterior = 10))
ex.grid <- as.matrix(expand.grid(seq(0,1,l=6), seq(0,1,l=6)))
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, prior =
                 prior.control(phi.discrete=seq(0, 2, l=3),
                 tausq.rel.discrete=seq(0, 2, l=3)),
                 output=output.control(n.post=100))
ex.bayes
plot(ex.data)
lines(ex.bayes, sum = mean) 
plot(ex.bayes)
lines(ex.bayes, summ="median", lty=2, post="par")
lines(ex.bayes, summ="mean", lwd=2, lty=2, post="par")
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(mgp = c(2,1,0))
image(ex.bayes, main="predicted values")
image(ex.bayes, val="variance", main="prediction variance")
image(ex.bayes, val= "simulation", number.col=1,
      main="a simulation from the \npredictive distribution")
image(ex.bayes, val= "simulation", number.col=2,
      main="another simulation from \nthe predictive distribution")

OC.test <- list(n.post = 10)

PC <- prior.control(phi.prior = c(.1, .2, .3, .2, ,.1, .1), phi.disc=seq(0.1, 0.6, l=6), tausq.rel.prior=c(.1, .4, .3, .2), tausq.rel.discrete=c(0,.1,.2,.3))
##ex.user <- krige.bayes(ex.data, prior = PC)  
ex.user <- krige.bayes(ex.data, prior = PC, output = OC.test)  
ex.user

# Simulating data at 2 different "times"
ap1 <- grf(50, cov.pars=c(1, .3))
ap2 <- grf(70, cov.pars=c(1, .3))
## A initial "usual" analysis
##ap1.kb <- krige.bayes(ap1)
ap1.kb <- krige.bayes(ap1, output = OC.test)
ap1.kb
## using the previous posterior as prior for next call
#ap2.kb <- krige.bayes(ap2, prior=post2prior(ap1.kb))
ap2.kb <- krige.bayes(ap2, prior=post2prior(ap1.kb), output = OC.test)
ap2.kb
##
## Another example with "user defined" prior
##
PC <- prior.control(phi.prior=c(.2,.3,.2,.1,.1,.1), phi.discrete = seq(0,.5,l=6))
##ap3.kb <- krige.bayes(ap1, prior = PC)
ap3.kb <- krige.bayes(ap1, prior = PC, output = OC.test)
ap3.kb
##
#ap4.kb <- krige.bayes(ap2, prior=post2prior(ap3.kb))
ap4.kb <- krige.bayes(ap2, prior=post2prior(ap3.kb), output = OC.test)
ap4.kb
##
## Now include tausq
##
PC <- prior.control(tausq.rel.prior = "uni", tausq.rel.discrete = seq(0, .5, l=6))
#ap5.kb <- krige.bayes(ap1)
ap5.kb <- krige.bayes(ap1, output = OC.test)
ap5.kb
## using the previous posterior as prior for next call
#ap6.kb <- krige.bayes(ap2, prior=post2prior(ap5.kb))
ap6.kb <- krige.bayes(ap2, prior=post2prior(ap5.kb), output = OC.test)
ap6.kb
##
##
##
PC <- prior.control(phi.prior=c(.2,.3,.2,.1,.1,.1), phi.discrete = seq(0,.5,l=6), tausq.rel.prior=c(.3,.4,.3), tausq.rel.discrete = c(0, .1, .2)) 
#ap7.kb <- krige.bayes(ap1, prior = PC)
ap7.kb <- krige.bayes(ap1, prior = PC, output = OC.test)
ap7.kb
#
#ap8.kb <- krige.bayes(ap2, prior=post2prior(ap7.kb))
ap8.kb <- krige.bayes(ap2, prior=post2prior(ap7.kb), output = OC.test)
ap8.kb

## with trend
data(s100)
prior2.b9 <- prior.control(beta.prior = "normal", beta = c(0,0,0),
beta.var = cbind(c(2,1.5,0),c(1.5,1.8,.2),c(0,0.2,1.5)),
phi.prior = "exponential", phi = 2.5, phi.discrete = c(2.5,3),
sigmasq.prior = "sc.inv.chisq", df.sigmasq = 5, sigmasq = 0.5)
#ap <- krige.bayes(s100,prior=prior2.b9, model=model.control(trend.d = "1st"))
ap <- krige.bayes(s100,prior=prior2.b9, model=model.control(trend.d = "1st"), output = OC.test)
ap

prior2.b9 <- prior.control(beta.prior = "normal", beta = c(0,0,0),
beta.var = cbind(c(2,1.5,0),c(1.5,1.8,0.5),c(0,0.5,1.5)),
phi.prior = "fixed", phi = 2.5,sigmasq.prior = "sc.inv.chisq",
df.sigmasq = 5, sigmasq = 0.5)
#ap <- krige.bayes(s100,prior=prior2.b9, model=model.control(trend.d="1st"))
ap <- krige.bayes(s100,prior=prior2.b9, model=model.control(trend.d="1st"), output = OC.test)
ap

prior2.b9 <- prior.control(beta.prior = "normal", beta = 0,beta.var =1,phi.prior = "fixed", phi = 2.5,sigmasq.prior ="sc.inv.chisq",df.sigmasq = 5, sigmasq = 0.5)
#ap <- krige.bayes(s100,prior=prior2.b9)
ap <- krige.bayes(s100,prior=prior2.b9, output = OC.test)
ap

par(op)
