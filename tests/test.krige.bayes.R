require(geoR)
options(digits = 3, width = 80)
set.seed(25)
## Simulating data at 2 different "times"
ap1 <- grf(50, cov.pars=c(1, .3))
ap2 <- grf(70, cov.pars=c(1, .3))
## A initial "usual" analysis
ap1.kb <- krige.bayes(ap1)
ap1.kb$post
## using the previous posterior as prior for next call
ap2.kb <- krige.bayes(ap2, prior=post2prior(ap1.kb))
ap2.kb$post
##
## Another example with "user defined" prior
##
PC <- prior.control(phi.prior=c(.2,.3,.2,.1,.1,.1), phi.discrete = seq(0,.5,l=6))
ap3.kb <- krige.bayes(ap1, prior = PC)
ap3.kb$post
##
ap4.kb <- krige.bayes(ap2, prior=post2prior(ap3.kb))
ap4.kb$post
##
## Now include tausq
##
PC <- prior.control(tausq.rel.prior = "uni", tausq.rel.discrete = seq(0, .5, l=6))
ap5.kb <- krige.bayes(ap1)
ap5.kb$post
## using the previous posterior as prior for next call
ap6.kb <- krige.bayes(ap2, prior=post2prior(ap5.kb))
ap6.kb$post
##
##
##
PC <- prior.control(phi.prior=c(.2,.3,.2,.1,.1,.1), phi.discrete = seq(0,.5,l=6), tausq.rel.prior=c(.3,.4,.3), tausq.rel.discrete = c(0, .1, .2)) 
ap7.kb <- krige.bayes(ap1, prior = PC)
ap7.kb$post
#
ap8.kb <- krige.bayes(ap2, prior=post2prior(ap7.kb))
ap8.kb$post

