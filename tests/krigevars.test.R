require(geoR)
options(digits = 3, width = 80)
set.seed(30)
data(s100)
pred.grid <- expand.grid((1:4)/4,(1:4)/4)
##
##
pr0 <- prior.control(beta.prior="fixed", beta=0, sigmasq.prior = "fixed", sigmasq = 0.7323, phi.prior = "fixed", phi = 0.2)
test0 <- krige.bayes(s100, prior = pr0, locations = pred.grid)

pr1 <- prior.control(beta.prior="normal", beta=0, beta.var.std=0.01, sigmasq.prior = "fixed", sigmasq = 0.7323, phi.prior = "fixed", phi = 0.2)
test1 <- krige.bayes(s100, prior = pr1, locations = pred.grid)

pr2 <- prior.control(beta.prior="normal", beta=0, beta.var.std=10, sigmasq.prior = "fixed", sigmasq = 0.7323, phi.prior = "fixed", phi = 0.2)
test2 <- krige.bayes(s100, prior = pr2, locations = pred.grid)

pr3 <- prior.control(beta.prior="normal", beta=0, beta.var.std=100, sigmasq.prior = "fixed", sigmasq = 0.7323, phi.prior = "fixed", phi = 0.2)
test3 <- krige.bayes(s100, prior = pr3, locations = pred.grid)

pr4 <- prior.control(beta.prior="normal", beta=0, beta.var.std=10000, sigmasq.prior = "fixed", sigmasq = 0.7323, phi.prior = "fixed", phi = 0.2)
test4 <- krige.bayes(s100, prior = pr4, locations = pred.grid)

pr5 <- prior.control(beta.prior="flat", sigmasq.prior = "fixed", sigmasq = 0.7323, phi.prior = "fixed", phi = 0.2)
test5 <- krige.bayes(s100, prior = pr5, locations = pred.grid)

apm0 <- cbind(test0$pred[[1]],test1$pred[[1]],test2$pred[[1]],test3$pred[[1]],test4$pred[[1]],test5$pred[[1]])
apv0 <- cbind(test0$pred[[2]],test1$pred[[2]],test2$pred[[2]],test3$pred[[2]],test4$pred[[2]],test5$pred[[2]])
apm0
apv0
##
##
##
pr0 <- prior.control(beta.prior="fixed", beta=0, sigmasq.prior = "rec",  phi.prior = "fixed", phi = 0.2)
test0 <- krige.bayes(s100, prior = pr0, locations = pred.grid)

pr1 <- prior.control(beta.prior="normal", beta=0, beta.var.std=0.01, sigmasq.prior = "rec", phi.prior = "fixed", phi = 0.2)
test1 <- krige.bayes(s100, prior = pr1, locations = pred.grid)

pr2 <- prior.control(beta.prior="normal", beta=0, beta.var.std=10, sigmasq.prior = "rec", phi.prior = "fixed", phi = 0.2)
test2 <- krige.bayes(s100, prior = pr2, locations = pred.grid)

pr3 <- prior.control(beta.prior="normal", beta=0, beta.var.std=100, sigmasq.prior = "rec", phi.prior = "fixed", phi = 0.2)
test3 <- krige.bayes(s100, prior = pr3, locations = pred.grid)

pr4 <- prior.control(beta.prior="normal", beta=0, beta.var.std=10000, sigmasq.prior = "rec", phi.prior = "fixed", phi = 0.2)
test4 <- krige.bayes(s100, prior = pr4, locations = pred.grid)

pr5 <- prior.control(beta.prior="flat", sigmasq.prior = "rec", phi.prior = "fixed", phi = 0.2)
test5 <- krige.bayes(s100, prior = pr5, locations = pred.grid)

apm <- cbind(test0$pred[[1]],test1$pred[[1]],test2$pred[[1]],test3$pred[[1]],test4$pred[[1]],test5$pred[[1]])
apv <-  cbind(test0$pred[[2]],test1$pred[[2]],test2$pred[[2]],test3$pred[[2]],test4$pred[[2]],test5$pred[[2]])
apm
apv

apv[,6]/apv[,4]

