require(geoR)
options(digits = 3, width = 80)
set.seed(30)
##rm(list=ls())
##
##
locs <- expand.grid((1:3)/3,(1:3)/3)
##
##
output.sk <- output.ok <- output.control(n.pred=1000, mes=F)
#output.sk <- output.ok <- output.control(mes=F)
output.kb <- output.control(n.pred=1000, signal=F, mess=F)
##
## Simple kriging
##
cat("==============\n")
cat("Simple kriging\n")
cat("==============\n")
##
##
ap <- grf(70, cov.pars=c(1, .3), mes=F)
ap.ml <- likfit(ap, ini=c(1, .3), mes=F)
ap.ml
summary(ap.ml)
krige.sk <- krige.control(type="sk", obj.model=ap.ml)
model.kb <- model.control()
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.sk01 <- krige.conv(ap, loc=locs, krige=krige.sk, output=output.sk)
ap.kb01 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output=output.kb)
ap.sk01[1:2]
ap.kb01$pred[1:2]
##
dm <- round(ap.sk01[[1]] - ap.kb01$pred[[1]], dig=8)
dv <- round(ap.sk01[[2]] - ap.kb01$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm1 sk mean not ok\n"), "dm1 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv1 sk var not ok\n"), "dv1 sk var OK\n"))
cat("-----------------------------\n")
##
##
ap$data <- exp(ap$data)
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0, mes=F)
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml)
model.kb <- model.control(lambda = 0)
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.sk02 <- krige.conv(ap, loc=locs, krige= krige.sk, output=output.sk)
ap.kb02 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.sk02[1:2]
ap.kb02$pred[1:2]
##
dm <- round(ap.sk02[[1]] - ap.kb02$pred[[1]], dig=8)
dv <- round(ap.sk02[[2]] - ap.kb02$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm2 sk mean not ok\n"), "dm2 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv2 sk var not ok\n"), "dv2 sk var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.5, mes=F)
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml)
model.kb <- model.control(lambda = 0.5)
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.sk03 <- krige.conv(ap, loc=locs, krige=krige.sk, output=output.sk)
ap.kb03 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.sk03[1:2]
ap.kb03$pred[1:2]
##
dm <- round(ap.sk03[[1]] - ap.kb03$pred[[1]], dig=8)
dv <- round(ap.sk03[[2]] - ap.kb03$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm3 sk mean not ok\n"), "dm3 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv3 sk var not ok\n"), "dv3 sk var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.25, mes=F)
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml)
model.kb <- model.control(lambda = 0.25)
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
seed.now <- .Random.seed
ap.sk04 <- krige.conv(ap, loc=locs, krige=krige.control(type="sk", obj.model=ap.ml), output=output.sk)
.Random.seed <- seed.now
ap.kb04 <- krige.bayes(ap, loc=locs, model = model.kb, prior=prior.kb, output=output.kb)
ap.sk04[1:2]
ap.kb04$pred[1:2]
##
dm <- round(ap.sk04[[1]] - ap.kb04$pred[[1]], dig=8)
dv <- round(ap.sk04[[2]] - ap.kb04$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm4 sk mean not ok\n"), "dm4 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv4 sk var not ok\n"), "dv4 sk var OK\n"))
cat("-----------------------------\n")
##
## Ordinary kriging
##
cat("================\n")
cat("Ordinary kriging\n")
cat("================\n")
##
##
ap$data <- log(ap$data)
##
##
ap.ml <- likfit(ap, ini=c(1, .3), mes=F)
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml)
model.kb <- model.control()
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.ok01 <- krige.conv(ap, loc=locs, krige=krige.ok, output=output.ok)
ap.kb01 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output=output.kb)
ap.ok01[1:2]
ap.kb01$pred[1:2]
##
dm <- round(ap.ok01[[1]] - ap.kb01$pred[[1]], dig=8)
dv <- round(ap.ok01[[2]] - ap.kb01$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm1 ok mean not ok\n"), "dm1 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv1 ok var not ok\n"), "dv1 ok var OK\n"))
cat("-----------------------------\n")
##
##
ap$data <- exp(ap$data)
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0, mes=F)
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml)
model.kb <- model.control(lambda = 0)
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.ok02 <- krige.conv(ap, loc=locs, krige= krige.ok, output=output.ok)
ap.kb02 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.ok02[1:2]
ap.kb02$pred[1:2]
##
dm <- round(ap.ok02[[1]] - ap.kb02$pred[[1]], dig=8)
dv <- round(ap.ok02[[2]] - ap.kb02$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm2 ok mean not ok\n"), "dm2 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv2 ok var not ok\n"), "dv2 ok var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.5, mes=F)
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml)
model.kb <- model.control(lambda = 0.5)
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.ok03 <- krige.conv(ap, loc=locs, krige=krige.ok, output=output.ok)
ap.kb03 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.ok03[1:2]
ap.kb03$pred[1:2]
##
dm <- round(ap.ok03[[1]] - ap.kb03$pred[[1]], dig=8)
dv <- round(ap.ok03[[2]] - ap.kb03$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm3 ok mean not ok\n"), "dm3 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv3 ok var not ok\n"), "dv3 ok var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.25, mes=F)
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml)
model.kb <- model.control(lambda = 0.25)
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
seed.now <- .Random.seed
ap.ok04 <- krige.conv(ap, loc=locs, krige=krige.control(type="ok", obj.model=ap.ml), output=output.ok)
.Random.seed <- seed.now
ap.kb04 <- krige.bayes(ap, loc=locs, model = model.kb, prior=prior.kb, output=output.kb)
ap.ok04[1:2]
ap.kb04$pred[1:2]
##
dm <- round(ap.ok04[[1]] - ap.kb04$pred[[1]], dig=8)
dv <- round(ap.ok04[[2]] - ap.kb04$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm4 ok mean not ok\n"), "dm4 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv4 ok var not ok\n"), "dv4 ok var OK\n"))
cat("-----------------------------\n")
##
##
## With trend model
##
##
cat("=========================\n")
cat("Simple kriging with trend\n")
cat("=========================\n")
##
##
ap <- grf(70, cov.pars=c(1, .3), mes=F)
ap$data <- ap$data - 2 * ap$coords[,1] + 3 * ap$coords[,2] + rnorm(70, sd=.5)
##
##
model.kb <- model.control(trend.d = "1st", trend.l = "1st")
##
## Simple kriging
##
ap.ml <- likfit(ap, ini=c(1, .3), mes=F, trend = "1st")
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml, trend.d = "1st", trend.l = "1st")
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.sk01 <- krige.conv(ap, loc=locs, krige=krige.sk, output=output.sk)
ap.kb01 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output=output.kb)
ap.sk01[1:2]
ap.kb01$pred[1:2]
##
dm <- round(ap.sk01[[1]] - ap.kb01$pred[[1]], dig=8)
dv <- round(ap.sk01[[2]] - ap.kb01$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm1 sk mean not ok\n"), "dm1 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv1 sk var not ok\n"), "dv1 sk var OK\n"))
cat("-----------------------------\n")
##
##
ap$data <- exp(ap$data)
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0, mes=F, trend = "1st")
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml, trend.d = "1st", trend.l = "1st")
model.kb$lambda <- 0
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.sk02 <- krige.conv(ap, loc=locs, krige= krige.sk, output=output.sk)
ap.kb02 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.sk02[1:2]
ap.kb02$pred[1:2]
##
dm <- round(ap.sk02[[1]] - ap.kb02$pred[[1]], dig=8)
dv <- round(ap.sk02[[2]] - ap.kb02$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm2 sk mean not ok\n"), "dm2 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv2 sk var not ok\n"), "dv2 sk var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.5, mes=F, trend = "1st")
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml, trend.d = "1st", trend.l = "1st")
model.kb$lambda <- 0.5
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.sk03 <- krige.conv(ap, loc=locs, krige=krige.sk, output=output.sk)
ap.kb03 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.sk03[1:2]
ap.kb03$pred[1:2]
##
dm <- round(ap.sk03[[1]] - ap.kb03$pred[[1]], dig=8)
dv <- round(ap.sk03[[2]] - ap.kb03$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm3 sk mean not ok\n"), "dm3 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv3 sk var not ok\n"), "dv3 sk var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.25, mes=F, trend = "1st")
ap.ml
krige.sk <- krige.control(type="sk", obj.model=ap.ml, trend.d = "1st", trend.l = "1st")
model.kb$lambda <- 0.25
prior.kb <- prior.control(beta.prior="fixed",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], beta=ap.ml$beta, sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
seed.now <- .Random.seed
ap.sk04 <- krige.conv(ap, loc=locs, krige=krige.sk, output=output.sk)
.Random.seed <- seed.now
ap.kb04 <- krige.bayes(ap, loc=locs, model = model.kb, prior=prior.kb, output=output.kb)
ap.sk04[1:2]
ap.kb04$pred[1:2]
##
dm <- round(ap.sk04[[1]] - ap.kb04$pred[[1]], dig=8)
dv <- round(ap.sk04[[2]] - ap.kb04$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm4 sk mean not ok\n"), "dm4 sk mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv4 sk var not ok\n"), "dv4 sk var OK\n"))
cat("-----------------------------\n")
##
##
##
cat("=============================\n")
##
## Ordinary kriging
##
cat("===========================\n")
cat("Ordinary kriging with trend\n")
cat("===========================\n")
##
ap$data <- log(ap$data)
##
##
model.kb <- model.control(trend.d = "1st", trend.l = "1st")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), mes=F, trend = "1st")
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml, trend.d = "1st", trend.l = "1st")
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.ok01 <- krige.conv(ap, loc=locs, krige=krige.ok, output=output.ok)
ap.kb01 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output=output.kb)
ap.ok01[1:2]
ap.kb01$pred[1:2]
##
dm <- round(ap.ok01[[1]] - ap.kb01$pred[[1]], dig=8)
dv <- round(ap.ok01[[2]] - ap.kb01$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm1 ok mean not ok\n"), "dm1 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv1 ok var not ok\n"), "dv1 ok var OK\n"))
cat("-----------------------------\n")
##
##
ap$data <- exp(ap$data)
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0, mes=F, trend = "1st")
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml,trend.d = "1st", trend.l = "1st")
model.kb$lambda  <- 0
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.ok02 <- krige.conv(ap, loc=locs, krige= krige.ok, output=output.ok)
ap.kb02 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.ok02[1:2]
ap.kb02$pred[1:2]
##
dm <- round(ap.ok02[[1]] - ap.kb02$pred[[1]], dig=8)
dv <- round(ap.ok02[[2]] - ap.kb02$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm2 ok mean not ok\n"), "dm2 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv2 ok var not ok\n"), "dv2 ok var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.5, mes=F, trend = "1st")
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml,trend.d = "1st", trend.l = "1st")
model.kb$lambda <- 0.5
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
ap.ok03 <- krige.conv(ap, loc=locs, krige=krige.ok, output=output.ok)
ap.kb03 <- krige.bayes(ap, loc=locs, prior=prior.kb, model = model.kb, output = output.kb)
ap.ok03[1:2]
ap.kb03$pred[1:2]
##
dm <- round(ap.ok03[[1]] - ap.kb03$pred[[1]], dig=8)
dv <- round(ap.ok03[[2]] - ap.kb03$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm3 ok mean not ok\n"), "dm3 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv3 ok var not ok\n"), "dv3 ok var OK\n"))
cat("-----------------------------\n")
##
##
ap.ml <- likfit(ap, ini=c(1, .3), lambda=0.25, mes=F, trend = "1st")
ap.ml
krige.ok <- krige.control(type="ok", obj.model=ap.ml, trend.d = "1st", trend.l = "1st")
model.kb$lambda <- 0.25
prior.kb <- prior.control(beta.prior="flat",sigmasq.prior="fixed",phi.prior="fixed", tausq.rel = ap.ml$nug/ap.ml$cov.pars[1], sigmasq=ap.ml$cov.pars[1], phi=ap.ml$cov.pars[2])
##
seed.now <- .Random.seed
ap.ok04 <- krige.conv(ap, loc=locs, krige=krige.ok, output=output.ok)
.Random.seed <- seed.now
ap.kb04 <- krige.bayes(ap, loc=locs, model = model.kb, prior=prior.kb, output=output.kb)
ap.ok04[1:2]
ap.kb04$pred[1:2]
##
dm <- round(ap.ok04[[1]] - ap.kb04$pred[[1]], dig=8)
dv <- round(ap.ok04[[2]] - ap.kb04$pred[[2]], dig=8)
cat("-----------------------------\n")
cat(ifelse(any(dm) > 0, stop("dm4 ok mean not ok\n"), "dm4 ok mean OK\n"))
cat(ifelse(any(dv) > 0, stop("dv4 ok var not ok\n"), "dv4 ok var OK\n"))
cat("-----------------------------\n")
##
##
##
## Other options 
##
output.kb <- output.control(n.pred=500, mess=F, quant=T, thres=c(.5, 1.5))
##
## fixed phi
##
cat("=====================================\n")
cat("examples with fixed phi\n")
cat("-------------------------------------\n")
prior.kb <- prior.control(phi.prior="fixed", phi=0.3)
ap.kb05 <- krige.bayes(ap, loc=locs, prior=prior.kb, output = output.kb)
cat("example 1 run ok\n")
ap.kb05$pred[1:2]
##
prior.kb <- prior.control(phi.prior="fixed", phi=0.3, beta.prior="normal", beta=-1, beta.var=0.2, sigmasq.prior="sc.inv.chisq", sigmasq=1.5, df.sig=12)
ap.kb06 <- krige.bayes(ap, loc=locs, prior=prior.kb, output = output.kb)
cat("example 2 run ok\n")
ap.kb06$pred[1:2]
##
cat("=====================================\n")
##
## All random
##
cat("\n", "\n")
cat("=====================================\n")
cat("examples with fixed phi and/or tausq.rel\n")
cat("-------------------------------------\n")
prior.kb <- prior.control(phi.discrete=seq(0,1,l=11), tausq.rel.discrete= seq(0,.2, l=5))
ap.kb07 <- krige.bayes(ap, loc=locs, prior=prior.kb, output = output.kb)
cat("example 1 run ok\n")
ap.kb07$pred[1:2]
##
prior.kb <- prior.control(phi.discrete=seq(0,1,l=11))
ap.kb08 <- krige.bayes(ap, loc=locs, prior=prior.kb, output = output.kb)
cat("example 2 run ok\n")
ap.kb08$pred[1:2]
##
prior.kb <- prior.control(phi.discrete=seq(0,1,l=11), tausq.rel.discrete= seq(0,.2, l=5), beta.prior="normal", beta=-1, beta.var=0.2, sigmasq.prior="sc.inv.chisq", sigmasq=1.5, df.sig=12)
ap.kb09 <- krige.bayes(ap, loc=locs, prior=prior.kb, output = output.kb)
cat("example 3 run ok\n")
ap.kb09$pred[1:2]
##
prior.kb <- prior.control(phi.discrete=seq(0,1,l=11), beta.prior="normal", beta=-1, beta.var=0.2, sigmasq.prior="sc.inv.chisq", sigmasq=1.5, df.sig=12)
ap.kb10 <- krige.bayes(ap, loc=locs, prior=prior.kb, output = output.kb)
cat("example 4 run ok\n")
ap.kb10$pred[1:2]
##
cat("=====================================\n")


