require(geoR)
options(digits = 3, width = 80)
set.seed(50)
##
## Tests for coincident locations
##    functions:  krige.bayes, krige.conv and ksline  
##
## Simulating some data
##
ap <- grf(20, cov.pars=c(3, .25), mess=F)
ap$coords <- round(ap$coords, dig=3)
## defining prediction locations with 2 points coinciding with data locations
loci <- matrix(0, ncol=2, nrow=6)
loci[1,] <- c(.2,.3)
loci[2,] <- ap$coords[10,]
loci[3,] <- c(.5,.5)
loci[4,] <- c(.5,.8)
loci[5,] <- ap$coords[15,]
loci[6,] <- c(.7,.9)
##
coinc.data <- round(ap$data[c(10, 15)], dig=4) 
## just checking
ap1 <- loccoords(ap$coords, loci)
if(sum(ap1==0) !=2)
  stop("problems with example, more or less than 2 coincident points")
##
## testing krige.conv
##
out.kc <- output.control(mes=F)
##
##
cat("===================================================\n")
cat("TESTING KRIGE.CONV\n")
kc <- krige.conv(ap, loc=loci, krige=krige.control(cov.pars=c(1,.25)), output=out.kc)
lapply(kc[1:2], round, dig=4)
kc.pred <- round(kc$pred[c(2,5)], dig=4)
kc.var <- round(kc$krige.var[c(2,5)], dig=4)
##
{
  if(any((kc.pred-coinc.data) !=0)){
    cat("WARNING: predictions DO NOT coincide with data\n")
  }
  else{
    cat("predictions coincide with data, ok\n")
  }
}
{
  if(any(kc.var != 0)){
    cat("WARNING: kriging variance NOT equals to zero\n")
  }
  else{
    cat("kriging variance equals zero, ok\n")
  }
}
cat("===================================================\n")
##
## testing krige.bayes
##
cat("\n")
cat("\n")
cat("===================================================\n")
cat("TESTING KRIGE.BAYES\n")
cat("\n")
cat("----------------------\n")
cat("Testing with fixed phi\n")
cat("----------------------\n")
kb1 <- krige.bayes(ap, loc=loci, prior=prior.control(
                                   sigmasq.prior="fixed", sigmasq=1,
                                   phi.prior="fixed", phi=10/3),
                   output=output.control(mess=F, signal=F))
kb1$post
kb1$pred$mean
kb1$pred$var
kb1.pred <- round(kb1$pred$mean[c(2,5)], dig=4)
kb1.var <- round(kb1$pred$variance[c(2,5)], dig=4)
##
{
  if(any((kb1.pred-coinc.data) !=0)){
    cat("WARNING: predictions DO NOT coincide with data\n")}
  else{cat("predictions coincide with data, ok\n")}
}
{
  if(any(kb1.var != 0)){
    cat("WARNING: kriging variance NOT equals to zero\n")}
  else{cat("kriging variance equals zero, ok\n")}
}
cat("\n")
cat("-------------------------\n")
cat("Testing with variable phi\n")
cat("-------------------------\n")
kb2 <- krige.bayes(ap, loc=loci, prior=prior.control(phi.disc=seq(0,2,l=21)),output=output.control(mess=F, n.pred=10))
##
kb2.mom <- round(kb2$pred$mean[c(2,5)], dig=4)
kb2.mom.var <- round(kb2$pred$variance[c(2,5)], dig=4)
kb2.sim <- round(kb2$pred$simula[c(2,5),], dig=4)
kb2.mea <- round(kb2$pred$mean.sim[c(2,5)], dig=4)
kb2.var <- round(kb2$pred$variance.sim[c(2,5)], dig=4)
##
kb2.mom
kb2.mom.var
kb2.sim
kb2.mea
kb2.var
{
  if(any((kb2.mom-coinc.data) !=0)){
    cat("WARNING: moments DO NOT coincide with data\n")}
  else{cat("moments coincide with data, ok\n")}
}
{
  if(any(kb2.mom.var > 0)){
    cat("WARNING: moments variance DO NOT coincide with data\n")}
  else{cat("moments variance equals zero, ok\n")}
}
{
  if(any((kb2.sim-coinc.data) !=0)){
    cat("WARNING: simulations DO NOT coincide with data\n")}
  else{cat("simulations coincide with data, ok\n")}
}
{
  if(any((kb2.mea-coinc.data) !=0)){
    cat("WARNING: means DO NOT coincide with data\n")}
  else{cat("means coincide with data, ok\n")}
}
{
  if(any(kb2.var > 0)){
    cat("WARNING: variances DO NOT coincide with data\n")}
  else{cat("variances equal zero, ok\n")}
}
cat("===================================================\n")
