require(geoR)
options(digits = 3, width = 80)
set.seed(21)
##
## Data from Isaaks & Srisvastava, 1989
##    pages 295-315
##
## Tests results for  Ordinary Kriging
##       functions: krige.conv, ksline, krige.bayes
##
"isaaks" <-
  structure(list(coords = matrix(c(61, 139, 63, 140, 64, 129, 68, 128, 71,
                   140,73, 141, 75, 128), ncol=2, byrow=T),
                 data=c(477, 696, 227, 646, 606, 791, 783),
                 x0=t(c(65, 137))
                 ),
            class = "geodata")
##
## True results
##
## Model: exponential, nugget=0, sigmasq=10, phi=10/3, aniso=c(0,1) 
Fig12.03 <- c(593, 8.96)
##
## Model: exponential, nugget=0, sigmasq=20, phi=10/3, aniso=c(0,1) 
Fig12.06 <- c(593, 17.91)
##
## Model: gaussian, nugget=0, sigmasq=10, phi=10/sqrt(3), aniso=c(0,1) 
Fig12.08 <- c(559, 4.78)
##
## Model: exponential, nugget=5, sigmasq=5, phi=10/3, aniso=c(0,1) 
Fig12.10 <- c(597, 10.31)
##
## Model: exponential, nugget=0, sigmasq=10, phi=20/3, aniso=c(0,1) 
Fig12.12 <- c(572, 5.76)
##
## Model: exponential, nugget=0, sigmasq=10, phi=5/3, aniso=c(3*pi/4,2) 
Fig12.16 <- c(620, 6.30)
##
## Model: exponential, nugget=0, sigmasq=10, phi=5/3, aniso=c(pi/4,2) 
Fig12.18 <- c(574, 8.16)
##
## Model: exponential, nugget=0, sigmasq=10, phi=5/3, aniso=c(3*pi/4,10) 
Fig12.20 <- c(654, 2.83)
##
isaaks.true <- rbind(Fig12.03, Fig12.06, Fig12.08, Fig12.10, Fig12.12, Fig12.16, Fig12.18, Fig12.20)
##
## Testing kriging.conv
##
out.kc <- output.control(mes = F)
##
##
kc12.03 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(10, 10/3)),
                             output = out.kc)[1:2])
##
kc12.06 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(20, 10/3)),
                             output = out.kc)[1:2])
##
kc12.08 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(10, 10/sqrt(3)),
                               cov.model="gaussian"),
                             output = out.kc)[1:2])
##
kc12.10 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(5, 10/3),
                               nugget=5),
                             output = out.kc)[1:2])
##
kc12.12 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(10, 20/3)),
                             output = out.kc)[1:2])
##
kc12.16 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(10, 10/3),
                               aniso.pars=c(3*pi/4, 2)),
                             output = out.kc)[1:2])
##
kc12.18 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(10, 10/3),
                               aniso.pars=c(pi/4, 2)),
                             output = out.kc)[1:2])
##
kc12.20 <- unlist(krige.conv(isaaks, loc=isaaks$x0,
                             krige=krige.control(cov.pars=c(10, 10/3),
                               aniso.pars=c(3*pi/4, 10)),
                             output = out.kc)[1:2])
##
isaaks.kc <- rbind(kc12.03, kc12.06, kc12.08, kc12.10, kc12.12, kc12.16, kc12.18, kc12.20)
isaaks.kc[,1] <- round(isaaks.kc[,1])
isaaks.kc[,2] <- round(isaaks.kc[,2], dig=2)
isaaks.kc
##
{
  if(all((isaaks.kc - isaaks.true) == 0)){
    cat("------------------\n")
    cat("krige.conv is ok\n")
    cat("------------------\n")
  }
  else{
    cat("----------------------------\n")
    cat("WARNING: krige.conv not ok\n")
    cat("----------------------------\n")
    print(isaaks.kc - isaaks.true)
  }
}
##
## Testing krige.bayes
##
kb12.03 <- krige.bayes(isaaks, loc=isaaks$x0, prior=prior.control(sigmasq.prior="fixed", sigmasq=10, phi.prior="fixed", phi=10/3),output=out.kc)
##
kb12.06 <- krige.bayes(isaaks, loc=isaaks$x0,prior=prior.control(sigmasq.prior="fixed", sigmasq=20, phi.prior="fixed", phi=10/3), output=out.kc)
##
kb12.08 <- krige.bayes(isaaks, loc=isaaks$x0,model=model.control(cov.model="gaussian"),prior=prior.control(sigmasq.prior="fixed", sigmasq=10, phi.prior="fixed", phi=10/sqrt(3)), output=out.kc)
##
kb12.10 <- krige.bayes(isaaks, loc=isaaks$x0,prior=prior.control(sigmasq.prior="fixed",sigmasq=5, tausq.rel=1,phi.prior="fixed", phi=10/3),output=output.control(mess=F, signal=F))
##
kb12.12 <- krige.bayes(isaaks, loc=isaaks$x0, prior=prior.control(sigmasq.prior="fixed", sigmasq=10, phi.prior="fixed", phi=20/3),output=out.kc)
##
kb12.16 <- krige.bayes(isaaks, loc=isaaks$x0,model=model.control(aniso.pars=c(3*pi/4, 2)),prior=prior.control(sigmasq.prior="fixed", sigmasq=10,phi.prior="fixed", phi=10/3),output=out.kc)
##
kb12.18 <- krige.bayes(isaaks, loc=isaaks$x0, model=model.control(aniso.pars=c(pi/4, 2)), prior=prior.control(sigmasq.prior="fixed", sigmasq=10, phi.prior="fixed", phi=10/3), output=out.kc)
##
kb12.20 <- krige.bayes(isaaks, loc=isaaks$x0,model=model.control(aniso.pars=c(3*pi/4, 10)),prior=prior.control(sigmasq.prior="fixed", sigmasq=10,phi.prior="fixed", phi=10/3), output=out.kc)
##
isaaks.kb <- rbind(unlist(kb12.03$pred[1:2]), unlist(kb12.06$pred[1:2]),
                   unlist(kb12.08$pred[1:2]), unlist(kb12.10$pred[1:2]),
                   unlist(kb12.12$pred[1:2]), unlist(kb12.16$pred[1:2]),
                   unlist(kb12.18$pred[1:2]), unlist(kb12.20$pred[1:2]))
isaaks.kb[,1] <- round(isaaks.kb[,1])
isaaks.kb[,2] <- round(isaaks.kb[,2], dig=2)
isaaks.kb
##
{
  if(all((isaaks.kb - isaaks.true) == 0)){
    cat("------------------\n")
    cat("krige.bayes is ok\n")
    cat("------------------\n")
  }
  else{
    cat("---------------------------\n")
    cat("WARNING: krige.bayes not ok\n")
    cat("---------------------------\n")
    print(isaaks.kb - isaaks.true)
  }
}
##
## Testing ksline
##
ks12.03 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(10, 10/3),
                         mess=F)[1:2])
ks12.06 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(20, 10/3),
                         mess=F)[1:2])
ks12.08 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(10, 10/sqrt(3)),
                         cov.model="gaussian", mess=F)[1:2])
ks12.10 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(5, 10/3),
                         nugget=5, mess=F)[1:2])
ks12.12 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(10, 20/3),
                         mess=F)[1:2])
ks12.16 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(10, 10/3),
                         aniso=c(3*pi/4, 2), mess=F)[1:2])
ks12.18 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(10, 10/3),
                         aniso=c(pi/4, 2), mess=F)[1:2])
ks12.20 <- unlist(ksline(isaaks, loc=isaaks$x0, cov.pars=c(10, 10/3),
                         aniso=c(3*pi/4, 10), mess=F)[1:2])
##
isaaks.ks <- rbind(ks12.03, ks12.06, ks12.08, ks12.10,
                   ks12.12, ks12.16, ks12.18, ks12.20)
##
isaaks.ks[,1] <- round(isaaks.ks[,1])
isaaks.ks[,2] <- round(isaaks.ks[,2], dig=2)
isaaks.ks
##
{
  if(all((isaaks.ks - isaaks.true) == 0)){
    cat("------------\n")
    cat("ksline is ok\n")
    cat("------------\n")
  }
  else{
    cat("----------------------\n")
    cat("WARNING: ksline not ok\n")
    cat("----------------------\n")
    print(isaaks.ks - isaaks.true)
  }
}
##
##q(save=FALSE)
#ap <- function(){
#  a <- 3 
#  b <- 5
#  print(ls(env=sys.frame(1)))
#  ap1 <- function()
#    {
#      cc <- 4
#      print(ls(env=sys.frame(1)))
#      print(ls(env=sys.frame(2)))
#      remove("a", envir=sys.frame(1))
#      print(ls(env=sys.frame(2)))
#      return(invisible())
#    }
#  ap1()
#  print(ls())
#  return(invisible())
#}  
#ap()
##
## Example 12.03, pages 290-296,  via matrix operations
##
V <- varcov.spatial(isaaks$coords, cov.pars=c(10, 10/3), nugget=0)
d0 <- loccoords(isaaks$coords, isaaks$x0)
v0 <- cov.spatial(d0, cov.pars=c(10, 10/3))
V <- V$varcov
V <- rbind(V,1)
V <- cbind(V,1)
V[8,8] <- 0
v0 <- c(v0, 1)
res <- c(round(solve(V,v0) %*% c(isaaks$data,0)),round(10-solve(V,v0) %*% v0, dig=2))
res
##
if(all((res-Fig12.03)==0)){
  cat("-----------------------------------\n")
  cat("matrices operations for Fig12.03 ok\n")
  cat("-----------------------------------\n")
}
if(any((res-Fig12.03)!=0)){
  cat("-----------------------------------\n")
  cat("WARNING: matrices operations NOT ok\n")
  cat("-----------------------------------\n")
  print(res)
}
##
## Example 12.10, pag 306,  via matrix operations
##
V <- varcov.spatial(isaaks$coords, cov.pars=c(5, 10/3), nugget=5)
d0 <- loccoords(isaaks$coords, isaaks$x0)
v0 <- cov.spatial(d0, cov.pars=c(5, 10/3))
V <- V$varcov
V <- rbind(V,1)
V <- cbind(V,1)
V[8,8] <- 0
v0 <- c(v0, 1)
res <- c(round(solve(V,v0) %*% c(isaaks$data,0)),round(10-solve(V,v0) %*% v0, dig=2))
res
##
if(all((res-Fig12.10)==0)){
  cat("-----------------------------------\n")
  cat("matrices operations for Fig12.10 ok\n")
  cat("------------------------------------\n")
}
if(any((res-Fig12.10)!=0)){
  cat("------------------------------------------------\n")
  cat("WARNING: matrices operations for Fig12.03 NOT ok\n")
  cat("------------------------------------------------\n")
  print(res)
}


