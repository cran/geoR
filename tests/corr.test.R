require(geoR)
options(digits = 3, width = 80)
tc <- function(x, phi, kappa, cov.model)
{
  res <- rep(0.0, length(x))
  .C("veccorrval",
     as.double(phi),
     as.double(kappa),
     as.double(x),
     as.integer(length(x)),
     cornr = as.integer(cor.number(cov.model)),
     out = as.double(res))$out
}

x <- 1:100
for(cm in c("exponential", "matern", "gaussian",
    "spherical", "circular", "cubic", "wave", "power", "powered.exponential",
    "cauchy", "gneiting", "pure.nugget"))
  {
    kp <- 1.9
    x1 <- tc(x, 30, kp, cm)
    x2 <- cov.spatial(x, cm, cov.pars=c(1, 30), kappa = kp)
    if(all(x1-x2) < 1e-8)
      cat(paste(cm, "OK\n"))
    else{
      cat(cm)
      cat("\n")
    }
    print(x1)
    print(x2)
  }

