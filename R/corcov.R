##
## Correlations and covariances for the package geoR
## -------------------------------------------------
##
## Includes functions to compute cor. and cov,
## vectors and matrices and related operations
## 

"matern" <-
  function (u, phi, kappa) 
{
  if(is.vector(u)) names(u) <- NULL
  if(is.matrix(u)) dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0,
                 (((2^(-(kappa-1)))/gamma(kappa)) *
                  (uphi^kappa) *
                  besselK(x=uphi, nu=kappa)), 1)    
  uphi[u > 600*phi] <- 0 
  return(uphi)
}

"cor.number" <- 
  function(cov.model= c("exponential", "matern", "gaussian",
             "spherical", "circular", "linear", "cubic", "wave", "power",
             "powered.exponential", "cauchy", "gneiting",
             "gneiting.matern", "pure.nugget"))
{
###	WARNING: codes below MUST be the same as in the C code
###              "cor_diag"
  cov.model <- match.arg(cov.model)
  cornumber <- switch(cov.model,
                      pure.nugget = as.integer(1),
                      exponential = as.integer(2),
                      spherical = as.integer(3),
                      gaussian = as.integer(4),
                      wave = as.integer(5),
                      cubic = as.integer(6),
                      power = as.integer(7),
                      powered.exponential = as.integer(8),
                      cauchy = as.integer(9),
                      gneiting = as.integer(10),
                      circular = as.integer(11),
                      matern = as.integer(12),
                      gneiting.matern = as.integer(13),
                      stop("wrong or no specification of cov.model")
                      )
  return(cornumber)
}

"cov.spatial" <-
  function(obj, cov.model = 'matern', cov.pars = stop("no cov.pars argument provided"),
           kappa = 0.5)
{
  ## extracting covariance paramters
  if(is.vector(cov.pars)) sigmasq <- cov.pars[1]
  else sigmasq <- cov.pars[, 1]
  if(is.vector(cov.pars)) phi <- cov.pars[2]
  else phi <- cov.pars[, 2]
  if(is.null(kappa)) kappa <- NA
  ## checking for nested models 
  if(is.vector(cov.pars)) ns <- 1
  else{
    ns <- nrow(cov.pars)
    if(length(cov.model) == 1) cov.model <- rep(cov.model, ns)
    if(length(kappa) == 1) kappa <- rep(kappa, ns)
  }
  if(length(cov.model) != ns) stop('wrong length for cov.model')
  if(length(kappa) != ns) stop('wrong length for kappa')
  ##
  cov.model <- sapply(cov.model, match.arg, c("matern", "exponential", "gaussian",
                                              "spherical", "circular", "cubic", "wave",
                                              "linear", "power", "powered.exponential", "cauchy",
                                              "gneiting", "gneiting.matern", "pure.nugget"))
  ## settings for power model (do not reverse order of the next two lines!)
  phi[cov.model == "linear"] <- 1
  cov.model[cov.model == "linear"] <- "power"
  ## checking input for cov. models with extra parameter(s)
  if(any(cov.model == 'gneiting.model') && ns > 1)
    stop('nested models including the gneiting.matern are not implemented') 
    for(i in 1:ns){
    if(cov.model[i] == "matern" | cov.model[i] == "powered.exponential" | 
       cov.model[i] == "cauchy" | cov.model[i] == "gneiting.matern"){
      if(is.na(kappa[i]))
        stop("for matern, powered.exponential, cauchy and gneiting.matern covariance functions the parameter kappa must be provided")
      if(cov.model[i] == "gneiting.matern" & length(kappa) != 2*ns)
        stop("gneiting.matern correlation function model requires a vector with 2 parameters in the argument kappa")
      if((cov.model[i] == "matern" | cov.model[i] == "powered.exponential" | 
          cov.model[i] == "cauchy") & length(kappa) != 1*ns)
        stop("kappa must have 1 parameter for this correlation function")
      if(cov.model[i] == "matern" & kappa[i] == 0.5) cov.model[i] == "exponential"
    }
    if(cov.model[i] == "power")
      if(any(phi[i] >= 2) | any(phi[i] <= 0))
        stop("for power model the phi parameters must be in the interval ]0,2[")
  }
  ##
  ## computing correlations/covariances
  ##
  covs <- array(0, dim = dim(obj))
  for(i in 1:ns) {
    if(phi[i] < 1e-12)
      cov.model[i] <- "pure.nugget"
    cov.values <- switch(cov.model[i],
                         pure.nugget = rep(0, length(obj)),
                         wave = (1/obj) * (phi[i] * sin(obj/phi[i])),
                         exponential = exp( - (obj/phi[i])),
                         matern = matern(u = obj, phi = phi[i], kappa = kappa[i]),
                         gaussian = exp( - ((obj/phi[i])^2)),
                         spherical = ifelse(obj < phi[i], (1 - 1.5 * (obj/phi[i]) +
                           0.5 * (obj/phi[i])^3), 0),
                         circular = {
                           obj.sc <- obj/phi[i];
                           obj.sc[obj.sc > 1] <- 1;
                           ifelse(obj < phi[i], (1 - (2 * ((obj.sc) *
                                                           sqrt(1 - ((obj.sc)^2)) +
                                                           asin(obj.sc)))/pi), 0)
                         },
                         cubic = {
                           obj.sc <- obj/phi[i];
                           ifelse(obj < phi[i], (1 - (7 * (obj.sc^2) -
                                                      8.75 * (obj.sc^3) +
                                                      3.5 * (obj.sc^5) -
                                                      0.75 * (obj.sc^7))), 0)
                         },
                         power = (obj)^phi,
                         powered.exponential = exp( - ((obj/phi[i])^kappa[i])),
                         cauchy = (1 + (obj/phi[i])^2)^(-kappa[i]),
                         gneiting = {
                           obj.sc <- obj/phi[i];
                           t2 <- (1 - obj.sc);
                           t2 <- ifelse(t2 > 0, (t2^8), 0);
                           (1 + 8 * obj.sc + 25 * (obj.sc^2) + 32 * (obj.sc^
                                                                     3)) * t2
                         },
                         gneiting.matern = { 
                           obj.sc <- obj/(phi[i] * kappa[2]);
                           t2 <- (1 - obj.sc);
                           t2 <- ifelse(t2 > 0, (t2^8), 0);
                    cov.values <- (1 + 8 * obj.sc + 25 * (obj.sc^2) + 32 * (obj.sc^3)) * t2;
                    cov.values * matern(u = obj, phi = phi[i], kappa = kappa[1])
                  },
                  stop("wrong or no specification of cov.model")
                  )
    cov.values <- sigmasq[i] * cov.values
    covs <- covs + cov.values
  }
  if(all(cov.model == "power")) covs <- max(covs) - covs
  else covs[obj < 1e-15] <- sum(sigmasq)
  return(covs)
}

"varcov.spatial" <-
  function(coords = NULL, dists.lowertri = NULL, cov.model = "matern",
            kappa = 0.5, nugget = 0, cov.pars = stop("no cov.pars argument"), 
            inv = FALSE, det = FALSE,
            func.inv = c("cholesky", "eigen", "svd", "solve"),
            scaled = FALSE, only.decomposition = FALSE, 
            sqrt.inv = FALSE, try.another.decomposition = TRUE,
            only.inv.lower.diag = FALSE) 
{
  if (is.R()) require(mva)
  ##
##  op.sem <- options()$show.error.message
##  options(show.error.message = FALSE)
##  on.exit(options(show.error.message = op.sem))
  require(methods)
  if(!exists("trySilent")){
    error.now <- options()$show.error.messages
    if (is.null(error.now) | error.now) 
      on.exit(options(show.error.messages = TRUE))
    options(show.error.messages = FALSE)
  }
  ##
  func.inv <- match.arg(func.inv)
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "linear", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (only.inv.lower.diag)  inv <- TRUE
  if (is.null(coords) & is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) & !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if (!is.null(coords))  n <- nrow(coords)
  if (!is.null(dists.lowertri))
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
  }
##  print(c(tausq=tausq, sigmasq=sigmasq, phi=phi, kappa=kappa))
  if (!is.null(coords)) dists.lowertri <- as.vector(dist(coords))
  if (round(1e+12 * min(dists.lowertri)) == 0) 
    warning("Two or more pairs of data at coincident (or very close) locations. \nThis may cause crashes in some matrices operations.\n")
  varcov <- matrix(0, n, n)
  if (scaled) {
    if (all(phi < 1e-12)) 
      varcov <- diag(x = (1 + (tausq/sum(sigmasq))), n)
    else {
      if (is.vector(cov.pars)) cov.pars.sc <- c(1, phi)
      else cov.pars.sc <- cbind(1, phi)
      covvec <- cov.spatial(obj = dists.lowertri, cov.model = cov.model, 
                            kappa = kappa, cov.pars = cov.pars.sc)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      if (is.R()) remove("covvec")
      else remove("covvec", frame = sys.nframe())
      diag(varcov) <- 1 + (tausq/sum(sigmasq))
    }
  }
  else {
    if (all(sigmasq < 1e-10) | all(phi < 1e-10)) {
      varcov <- diag(x = (tausq + sum(sigmasq)), n)
    }
    else {
      covvec <- cov.spatial(obj = dists.lowertri, cov.model = cov.model, 
                            kappa = kappa, cov.pars = cov.pars)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      if (is.R()) remove("covvec")
      else remove("covvec", frame = sys.nframe())
      diag(varcov) <- tausq + sum(sigmasq)
    }
  }
  if (inv | det | only.decomposition | sqrt.inv | only.inv.lower.diag) {
    if (func.inv == "cholesky") {
      if(exists("trySilent"))
        varcov.sqrt <- trySilent(chol(varcov))
      else
        varcov.sqrt <- try(chol(varcov))
      if (inherits(varcov.sqrt, "try-error")) {
        if (try.another.decomposition){
          cat("trying another decomposition (svd)\n")
          func.inv <- "svd"
        }
        else {
          print(varcov.sqrt[1])
          stop()
        }
      }
      else {
        if (only.decomposition | inv) 
          if (is.R()) remove("varcov")
          else remove("varcov", frame = sys.nframe())
        if (only.decomposition == FALSE) {
          if (det) cov.logdeth <- sum(log(diag(varcov.sqrt)))
          if (sqrt.inv) inverse.sqrt <- solve(varcov.sqrt)
          if (inv) {
            if (is.R()) {
              invcov <- chol2inv(varcov.sqrt)
              if (!sqrt.inv)
                remove("varcov.sqrt")
            }
            else {
              invcov.sqrt <- solve.upper(varcov.sqrt)
              invcov <- invcov.sqrt %*% t(invcov.sqrt)
              if (!sqrt.inv) 
                remove("varcov.sqrt", frame = sys.nframe())
            }
          }
        }
      }
    }
    if (func.inv == "svd") {
      varcov.svd <- svd(varcov, nv = 0)
      if(exists("trySilent"))
        cov.logdeth <- trySilent(sum(log(sqrt(varcov.svd$d))))
      else
        cov.logdeth <- try(sum(log(sqrt(varcov.svd$d))))
      if (inherits(cov.logdeth, "try-error")) {
        if (try.another.decomposition){
          cat("trying another decomposition (eigen)\n")
          func.inv <- "eigen"
        }
        else {
          print(cov.logdeth[1])
          stop()
        }
      }
      else {
        if (only.decomposition | inv) 
          if (is.R())  remove("varcov")
          else remove("varcov", frame = sys.nframe())
        if (only.decomposition) 
          varcov.sqrt <- t(varcov.svd$u %*% (t(varcov.svd$u) * 
                                             sqrt(varcov.svd$d)))
        if (inv) {
          invcov <- t(varcov.svd$u %*% (t(varcov.svd$u) * 
                                        (1/varcov.svd$d)))
        }
        if (sqrt.inv) 
          inverse.sqrt <- t(varcov.svd$u %*% (t(varcov.svd$u) * 
                                              (1/sqrt(varcov.svd$d))))
      }
    }
    if (func.inv == "solve") {
      if (det) 
        stop("the option func.inv == \"solve\" does not allow computation of determinants. \nUse func.inv = \"chol\",\"svd\" or \"eigen\"\n")
      if(exists("trySilent"))
        invcov <- trySilent(solve(varcov))
      else
        invcov <- try(solve(varcov))
      if (inherits(cov.logdeth, "try-error")) {
        if (try.another.decomposition) 
          func.inv <- "eigen"
        else {
          print(invcov[1])
          stop()
        }
      }
      if (is.R()) remove("varcov")
      else remove("varcov", frame = sys.nframe())
    }
    if (func.inv == "eigen") {
      if(exists("trySilent")){
        varcov.eig <- trySilent(eigen(varcov, symmetric = TRUE))
        cov.logdeth <- trySilent(sum(log(sqrt(varcov.eig$val))))
      }
      else{
        varcov.eig <- try(eigen(varcov, symmetric = TRUE))
        cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))))
      }
      if (inherits(cov.logdeth, "try.error") | inherits(varcov.eig, "try-error")) {
        diag(varcov) <- 1.0001 * diag(varcov)
        if(exists("trySilent")){
          varcov.eig <- trySilent(eigen(varcov, symmetric = TRUE))
          cov.logdeth <- trySilent(sum(log(sqrt(varcov.eig$val))))
        }
        else{
          varcov.eig <- try(eigen(varcov, symmetric = TRUE))
          cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))))
        }
        if (inherits(cov.logdeth, "try.error") | inherits(varcov.eig, "try-error")) {
          return(list(crash.parms = c(tausq=tausq, sigmasq=sigmasq, phi=phi, kappa=kappa)))
        }
      }
      else {
        if (only.decomposition | inv) 
          if (is.R()) remove("varcov")
          else remove("varcov", frame = sys.nframe())
        if (only.decomposition) 
          varcov.sqrt <- (varcov.eig$vec %*% diag(sqrt(varcov.eig$val)) %*% 
                          t(varcov.eig$vec))
        if (inv) 
          invcov <- (varcov.eig$vec %*% diag(1/varcov.eig$val) %*% 
                     t(varcov.eig$vec))
        if (sqrt.inv) 
          inverse.sqrt <- (varcov.eig$vec %*% diag(1/sqrt(varcov.eig$val)) %*% 
                           t(varcov.eig$vec))
      }
    }
  }
  if (only.decomposition == FALSE) {
    if (det) {
      if (inv) {
        if (only.inv.lower.diag) 
          result <- list(lower.inverse = invcov[lower.tri(invcov)], 
                         diag.inverse = diag(invcov), log.det.to.half = cov.logdeth)
        else result <- list(inverse = invcov, log.det.to.half = cov.logdeth)
      }
      else {
        result <- list(varcov = varcov, log.det.to.half = cov.logdeth)
      }
      if (sqrt.inv) 
        result$sqrt.inverse <- inverse.sqrt
    }
    else {
      if (inv) {
        if (only.inv.lower.diag) 
          result <- list(lower.inverse = invcov[lower.tri(invcov)], 
                         diag.inverse = diag(invcov))
        else {
          if (sqrt.inv) 
            result <- list(inverse = invcov, sqrt.inverse = inverse.sqrt)
          else result <- list(inverse = invcov)
        }
      }
      else result <- list(varcov = varcov)
    }
  }
  else result <- list(sqrt.varcov = varcov.sqrt)
  result$crash.parms <- NULL
  return(result)
}


