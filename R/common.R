"cov.spatial" <-
  function(obj, cov.model = c("matern", "exponential", "gaussian",
                  "spherical", "circular", "cubic", "wave",
                  "power", "powered.exponential", "cauchy",
                  "gneiting", "gneiting.matern", "pure.nugget"),
           cov.pars = stop("no cov.pars argument provided"),
           kappa = 0.5)
{
  ##
  ## checking/reading input
  ##
  cov.model <- match.arg(cov.model)
  if(cov.model == "matern" | cov.model == "powered.exponential" | 
     cov.model == "cauchy" | cov.model == "gneiting.matern"){
    if(is.null(kappa))
      stop("for matern, powered.exponential, cauchy and gneiting.matern covariance functions the parameter kappa must be provided")
    if(cov.model == "gneiting.matern" & length(kappa) != 2)
      stop("gneiting.matern correlation function model requires a vector with 2 parameters in the argument kappa")
    if((cov.model == "matern" | cov.model == "powered.exponential" | 
        cov.model == "cauchy") & length(kappa) != 1)
      stop("for this choice of  correlation function model kappa should be a scalar parameter")
    if(cov.model == "matern" & kappa == 0.5)
      cov.model == "exponential"
  }
  ##
  if(is.vector(cov.pars))
    sigmasq <- cov.pars[1]
  else sigmasq <- cov.pars[, 1]
  if(is.vector(cov.pars))
    phi <- cov.pars[2]
  else phi <- cov.pars[, 2]
  if(is.vector(cov.pars))
    ns <- 1
  else ns <- nrow(cov.pars)
  covs <- array(0, dim = dim(obj))
  ##
  if(cov.model == "power")
    if(any(phi >= 2) | any(phi <= 0))
      stop("for power model cov.pars[2] must be in the interval ]0,2[")
  ##
  ## computing correlations/covariances
  ##
  for(i in 1:ns) {
    if(phi[i] < 1e-12)
      cov.model <- "pure.nugget"
    cov.values <- switch(cov.model,
                         pure.nugget = rep(0, length(obj)),
                         wave = (1/obj) * (phi[i] * sin(obj/phi[i])),
                         exponential = exp( - (obj/phi[i])),
                         matern = matern(u = obj, phi = phi[i], kappa = kappa),
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
                         powered.exponential = {
                           if(kappa > 2 | kappa <= 0)
                             stop("for power exponential correlation model the parameter kappa must be in the intervel (0,2]"
                                  );
                           exp( - ((obj/phi[i])^kappa))
                         },
                         cauchy = (1 + (obj/phi[i])^2)^(-kappa),
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
  if(cov.model == "power") covs <- max(covs) - covs
  else covs[obj < 1e-15] <- sum(sigmasq)
  return(covs)
}

"cor.number" <- 
  function(cov.model= c("exponential", "matern", "gaussian",
             "spherical", "circular", "cubic", "wave", "power",
             "powered.exponential", "cauchy", "gneiting",
             "gneiting.matern", "pure.nugget"))
{
###	WARNING: codes above must be the same as in the C code
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

"read.geodata" <-
  function(file, header = FALSE, coords.col= 1:2, data.col = 3,
           data.names = NULL, covar.col = NULL, covar.names = "header",
           realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar"),
           rep.data.action, rep.covar.action, ...)
{
  call.fc <- match.call()
  ##
  obj <- read.table(file = file, header = header, ...)
  if(covar.names == "header"){
    if(!is.null(covar.col)){
      col.names <- names(obj)
      covar.names <- col.names[covar.col]
    }
    else covar.names <- NULL
  }
  ##
  if(missing(rep.data.action)) rep.data.action <- "none"
  if(!is.function(rep.data.action))
    rep.data.action <- match.arg(rep.data.action, choices = c("none", "first")) 
  if(missing(rep.covar.action)) rep.covar.action <- rep.data.action
  if(!is.function(rep.covar.action))
    rep.covar.action <- match.arg(rep.covar.action, choices = c("none", "first")) 
  ##
  res <- as.geodata(obj = obj, coords.col = coords.col, data.col = data.col,
                    covar.col = covar.col, covar.names = covar.names,
                    realisations = realisations, rep.data.action = rep.data.action,
                    rep.covar.action = rep.covar.action)
  res$call <- call.fc
  return(res)
}

"as.geodata" <-
  function(obj, coords.col = 1:2, data.col = 3, data.names = NULL, 
           covar.col = NULL, covar.names = "obj.names", realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar", "none"),
           rep.data.action, rep.covar.action)
{
  if(!is.matrix(obj) & !is.data.frame(obj))
    stop("object must be a matrix or data.frame")
  if(!is.null(data.names) & length(data.col) < 2)
    stop("data.names allowed only if there is more than 1 column of data")
  res <- list()
  ##
  ## testing for NA's setting the coordinates of the data locations
  ##
  if(any(is.na(obj[,coords.col]))){
    warning("NA's not allowed in the coordinates")
    obj <- obj[complete.cases(obj),]
    warning("eliminating rows with NA's")
  }
  res$coords <- as.matrix(obj[,coords.col])
  ##
  ## setting the data
  ##
  res$data <- as.matrix(obj[,data.col])
  if(length(data.col) == 1) res$data <- as.vector(res$data)
  else
    if(!is.null(data.names)) colnames(res$data) <- data.names
  ##
  ## setting the covariates, if the case 
  ##
  if(!is.null(covar.col)){
    res[[3]] <- as.matrix(obj[,covar.col])
    if(covar.names == "obj.names"){
      if(is.matrix(obj))      
        col.names <- dimnames(obj)[2]
      if(is.data.frame(obj))      
        col.names <- names(obj)
    }
    if(length(covar.col) == 1){
      if(covar.names == "obj.names"){
        if(is.null(col.names)) names(res)[3] <- "covariate"
        else names(res)[3] <- col.names[covar.col]
      }
      else
        names(res)[3] <- covar.names
    }
    else{
      names(res)[3] <- "covariate"
      if(covar.names == "obj.names")
        if(is.null(col.names)) colnames(res[[3]]) <- paste("covar", 1:length(covar.col), sep="")
        else  colnames(res[[3]]) <- col.names[covar.col]
      else
        colnames(res[[3]]) <- covar.names
    }
    res[[3]] <- as.data.frame(res[[3]])
  }
  ##
  ## Dealing with NA's
  ##
  if(na.action != "none"){
    if(na.action == "ifany")
      na.data <- na.covar <- TRUE
    if(na.action == "ifdata")
      {na.data <- TRUE; na.covar <- FALSE}
    if(na.action == "ifcovar")
      {na.data <- FALSE; na.covar <- TRUE}
    not.na <- function(x) !any(is.na(x))
    if(na.data){
      ind <- apply(as.matrix(res$data), 1, not.na)
      if(!all(ind)){
        res$coords <- res$coords[ind,]
        res$data <- drop(as.matrix(res$data)[ind,])
        if(!is.null(covar.col))
          res[[3]] <- drop(as.matrix(res[[3]][ind,]))
        cat(paste("as.geodata:", sum(!ind), "points removed due to NA in the data\n")) 
      }
    }
    if(!is.null(covar.col) && na.covar){
      ind <- apply(as.matrix(res[[3]]), 1, not.na)
      if(!all(ind)){
        res$coords <- res$coords[ind,]
        res$data <- drop(as.matrix(res$data)[ind,])
        if(!is.null(covar.col))
          res[[3]] <- drop(as.matrix(res[[3]][ind,]))
        cat(paste("as.geodata:", sum(!ind), "points removed due to NA in the covariate(s)\n")) 
      }
    }
  }
  ##
  ## Checking whether there are data from different realisations
  ##
  if(is.null(realisations)) realisations <- as.factor(rep(1, nrow(res$coords)))
  else{
    if(is.numeric(realisations) && length(realisations) == 1)
      realisations <- as.factor(obj[,realisations])
    res$realisations <- realisations
  }
  if(length(realisations) != nrow(res$coords))
    stop("realisations and coords have incompatible dimensions")
  ##
  ## Checking whether there are data at coincident locations
  ## and dealing with this acoording to the value of the argument
  ## rep.data.action 
  ##
  if(missing(rep.data.action)) rep.data.action <- "none"
  if(!is.function(rep.data.action))
    rep.data.action <- match.arg(rep.data.action, choices = c("none", "first")) 
  if(missing(rep.covar.action)) rep.covar.action <- rep.data.action
  if(!is.function(rep.covar.action))
    rep.covar.action <- match.arg(rep.covar.action, choices = c("none", "first")) 
  require(mva)
  if(is.function(rep.data.action) || rep.data.action == "first"){
    rep.lev <- as.character(paste("x",res$coords[,1],"y",res$coords[,2], sep=""))
    rep.dup <- duplicated(rep.lev)
    if(sum(rep.dup) > 0)
      cat(paste("as.geodata:", sum(rep.dup), "redundant locations found"))
    res$coords <- res$coords[!rep.dup,]
    measure.var.f <- function(x) return(summary(lm(x ~ as.factor(rep.lev)))$sigma^2)
    res$m.var <- drop(apply(as.matrix(res$data),2,measure.var.f))
    rep.action.f <- function(x, rep.action){ 
      if(!is.function(rep.action) && rep.action == "first")
        return(x[!rep.dup])
      else
        return((as.vector(by(x, rep.lev, rep.action))[codes(factor(rep.lev))])[!rep.dup])
    }
    res$data <- drop(apply(as.matrix(res$data), 2, rep.action.f, rep.action=rep.data.action))
    if(!is.null(covar.col))
      res[[3]] <- drop(apply(as.matrix(res[[3]]), 2, rep.action.f, rep.action=rep.covar.action))
    if(!is.null(res$realisations))
      res$realisations <- res$realisations[!rep.dup]
  }
  else{
    check.coincide <- function(x){sum(dist(x) < 1e-12) > 0}
    any.coincide <- lapply(split(as.data.frame(res$coords), realisations), check.coincide)
    any.coincide <- as.vector(unlist(any.coincide))
    if(sum(any.coincide) > 0)
      cat("WARNING: there are data at coincident locations, some of the geoR's functions will not work.\n")
  }
  ##
  class(res) <- "geodata"
  return(res)
}

"coords.aniso" <- 
  function(coords, aniso.pars, reverse=FALSE)
{
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if(length(aniso.pars) != 2)
    stop("argument aniso.pars must be a vector with 2 elementsm the anisotropy angle and anisotropy ratio, respectively")
  psiA <- aniso.pars[1]
  psiR <- aniso.pars[2]
  if(psiR < 1){
    psiR <- round(psiR, dig=8)
    if(psiR < 1)
      stop("anisotropy ratio must be greater than 1")
  }
  rm <- matrix(c(cos(psiA), -sin(psiA),
                 sin(psiA), cos(psiA)),
               ncol = 2)
  tm <- diag(c(1, 1/psiR))
  if(reverse)
    coords.mod <- coords %*% solve(rm %*% tm)
  else
    coords.mod <- coords %*% rm %*% tm
  return(coords.mod)
}

#"dist0.krige" <-
#function (x0, coords) 
#{
#  if (length(x0) != 2) 
#    stop(paste("length of x0 is", length(x0), "(it must be 2)"))
#  coords[, 1] <- coords[, 1] - x0[1]
#  coords[, 2] <- coords[, 2] - x0[2]
#  return(sqrt(coords[, 1]^2 + coords[, 2]^2))
#}


"polygrid" <- 
  function(xgrid, ygrid, poly, vec.inout = FALSE)
{
  if(is.R()){
    if(require(splancs) == FALSE){
      cat("ERROR: cannot run the function\n")
      cat("package \"splancs\" should be installed/loaded")
      return(invisible())
    }
  }
  else library(splancs)
  if(exists("inout")){
    xygrid <- expand.grid(x = xgrid, y = ygrid)
    ind <- as.vector(inout(pts=xygrid, poly=poly))
    xypoly <- xygrid[ind == TRUE,  ]
    if(vec.inout == FALSE)
      return(xypoly)
    else return(list(xypoly = xypoly, vec.inout = ind))
  }
  else{
    cat("ERROR: cannot run the function\n")
    cat("package \"splancs\" should be installed/loaded")
    return(invisible())
  }
}

#"variog.env" <-
#  function (x.variog, coords, model.pars, nsim = 99, messages.screen = TRUE)  
#{
#  cat("This function has been made obsolete\n")
#  cat("There are now two functions for variogram envelops:\n")
#  cat(" - variog.env.model:\n")
#  cat("       the same as the previous variog.env, based on the model")
#  cat(" - variog.env.mc:\n")
#  cat("       the new one based on permutations of the data")
#  return(invisible())
#}

"varcov.spatial" <-
  function (coords = NULL, dists.lowertri = NULL, cov.model = "matern",
            kappa = 0.5, nugget = 0, cov.pars = stop("no cov.pars argument"), 
            inv = FALSE, det = FALSE,
            func.inv = c("cholesky", "eigen", "svd", "solve"),
            scaled = FALSE, only.decomposition = FALSE, 
            sqrt.inv = FALSE, try.another.decomposition = TRUE,
            only.inv.lower.diag = FALSE) 
{
  if (is.R()) 
    require(mva)
  ##
  op.sem <- options()$show.error.message
  options(show.error.message = FALSE)
  on.exit(options(show.error.message = op.sem))
  ##
  func.inv <- match.arg(func.inv)
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (only.inv.lower.diag) 
    inv <- TRUE
  if (is.null(coords) & is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) & !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if (!is.null(coords)) 
    n <- nrow(coords)
  if (!is.null(dists.lowertri)) {
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  }
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
  if (!is.null(coords)) {
    dists.lowertri <- as.vector(dist(coords))
  }
  if (round(1e+12 * min(dists.lowertri)) == 0) 
    warning("Two or more pairs of data at coincident (or very close) locations. \nThis may cause crashes in some matrices operations.\n")
  varcov <- matrix(0, n, n)
  if (scaled) {
    if (all(phi < 1e-12)) 
      varcov <- diag(x = (1 + (tausq/sum(sigmasq))), n)
    else {
      if (is.vector(cov.pars)) 
        cov.pars.sc <- c(1, phi)
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
      varcov.sqrt <- try(chol(varcov))
      if (inherits(varcov.sqrt, "try-error")) {
        if (try.another.decomposition) 
          func.inv <- "eigen"
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
          if (det) 
            cov.logdeth <- sum(log(diag(varcov.sqrt)))
          if (sqrt.inv) 
            inverse.sqrt <- solve(varcov.sqrt)
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
      cov.logdeth <- try(sum(log(sqrt(varcov.svd$d))))
      if (inherits(cov.logdeth, "try-error")) {
        if (try.another.decomposition) 
          func.inv <- "eigen"
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
      varcov.eig <- try(eigen(varcov, symmetric = TRUE))
      cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))))
      if (inherits(cov.logdeth, "try.error") | inherits(varcov.eig, "try-error")) {
        diag(varcov) <- 1.0001 * diag(varcov)
        varcov.eig <- try(eigen(varcov, symmetric = TRUE))
        cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))))
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

"set.coords.lims" <-
  function(coords, xlim, ylim)
{
  coords.lims <- apply(coords, 2, range)
  if(!missing(xlim) && is.numeric(xlim)) coords.lims[,1] <- xlim[order(xlim)]
  if(!missing(ylim) && is.numeric(ylim)) coords.lims[,2] <- ylim[order(ylim)]
  coords.diff <- diff(coords.lims)
  if (coords.diff[1] != coords.diff[2]) {
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[, ind.min] <-
      coords.lims[, ind.min] +
        c(-coords.diff.diff, coords.diff.diff)/2
  }
  return(coords.lims)
}

