"grf" <-
  function(n, grid = "irreg",
           nx = round(sqrt(n)), ny = round(sqrt(n)),
           xlims = c(0, 1), ylims = c(0, 1), nsim = 1, 
           cov.model = "matern",
           cov.pars = stop("covariance parameters (sigmasq and phi) needed"),
           kappa = 0.5,  nugget=0, lambda=1, aniso.pars = NULL,
           method = c("cholesky", "svd", "eigen", "circular.embedding"),
           messages.screen = TRUE)
{
  ##
  ## reading and checking input
  ##
  call.fc <- match.call()
  method <- match.arg(method)
  if((method == "circular.embedding") & messages.screen)
    cat("grf: for simulation of fiends with large number of points the consider the packages RandomFields should be considered.\n") 
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(!is.null(kappa))
    if(cov.model == "matern" & kappa == 0.5) cov.model <- "exponential"
  rseed <- .Random.seed
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    nst <- 1
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
    nst <- nrow(cov.pars)
  }
  sill.total <- tausq + sum(sigmasq)
  messa <- grf.aux1(nst, nugget, sigmasq, phi, kappa, cov.model)
  if (messages.screen) {
    cat(messa$nst)
    cat(messa$nugget)
    cat(messa$cov.structures)
    cat(paste("grf: decomposition algorithm used is: ", method, "\n"))
  }
  results <- list()
  ##
  ## defining the locations for the simulated data
  ##
  if (is.matrix(grid) | is.data.frame(grid)) {
    results$coords <- as.matrix(grid)
    if (messages.screen) 
      cat("grf: simulation(s) on a grid provided by the user\n")
  }
  else {
    if (grid == "irreg") {
      results$coords <- cbind(x = runif(n, xlims[1], xlims[2]),
                              y = runif(n, ylims[1], ylims[2]))
      if (messages.screen) 
        cat(paste("grf: simulation(s) on random locations with ", n, " points\n"))
    }
    else {
      results$coords <- as.matrix(expand.grid(x = seq(xlims[1], xlims[2], l = nx),
                                              y = seq(ylims[1], ylims[2], l = ny)))
      if (messages.screen) 
        cat(paste("grf: generating grid ", nx, " * ", ny, 
                  " with ", (nx*ny), " points\n"))
    }
  }
  n <- nrow(results$coords)
  ##
  ## transforming to the isotropic space 
  ##
  if(!is.null(aniso.pars)) {
    if(method == "circular.embedding")
      stop("anisotropic models not implemented for the circular embedding method. \nConsider using the package \"RandomFields")
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if(messages.screen)
      cat("grf: transforming to the isotropic space \n")
    results$coords <- coords.aniso(coords = results$coords,
                                   aniso.pars = aniso.pars)
  }
  ##
  ## simulating data at locations defined by the matrix results$coords
  ##
  if (all(phi) == 0) {
    results$data <- matrix(rnorm((n * nsim), mean = 0, sd = sqrt(sill.total)), 
                           nrow = n, ncol = nsim)
  }
  else {
    if (method == "circular.embedding") {
      if (grid == "irreg") 
        stop("Option for \"circular.embedding\" algorithm only allowed for regular grids. You might have to include the argument grid=\"reg\"")
      stepx <- (xlims[2] - xlims[1])/(nx - 1)
      stepy <- (ylims[2] - ylims[1])/(ny - 1)
      if (round(1e+08 * stepx) != round(1e+08 *stepy)) 
        stop("grf: distance between grid points must be the same in X and Y directions")
      temp <- list(n = n, nst = nst, sigmasq = sigmasq, 
                   xlims = xlims, ylims = ylims, stepx = stepx, 
                   cov.model = cov.model, phi = phi, kappa = kappa)
      if(messages.screen)
        cat("\ngrf: WARNING:\nmessages of the type mtot=XXXXX will appear on your screen. \nIf there are many (3 or more, say) or they run indefinitely, you should stop the simulation and try again with a different grid (e.g. try to add 1 point in each direction)\n")
      grf.aux3 <- function (nsim, temp) {
        realiz <- rep(0, temp$n)
        for (i in 1:temp$nst) {
          realiz <- realiz + sqrt(temp$sigmasq[i]) *
            grf.aux2(xlim = temp$xlims, 
                     ylim = temp$ylims, step = temp$stepx,
                     cov.model = temp$cov.model,
                     phi = temp$phi[i], kappa = temp$kappa)
          NULL
        }
        return(realiz)
      }      
      results$data <- apply(as.matrix(1:nsim), 1, grf.aux3, temp = temp)
      if (nugget != 0) {
        results$data <- results$data + matrix(rnorm((n * nsim), sd = sqrt(nugget)), 
                                              ncol = nsim)
      }
    }
    else{
      results$data <- matrix(rnorm((n * nsim)), nrow = n, ncol = nsim)
      cov.decomp <- t(varcov.spatial(coords = results$coords, 
                                   cov.model = cov.model, kappa = kappa,
                                   nugget = nugget, cov.pars = cov.pars, 
                                   only.decomposition = TRUE,
                                   func.inv = method)$sqrt.varcov)
      results$data <- cov.decomp %*% results$data
    }
    if (nsim == 1) 
      results$data <- as.vector(results$data)
  }
  ##
  ## transforming data (Box - Cox)
  ##
  if (lambda != 1){
    if (lambda != 0)
      results$data <- (results$data * lambda + 1)^(1/lambda)
    else
      results$data <- exp(results$data)
    messa$transformation <- paste("grf: Data transformed (Box-Cox), for lambda =", lambda)
    if (messages.screen) 
      cat(messa$transformation); cat("\n")
  }
  ##
  ## back-transformation to the anisotropic space 
  ##
  if(!is.null(aniso.pars)) {
    if(messages.screen)
      cat("grf: back-transforming to the anisotropic space \n")
    results$coords <- coords.aniso(coords = results$coords,
                                   aniso.pars = aniso.pars, reverse=TRUE)
  }
  else{aniso.pars <- "no anisotropy parameters provided/used"}
  ##
  ## preparing output
  ##
  if (messages.screen) 
    cat(paste("grf: End of simulation procedure. Number of realizations:",
              nsim, "\n"))
  results  <- c(results, list(cov.model = cov.model, 
                              nugget = nugget, cov.pars = cov.pars,
                              kappa = kappa, lambda = lambda,
                              aniso.pars = aniso.pars, method = method,
                              .Random.seed = rseed, messages = messa,
                              call = call.fc))
  class(results) <- c("grf", "geodata")
  return(results)
}

"grf.aux1" <-
  function (nst, nugget, sigmasq, phi, kappa, cov.model) 
{
  cov.nst <- paste("grf: process with ", nst, " covariance structure(s)\n")
  cov.nugget <- paste("grf: nugget effect is: tausq=", nugget,"\n")
  cov.message <- NULL
  for (i in 1:nst) {
    if (phi[i] == 0) 
      cov.message[i] <- paste("grf: covariance model", i, "is a pure nugget effect\n")
    else {
      if (cov.model == "matern" | cov.model == "powered.exponential" | 
          cov.model == "cauchy" | cov.model == "gneiting-matern") 
        cov.message[i] <- paste("grf: covariance model ", 
                                i, " is: ", cov.model, "(sigmasq=", sigmasq[i], 
                                ", phi=", phi[i], ", kappa = ", kappa, ")\n", sep = "")
      else cov.message[i] <- paste("grf: covariance model ", 
                                   i, " is: ", cov.model, "(sigmasq=", sigmasq[i], 
                                   ", phi=", phi[i], ")\n", sep = "")
    }
  }
  return(list(nst = cov.nst, nugget = cov.nugget, cov.structures = cov.message))
}
"grf.aux2" <-
  function (xlim, ylim, step, cov.model, phi, kappa = 0.5) 
{
  if(!is.null(kappa))
    if (cov.model == "matern" & kappa == 0.5) {
      cov.model <- "exponential"
    }
  cs <- switch(cov.model,
               spherical = c(1,1),
               gneiting.matern = c(2, kappa),
               powered.exponential = c(3,kappa),
               exponential = c(3,1),
               gaussian = c(3,2),
               matern = c(4, kappa),
               gneiting = c(5,1),
               cauchy = c(6,kappa),
               wave = c(7,1))
  covfct <- cs[1]
  shape <- cs[2]
  parameter <- c(phi, shape, 1)
  nx <- c(diff(xlim)/step + 1, diff(ylim)/step + 1)
  storage.mode(nx) <- "integer"
  ln <- 2
  storage.mode(ln) <- "integer"
  covnr <- covfct
  storage.mode(covnr) <- "integer"
  simustep <- step
  storage.mode(simustep) <- "double"
  param <- parameter
  storage.mode(param) <- "double"
  lparam <- length(parameter)
  storage.mode(lparam) <- "integer"
  res <- double(nx[1] * nx[2])
  mm <- integer(ln)
  x <- .C("woodandchan", covnr, nx, ln, simustep, param, lparam, 
          res = res, m = mm)$res
  cat("\n")
  return(x)
}


"likfit.nospatial" <-
  function(temp.list, ...)
{
  results <- list()
  z <- temp.list$z
  n <- temp.list$n
  beta.size <- temp.list$beta.size
  xmat <- temp.list$xmat
  txmat <- temp.list$txmat
  ixx <- solve(crossprod(xmat))
  if(temp.list$fix.lambda == FALSE){
    if (temp.list$minimisation.function == "nlm"){
      assign(".temp.lower.lambda",-2, pos=1)
      assign(".temp.upper.lambda", 2, pos=1)
      results <- nlm(proflik.lambda, 1, ...)
      if(exists(".temp.lambda")){
        results$lambda <- .temp.lambda
        remove(".temp.lambda", pos=1, inherits = TRUE)
      }
      else{
        results$lambda <- results$estimate
      }
      rm(.temp.lower.lambda, .temp.upper.lambda, inherits = TRUE, pos=1)
    }
    if (temp.list$minimisation.function == "nlmP"){
      results <- nlmP(proflik.lambda, 1, lower=-2, upper=2,...)  
      results$lambda <- results$estimate
    }
    if (temp.list$minimisation.function == "optim"){
      results <- optim(1, proflik.lambda, method="L-BFGS-B", lower=-2, upper=2,...)
      results$minimum <- results$value
      results$lambda <- results$par
    }
    if(results$lambda == 1) {
      temp.list$log.jacobian <- 0
    }
    else {
      if(any(z <= 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      if(any(z^(results$lambda - 1) <= 0))
        temp.list$log.jacobian <- log(prod(z^(results$lambda - 1)))
      else temp.list$log.jacobian <- sum(log(z^(results$lambda - 1)))
      if(results$lambda == 0)
        z <- log(z)
      else z <- ((z^results$lambda) - 1)/results$lambda
    }
  }
  else{
    results$lambda <- temp.list$lambda
    results$code <- 1
    if (temp.list$minimisation.function == "optim") results$convergence <- 0
  }
  ssres <- (z %*% (diag(n) - xmat %*%
                   solve(crossprod(xmat)) %*% txmat) %*% z)
  if(temp.list$method == "ML"){
    results$tausqhat <- ssres/n
    if(temp.list$fix.lambda)
      results$minimum <- as.vector(((n/2) * log(2 * pi) +
                          (n/2) * log(results$tausqhat) +
                          (n/2)  -
                          temp.list$log.jacobian))
  }
  if(temp.list$method == "RML") {
    results$tausqhat  <- (ssres/(n-beta.size))
    if(temp.list$fix.lambda){
      results$minimum <- as.vector((((n - beta.size)/2) * log(2 * pi) +
                          ((n - beta.size)/2) * log(results$tausqhat) +
                          (n/2) -
                          temp.list$log.jacobian
                          ))
    }
  }
  if (temp.list$minimisation.function == "optim") results$value <- results$minimum    
  return(results)
}

"loglik.spatial" <-
function(pars)
{
  tausq <- pars[1]
  sigmasq <- pars[2]
  sill.total <- tausq + sigmasq
  phi <- pars[3]
  lambda <- pars[4]
  z <- .temp.list$z
  n <- .temp.list$n
  if(.temp.list$fix.lambda == FALSE) {
    if(lambda == 1) {
      .temp.list$log.jacobian <- 0
    }
    else {
      if(any(z < 0))
        stop("Transformation option not allowed when there are zeros or negative data"
             )
      if(any(z^(lambda - 1) <= 0))
        .temp.list$log.jacobian <- log(prod(z^(lambda -
                                               1)))
      else .temp.list$log.jacobian <- sum(log(z^(lambda - 1)))
      if(lambda == 0)
        z <- log(z)
      else z <- ((z^lambda) - 1)/lambda
    }
  }
  beta.size <- .temp.list$beta.size
  kappa <- .temp.list$kappa
  covinf <- varcov.spatial(dists.lowertri = .temp.list$
                           dists.lowertri, cov.model = .temp.list$cov.model,
                           kappa = kappa, nugget = tausq,
                           cov.pars = c(sigmasq, phi), scaled = FALSE,
                           inv = TRUE, det = TRUE,
                           only.inv.lower.diag = TRUE)
  xix <- as.double(rep(0, beta.size*beta.size))
  xix <- .C("bilinearform_XAY",
            as.double(covinf$lower.inverse),
            as.double(covinf$diag.inverse),
            as.double(as.vector(.temp.list$xmat)),
            as.double(as.vector(.temp.list$xmat)),
            as.integer(beta.size),
            as.integer(beta.size),
            as.integer(n),
            res = xix)$res
  attr(xix, "dim") <- c(beta.size, beta.size)
  if(length(as.vector(xix)) == 1) {
    ixix <- 1/xix
    choldet <- 0.5 * log(xix)
  }
  else {
    chol.xix <- chol(xix)
    ixix <- chol2inv(chol.xix)
    choldet <- sum(log(diag(chol.xix)))
  }
  xiy <- as.double(rep(0, beta.size))
  xiy <- .C("bilinearform_XAY",
            as.double(covinf$lower.inverse),
            as.double(covinf$diag.inverse),
            as.double(as.vector(.temp.list$xmat)),
            as.double(as.vector(z)),
            as.integer(beta.size),
            as.integer(1),
            as.integer(n),
            res = xiy)$res
  beta.hat <- as.vector(ixix %*% xiy)
  yiy <- as.double(0.0)
  yiy <- .C("bilinearform_XAY",
            as.double(covinf$lower.inverse),
            as.double(covinf$diag.inverse),
            as.double(as.vector(z)),
            as.double(as.vector(z)),
            as.integer(1),
            as.integer(1),
            as.integer(n),
            res = yiy)$res
  ssresmat <- as.vector(yiy - 2*crossprod(beta.hat,xiy) +  beta.hat %*% xix %*% beta.hat)
  if(.temp.list$method == "ML") {
    loglik <- ( - (n/2) * log(2 * pi) -
               covinf$log.det.to.half -
               0.5 * ssresmat + 
               .temp.list$log.jacobian)
  }
  if(.temp.list$method == "RML") {
    xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE, only.values = TRUE)
    loglik <- ( - ((n - beta.size)/2) * log(2 * pi) +
               0.5 * sum(log(xx.eigen$values)) -
               covinf$log.det.to.half -
               (0.5) * ssresmat -
               choldet +
               .temp.list$log.jacobian)
  }
  return(as.vector(loglik))
}

"matern" <-
  function (u, phi, kappa) 
{
  if(is.vector(u)) names(u) <- NULL
  if(is.matrix(u)) dimnames(u) <- list(NULL, NULL)
  uvec <- u/phi
  ucov <- ifelse(uvec > 0, ((((2^(kappa-1))*gamma(kappa))^(-1))*(uvec^kappa)*besselK(x=uvec, nu=kappa)), 1)    
  ## The following needs more checking:
  ##        min.non.inf <- min(ucov[ucov != Inf])
  ##        dist.min.non.inf <- min(u[ucov <= 1.01 * min.non.inf])
  ##        ucov[u <= dist.min.non.inf] <- 0
  ##        print(c(phi, dist.min.non.inf, min.non.inf))        
  ##
  ucov[u > 600*phi] <- 0 
  return(ucov)
}


"nlmP" <-
  function(objfunc, params, lower = rep( -Inf, length(params)),
           upper = rep(+Inf, length(params)), ... )
{
  ## minimizer, using nlm with transformation of variables
  ##
  ## objfunc is a function to be optimised
  ## params is a starting value for the parameters
  Nparams <- length(params)
  if(length(lower) != Nparams)
    stop(" lower boundry different length than params")
  if(length(upper) != Nparams)
    stop(" upper boundry different length than params")
  checklimits <- upper - lower
  if(any(checklimits <= 0))
    stop(" bad boundries")
  if(any(params < lower))
    stop(" starting params too low")
  if(any(params > upper))
    stop(" starting params too high")
  
  bothlimQQ <- (lower != (-Inf)) & (upper != +Inf)
  loweronlyQQ <- (lower != (-Inf)) & (upper == +Inf)
  upperonlyQQ <- (lower == (-Inf)) & (upper != +Inf)
  ubothQQ <- upper[bothlimQQ]
  lbothQQ <- lower[bothlimQQ]
  dbothQQ <- ubothQQ - lbothQQ
  loneQQ <- lower[loweronlyQQ]
  uoneQQ <- upper[upperonlyQQ]
  
  .bounds.list <- list(bothlimQQ = bothlimQQ, 
                       loweronlyQQ = loweronlyQQ,
                       upperonlyQQ = upperonlyQQ,
                       ubothQQ = ubothQQ,
                       lbothQQ = lbothQQ,
                       dbothQQ = dbothQQ,
                       loneQQ = loneQQ,
                       uoneQQ = uoneQQ)
  
  assign(".objfuncQQ", objfunc, pos=1)
  assign(".bounds.list", .bounds.list, pos=1)
  
  ## reduce the parameter space by a scale to keep parameters
  ## away from the boundries
  
  normaltomad <- function(normalparamsX)
    {
      madparamsX <- normalparamsX
      if(any(.bounds.list$bothlimQQ)) {
        noughtone <- (normalparamsX[.bounds.list$bothlimQQ] -
                      .bounds.list$lbothQQ)/.bounds.list$dbothQQ
        madparamsX[.bounds.list$bothlimQQ] <- log(noughtone/(1 - noughtone))
      }
      
      if(any(.bounds.list$loweronlyQQ))
        madparamsX[.bounds.list$loweronlyQQ] <-
          log(normalparamsX[.bounds.list$loweronlyQQ] - .bounds.list$loneQQ)
      
      if(any(.bounds.list$upperonlyQQ))
        madparamsX[.bounds.list$upperonlyQQ] <-
          log(.bounds.list$uoneQQ - normalparamsX[.bounds.list$upperonlyQQ])
      
      return(madparamsX)
    }
  
  madtonormalQQ <<- function(madparamsX)
    {
      normalparamsX <- madparamsX
      
      if(any(.bounds.list$bothlimQQ)) {
###        madparamsX[((.bounds.list$bothlimQQ) & (madparamsX > 300))] <- 300
        emad <- exp(madparamsX[.bounds.list$bothlimQQ])
        normalparamsX[.bounds.list$bothlimQQ] <-
          .bounds.list$dbothQQ * (emad/(1 + emad)) + .bounds.list$lbothQQ
      }
      
      if(any(.bounds.list$loweronlyQQ)){
        normalparamsX[.bounds.list$loweronlyQQ] <-
          exp(madparamsX[.bounds.list$loweronlyQQ]) + .bounds.list$loneQQ
      }
      
      if(any(.bounds.list$upperonlyQQ))
        normalparamsX[.bounds.list$upperonlyQQ] <-
          - exp(madparamsX[.bounds.list$upperonlyQQ]) + .bounds.list$uoneQQ
      
      if(exists(".ind.prof.phi"))
        if(is.nan(normalparamsX[.ind.prof.phi]))
          normalparamsX[.ind.prof.phi] <- 0
      
      return(normalparamsX)
    }
  
  newobjfunc <- function(madparams) {
    normalparams <-  madtonormalQQ(madparams)
    
    .objfuncQQ(normalparams)
    
  }
  
  startmadparams <- normaltomad(params)
  result <- nlm(newobjfunc, startmadparams, ...)
  result$madestimate <- result$estimate
  result$estimate <- madtonormalQQ(result$madestimate)
  remove(".bounds.list", pos=1, inherits=T)
  remove(".objfuncQQ", pos=1, inherits=T)
  remove("madtonormalQQ", pos=1, inherits=T)
  
###  return(result, madtonormalQQ(normaltomad(params)),params)
  return(result)
}

"plot.variogram" <-
  function (obj, max.dist = max(obj$u), ylim = "default",
            scaled = FALSE,
            var.lines = FALSE, envelope.obj = NULL,
            bin.cloud = FALSE, type = NULL, ...) 
{
  if(is.null(type)){
    if (obj$output.type == "bin") type <- "b"
    if (obj$output.type == "smooth") type <- "l"
    if (obj$output.type == "cloud") type <- "p"
  }
  if (bin.cloud == TRUE && type != "b") 
    stop("plot.variogram: object must be a binned variogram with option bin.cloud=T")
  if (bin.cloud == TRUE && all(is.na(obj$bin.cloud))) 
    stop("plot.variogram: object must be a binned variogram with option bin.cloud=T")
  if (bin.cloud == TRUE && any(!is.na(obj$bin.cloud))) 
    boxplot(obj$bin.cloud, varwidth = TRUE, names = as.character(obj$u), 
            xlab = "midpoints of distance class", ylab = paste("variogram values / ", 
                                                    obj$estimator.type, "estimator"), ...)
  else {
    u <- obj$u[obj$u <= max.dist]
    obj$v <- as.matrix(obj$v)
    v <- obj$v[obj$u <= max.dist, 1]
    ymax <- max(obj$v[obj$u <= max.dist, ])
    if (!is.null(envelope.obj)) 
      ymax <- max(envelope.obj$v.upper, ymax)
    if (scaled) v <- v/obj$var.mark
    if (ylim == "default") 
      plot(u, v, xlim = c(0, max.dist), ylim = c(0, ymax), 
           xlab = "Distance", ylab = "Semi-variance", type = type, ...)
    else
      plot(u, v, xlim = c(0, max.dist), ylim = ylim, xlab = "Distance", 
           ylab = "Semi-variance", type = type, ...)
    if (var.lines) {
      if (scaled) abline(h = 1, lty = 3)
      else abline(h = obj$var.mark, lty = 3)
    }
    if (ncol(obj$v) > 2) {
      lines.f <- function(v){lines(u,v, type=type)} 
      apply(obj$v[,-1],2, lines.f)
    }
    if (!is.null(envelope.obj)) {
      lines(u, envelope.obj$v.lower, lty = 4)
      lines(u, envelope.obj$v.upper, lty = 4)
    }
  }
  return(invisible())
}

"rfm.bin" <-
  function (cloud, l = 15, uvec = "default", nugget.tolerance = 0, 
            estimator.type = c("classical", "robust"), bin.cloud = FALSE,
            max.dist) 
{
  if (all(uvec == "default")) 
    uvec <- seq(0, max(cloud$u), l = l)
  estimator.type <- match.arg(estimator.type)
  ##  if(nugget.tolerance > 0) {
  ##    dnug <- mean(cloud$u[cloud$u <= nugget.tolerance])
  ##    cloud$u[cloud$u <= nugget.tolerance] <- 0
  ##    uvec <- uvec[uvec > nugget.tolerance]
  ##  }
  ##  u <- c(0, uvec)
  ##  n <- length(u)
  if(all(uvec == "default"))
    uvec <- seq(0, max.dist, l = 15)
  ubin <- c(0, uvec)
  nvec <- length(ubin)
  d <- 0.5 * diff(ubin[2:nvec])
  bins.lim <- c(0, (ubin[2:(nvec - 1)] + d), (d[nvec - 2] + ubin[
                                                                 nvec]))
  if(uvec[1] == 0 & nugget.tolerance == 0)
    uvec[1] <- (bins.lim[1] + bins.lim[2])/2
  if(nugget.tolerance > 0) {
    bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >
                                                nugget.tolerance])
    uvec <- c(0, (bins.lim[ - (1:2)] - 0.5 * diff(bins.lim)[
                                                            -1]))
  }
  nbins <- nc <- length(bins.lim) - 1
  if(is.null(max.dist))
    max.dist <- max(bins.lim)
  min.dist <- min(cloud$u)
  ##	d <- 0.5 * (u[3:n] - u[2:(n - 1)])
  ##	low <- c(0, (u[2:(n - 1)] + d))
  ##	high <- c((u[3:n] - d), (d[n - 2] + u[n]))
  ##	nc <- n - 1
  if (!is.matrix(cloud$v)) {
    vbin <- rep(0, nc)
    nbin <- rep(0, nc)
    sdbin <- rep(0, nc)
    if (bin.cloud == TRUE) 
      bins.clouds <- list()
    for (i in 1:nc) {
      ind <- (cloud$u > bins.lim[i]) & (cloud$u <= bins.lim[i+1])
      vbin[i] <- mean(cloud$v[ind])
      if (bin.cloud == TRUE) 
        bins.clouds[[i]] <- cloud$v[ind]
      nbin[i] <- sum(ind)
      if (estimator.type == "robust") 
        vbin[i] <- ((vbin[i])^4)/(0.914 + (0.988/nbin[i]))
      if (nbin[i] > 0) 
        sdbin[i] <- sqrt(var(cloud$v[ind]))
      else sdbin[i] <- NA
      NULL
    }
    if (uvec[1] == 0) 
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (min.dist == 0) {
      ind <- (cloud$u == 0)
      n.zero <- sum(ind)
      v.zero <- mean(cloud$v[ind])
      if (bin.cloud == TRUE) {
        bins.clouds[2:(length(bins.clouds) + 1)] <- bins.clouds[1:nc]
        bins.clouds[[1]] <- cloud$v[ind]
      }
      if (estimator.type == "robust") 
        v.zero <- ((v.zero)^4)/(0.914 + (0.988/n.zero))
      if (n.zero > 0) 
        sd.zero <- sqrt(var(cloud$v[ind]))
      else sd.zero <- NA
      uvec <- c(0, uvec)
      vbin <- c(v.zero, vbin)
      nbin <- c(n.zero, nbin)
      sdbin <- c(sd.zero, sdbin)
    }
    u <- uvec[!is.na(vbin)]
    v <- vbin[!is.na(vbin)]
    n <- nbin[!is.na(vbin)]
    sd <- sdbin[!is.na(vbin)]
    if (bin.cloud == TRUE) 
      bins.clouds <- bins.clouds[!is.na(vbin)]
  }
  else {
    if (bin.cloud == TRUE) 
      stop("option bins.cloud=T allowed only for 1 variable")
    nvcols <- ncol(cloud$v)
    vbin <- matrix(0, nrow = nc, ncol = nvcols)
    nbin <- rep(0, nc)
    sdbin <- matrix(0, nrow = nc, ncol = nvcols)
    for (i in 1:nc) {
      ind <- (cloud$u >= bins.lim[i]) & (cloud$u < bins.lim[i+1])
      nbin[i] <- sum(ind)
      for (j in 1:nvcols) {
        vbin[i, j] <- mean(cloud$v[ind, j])
        if (estimator.type == "robust") 
          vbin[i, j] <- ((vbin[i, j])^4)/(0.914 + (0.988/nbin[i]))
        if (nbin[i] > 0) 
          sdbin[i, j] <- sqrt(var(cloud$v[ind, j]))
        else sdbin[i, j] <- NA
      }
      NULL
    }
    if (uvec[1] == 0) 
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (min.dist == 0) {
      v.zero <- rep(0, nvcols)
      n.zero <- rep(0, nvcols)
      sd.zero <- rep(0, nvcols)
      for (j in 1:nvcols) {
        ind <- (cloud$u == 0)
        n.zero[j] <- sum(ind)
        v.zero[j] <- mean(cloud$v[ind, j])
        if (estimator.type == "robust") 
          v.zero[j] <- ((v.zero[j])^4)/(0.914 + (0.988/n.zero[j]))
        if (n.zero[j] > 0) 
          sd.zero[j] <- sqrt(var(cloud$v[ind, j]))
        else sd.zero[j] <- NA
        uvec <- c(0, uvec)
        vbin <- rbind(v.zero, vbin)
        nbin <- c(n.zero, nbin)
        sdbin <- rbind(sd.zero, sdbin)
      }
    }
    u <- uvec[!is.na(vbin[, 1])]
    n <- nbin[!is.na(vbin[, 1])]
    v <- matrix(0, nrow = length(u), ncol = nvcols)
    sd <- matrix(0, nrow = length(u), ncol = nvcols)
    for (j in 1:nvcols) {
      v[, j] <- vbin[!is.na(vbin[, j]), j]
      sd[, j] <- sdbin[!is.na(vbin[, j]), j]
    }
  }
  if (nugget.tolerance > 0) {
    u[1] <- nugget.tolerance
  }
  result <- list(u = u, v = v, n = n, sd = sd, output.type = "bin", bins.lim = bins.lim)
  if (!is.matrix(cloud$v) && bin.cloud == TRUE) 
    result$bin.cloud <- bins.clouds
  if (!is.null(class(cloud))) 
    class(result) <- class(cloud)
  return(result)
}

"points.geodata" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            data.col = 1, pt.sizes = c("data.proportional",
                            "rank.proportional", "quintiles",
                            "quartiles", "deciles", "equal"),
            cex.min, cex.max, pch.seq, col.seq, add.to.plot = FALSE,
            round.quantiles = FALSE, graph.pars = FALSE, ...) 
{
  pt.sizes <- match.arg(pt.sizes)
  if (add.to.plot == FALSE) {
    coords.lims <- apply(coords, 2, range)
    coords.diff <- diff(coords.lims)
    if (coords.diff[1] != coords.diff[2]) {
      coords.diff.diff <- abs(diff(as.vector(coords.diff)))
      ind.min <- which(coords.diff == min(coords.diff))
      coords.lims[, ind.min] <- coords.lims[, ind.min] + 
        c(-coords.diff.diff, coords.diff.diff)/2
    }
    par(pty = "s")
    plot(apply(coords, 2, range), type = "n", xlim = coords.lims[, 
                                                1], ylim = coords.lims[, 2], ...)
  }
  if (missing(cex.min)) 
    cex.min <- 0.5
  if (missing(cex.max)) 
    cex.max <- 1.5
  if (is.matrix(data)) 
    data <- as.vector(data[, data.col])
  graph.list <- list()
  if (pt.sizes == "quintiles" | pt.sizes == "quartiles" | pt.sizes == 
      "deciles") {
    if (pt.sizes == "quintiles") {
      n.quant <- 5
      if (missing(col.seq)) 
        col.seq <- c("blue", "green", "yellow", "orange3", "red2") 
    }
    if (pt.sizes == "quartiles") {
      n.quant <- 4
      if (missing(col.seq)) 
        col.seq <- c("blue", "green", "yellow", "red") 
    }
    if (pt.sizes == "deciles") {
      n.quant <- 10
      if (missing(col.seq)) 
        col.seq <- rainbow(13)[10:1]
    }
    if (missing(pch.seq)) 
      pch.seq <- rep(21, n.quant)
    cex.pt <- seq(cex.min, cex.max, l = n.quant)
    data.quantile <- quantile(data, probs = seq(0, 1, by = (1/n.quant)))
    if (round.quantiles == TRUE) {
      data.quantile[1] <- floor(data.quantile[1])
      data.quantile[n.quant + 1] <- ceiling(data.quantile[n.quant + 1])
      data.quantile <- round(data.quantile)
    }
    graph.list$quantiles <- data.quantile
    graph.list$cex <- cex.pt
    graph.list$col <- col.seq
    graph.list$pch <- pch.seq
    graph.list$data.group <- cut(data, breaks=data.quantile, include.l=T)
    if (add.to.plot) 
      points(coords, pch = 21, cex = cex.pt[as.numeric(graph.list$data.group)], bg = col.seq[as.numeric(graph.list$data.group)], ...)
    else
      points(coords, pch = 21, cex = cex.pt[as.numeric(graph.list$data.group)], bg = col.seq[as.numeric(graph.list$data.group)])
  }
  else {
    if (missing(pch.seq)) 
      pch.seq <- 21
    if (missing(col.seq)) 
      col.seq <- 0
    n <- length(data)
    coords.order <- coords[order(data), ]
    data.order <- data[order(data)]
    if (pt.sizes == "rank.proportional") {
      data.quantile <- range(data.order)
      size <- seq(cex.min, cex.max, l = n)
      graph.list$cex <- range(size)
      graph.list$pch <- unique(range(pch.seq))
      graph.list$col <- col.seq
      if (length(col.seq) == 1) 
        col.seq <- rep(col.seq, n)
      for (i in 1:n) {
        if (add.to.plot) 
          points(coords.order[i, , drop = FALSE], cex = size[i], 
                 pch = pch.seq, bg = col.seq[i], ...)
        else points(coords.order[i, , drop = FALSE], 
                    cex = size[i], pch = pch.seq, bg = col.seq[i])
      }
    }
    if (pt.sizes == "data.proportional") {
      r.y <- range(data.order)
      size <- cex.min + ((data.order - r.y[1]) * (cex.max - 
                                                  cex.min))/(r.y[2] - r.y[1])
      graph.list$cex <- c(cex.min, cex.max)
      graph.list$pch <- unique(range(pch.seq))
      graph.list$col <- col.seq
      if (length(col.seq) == 1) 
        col.seq <- rep(col.seq, n)
      for (i in 1:n) {
        if (add.to.plot) 
          points(coords.order[i, , drop = FALSE], cex = size[i], 
                 pch = pch.seq, bg = col.seq[i], ...)
        else points(coords.order[i, , drop = FALSE], 
                    cex = size[i], pch = pch.seq, bg = col.seq[i])
      }
    }
    if (pt.sizes == "equal") {
      if (add.to.plot) 
        points(coords, pch = pch.seq, bg = col.seq, cex = cex.max, 
               ...)
      else points(coords, pch = pch.seq, bg = col.seq, 
                  cex = cex.max)
    }
  }
  if (graph.pars == TRUE) 
    return(graph.list)
  else return(invisible())
}

