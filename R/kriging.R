"krige.conv" <-
  function (geodata, coords=geodata$coords, data=geodata$data,
            locations, borders = NULL, krige, output)
{
  if(missing(geodata))
    geodata <- list(coords = coords, data = data)
  call.fc <- match.call()
  base.env <- sys.frame(sys.nframe())
  ##
  ## reading input
  ##
  if(missing(krige))
    krige <- krige.control()
  else{
    ##    if(is.null(class(krige)) || class(krige) != "krige.geoR"){
    if(length(class(krige)) == 0 || class(krige) != "krige.geoR"){
      if(!is.list(krige))
        stop("krige.conv: the argument krige only takes a list or an output of the function krige.control")
      else{
        krige.names <- c("type.krige","trend.d","trend.l","obj.model",
                         "beta","cov.model", "cov.pars",
                         "kappa","nugget","micro.scale","dist.epsilon",
                         "lambda","aniso.pars")
        krige.user <- krige
        krige <- list()
        if(length(krige.user) > 0){
          for(i in 1:length(krige.user)){
            n.match <- match.arg(names(krige.user)[i], krige.names)
            krige[[n.match]] <- krige.user[[i]]
          }
        }
        if(is.null(krige$type.krige)) krige$type.krige <- "ok"  
        if(is.null(krige$trend.d)) krige$trend.d <-  "cte"
        if(is.null(krige$trend.l)) krige$trend.l <-  "cte"
        if(is.null(krige$obj.model)) krige$obj.model <-  NULL
        if(is.null(krige$beta)) krige$beta <- NULL 
        if(is.null(krige$cov.model)) krige$cov.model <- "matern"  
        if(is.null(krige$cov.pars))
          stop("covariance parameters (sigmasq and phi) should be provided in cov.pars")
        if(is.null(krige$kappa)) krige$kappa <-  0.5
        if(is.null(krige$nugget)) krige$nugget <-  0
        if(is.null(krige$micro.scale)) krige$micro.scale <- 0  
        if(is.null(krige$dist.epsilon)) krige$dist.epsilon <-  1e-10
        if(is.null(krige$aniso.pars)) krige$aniso.pars <- NULL  
        if(is.null(krige$lambda)) krige$lambda <- 1 
        krige <- krige.control(type.krige = krige$type.krige,
                               trend.d = krige$trend.d,
                               trend.l = krige$trend.l,
                               obj.model = krige$obj.model,
                               beta = krige$beta,
                               cov.model = krige$cov.model,
                               cov.pars = krige$cov.pars,
                               kappa = krige$kappa,
                               nugget = krige$nugget,
                               micro.scale = krige$micro.scale,
                               dist.epsilon = krige$dist.epsilon, 
                               aniso.pars = krige$aniso.pars,
                               lambda = krige$lambda)
        
      }
    }
  }
  cov.model <- krige$cov.model
  kappa <- krige$kappa
  lambda <- krige$lambda
  beta <- krige$beta
  cov.pars <- krige$cov.pars
  nugget <- krige$nugget
  micro.scale <- krige$micro.scale
  aniso.pars <- krige$aniso.pars
  ##
  ## reading output options
  ##
  if(missing(output))
    output <- output.control()
  else{
    ##    if(is.null(class(output)) || class(output) != "output.geoR"){
    if(length(class(krige)) == 0 || class(output) != "output.geoR"){
      if(!is.list(output))
        stop("krige.conv: the argument output can take only a list or an output of the function output.control")
      else{
        output.names <- c("n.posterior","n.predictive","moments","n.back.moments","simulations.predictive",
                          "mean.var","quantile","threshold","signal","messages.screen")
        output.user <- output
        output <- list()
        if(length(output.user) > 0){
          for(i in 1:length(output.user)){
            n.match <- match.arg(names(output.user)[i], output.names)
            output[[n.match]] <- output.user[[i]]
          }
        }
        if(is.null(output$n.posterior)) output$n.posterior <- 1000 
        if(is.null(output$n.predictive)) output$n.predictive <- NULL
        if(is.null(output$moments)) output$moments <- TRUE
        if(is.null(output$n.back.moments)) output$n.back.moments <- 1000 
        if(is.null(output$simulations.predictive)){
          if(is.null(output$n.predictive)) output$simulations.predictive <- NULL
          else
            output$simulations.predictive <- ifelse(output$n.predictive > 0, TRUE, FALSE)
        }
        if(is.null(output$mean.var)) output$mean.var <- NULL
        if(is.null(output$quantile)) output$quantile <- NULL
        if(is.null(output$threshold)) output$threshold <- NULL
        if(is.null(output$signal)) output$signal <- NULL
        if(is.null(output$messages.screen)) output$messages.screen <- TRUE
        output <- output.control(n.posterior = output$n.posterior,
                                 n.predictive = output$n.predictive,
                                 moments = output$moments,
                                 n.back.moments = output$n.back.moments, 
                                 simulations.predictive = output$simulations.predictive,
                                 mean.var = output$mean.var,
                                 quantile = output$quantile,
                                 threshold = output$threshold,
                                 signal = output$signal,
                                 messages.screen = output$messages.screen)
      }
    }
  }
  signal <- ifelse(is.null(output$signal), FALSE, output$signal)
  messages.screen <- output$messages.screen
  n.predictive <- output$n.predictive
  n.back.moments <- output$n.back.moments
  ##
  n.predictive <- ifelse(is.null(n.predictive), 0, n.predictive)
  simulations.predictive <- ifelse(is.null(output$simulations.predictive), FALSE, TRUE)
  keep.simulations <- ifelse(is.null(output$keep.simulations), TRUE, FALSE)
  mean.estimator <- output$mean.estimator
  if(is.null(mean.estimator) & simulations.predictive) mean.estimator <- TRUE
  quantile.estimator <- output$quantile.estimator
  probability.estimator <- output$probability.estimator
  if(!is.null(probability.estimator)){
    if(length(probability.estimator) > 1 &
       length(probability.estimator) != nrow(locations))
      stop("krige.conv: probability.estimator must either have length 1, or have length = nrow(locations)\n")
  }
  if(simulations.predictive & n.predictive == 0) n.predictive <- 1000
  ##
  ## checking input
  ##
  if(krige$type.krige == "ok") beta.prior <- "flat"
  if(krige$type.krige == "sk") beta.prior <- "deg"
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning("krige.conv: coordinates provided as a vector, assuming one spatial dimension")
  }
  coords <- as.matrix(coords)
  if(is.vector(locations)) {
    if(length(locations) == 2) {
      locations <- t(as.matrix(locations))
      if(messages.screen) 
        warning("krige.conv: assuming that there is only 1 prediction point")
    }
    else{
      warning("krige.conv: locations provided as a vector, assuming one spatial dimension")
      locations <- as.matrix(cbind(locations, 0))
    }
  }
  else locations <- as.matrix(locations)
  ##
  ## selecting locations inside the borders 
  ##
  if(!is.null(borders)){
    locations <- locations.inside(locations, borders)
    if(nrow(locations) == 0)
      stop("\nkrige.conv: there are no prediction locations inside the borders")
    if(messages.screen)
      cat("krige.conv: results will be returned only for prediction locations inside the borders\n")
  }
  dimnames(coords) <- list(NULL, NULL)
  dimnames(locations) <- list(NULL, NULL)
  ##
  ## Checking for 1D prediction 
  ##
  if(length(unique(locations[,1])) == 1 | length(unique(locations[,2])) == 1)
    krige1d <- TRUE
  else krige1d <- FALSE
  ##
  ## building the trend matrix
  ##
  if(messages.screen){
    if(is.numeric(krige$trend.d))
      cat("krige.conv: model with covariates matrix provided by the user")
    else
      cat(switch(as.character(krige$trend.d)[1],
                 "cte" = "krige.conv: model with constant mean",
                 "1st" = "krige.conv: model with mean given by a 1st order polynomial on the coordinates",
                 "2nd" = "krige.conv: model with mean given by a 2nd order polynomial on the coordinates",
                 "krige.conv: model with mean defined by covariates provided by the user"))
    cat("\n")
  }
  trend.d <- unclass(trend.spatial(trend=krige$trend.d, geodata = geodata))
  if (nrow(trend.d) != nrow(coords)) 
      stop("coords and trend.d have incompatible sizes")
  beta.size <- ncol(trend.d)
  if(beta.prior == "deg")
    if(beta.size != length(beta))
      stop("size of mean vector is incompatible with trend specified") 
  trend.l <- unclass(trend.spatial(trend=krige$trend.l, geodata = list(coords = locations)))
  if (nrow(trend.l) != nrow(locations)) 
    stop("locations and trend.l have incompatible sizes")
  if(beta.size > 1)
    beta.names <- paste("beta", (0:(beta.size-1)), sep="")
  else beta.names <- "beta"
  ##
  ## Anisotropy correction (this should be placed AFTER trend.d/trend.l
  ##
  if(!is.null(aniso.pars)) {
#    if((abs(aniso.pars[1]) > 0.001) & (abs(aniso.pars[2] - 1) > 0.001)){
    if(abs(aniso.pars[2] - 1) > 0.0001){
      if(messages.screen)
        cat("krige.conv: anisotropy correction performed\n")
      coords <- coords.aniso(coords = coords, aniso.pars = aniso.pars)
      locations <- coords.aniso(coords = locations, aniso.pars = aniso.pars)
    }
  }
  ##
  ## Box-Cox transformation
  ##
  if(abs(lambda - 1) > 0.001) {
    if(messages.screen) cat("krige.conv: performing the Box-Cox data transformation\n")
    data <- BCtransform(data, lambda = lambda)$data
  }
  ## 
  ## setting covariance parameters
  ##
  if(is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    cpars <- c(1, phi)
  }
  else {
    stop("current version of krige.conv does not accept nested covariance models\n") 
    ##    sigmasq <- cov.pars[, 1]
    ##    phi <- cov.pars[, 2]
    ##    cpars <- cbind(1, phi)
  }
  ##  sill.partial <- micro.scale + sum(sigmasq)
  sill.partial <- sum(sigmasq)
  tausq.rel <- nugget/sum(sigmasq)
  tausq.rel.micro <- micro.scale/sum(sigmasq)
  n <- length(data)
  ni <- nrow(trend.l)
  ##
  ## starting kriging calculations
  ##
  kc <- list()
  invcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                           kappa = kappa, nugget = tausq.rel,
                           cov.pars = cpars, inv = TRUE,
                           only.inv.lower.diag = TRUE)
  ittivtt <- solve.geoR(bilinearformXAY(X = as.vector(trend.d),
                                        lowerA = as.vector(invcov$lower.inverse),
                                        diagA = as.vector(invcov$diag.inverse), 
                                        Y = as.vector(trend.d)))
  if(beta.prior == "flat"){
    beta.flat <- drop(ittivtt %*% bilinearformXAY(X = as.vector(trend.d),
                                                  lowerA = as.vector(invcov$lower.inverse),
                                                  diagA = as.vector(invcov$diag.inverse), 
                                                  Y = as.vector(data)))
  }
  v0 <- loccoords(coords = coords, locations = locations)
  if(n.predictive > 0){
    ## checking if there are data points coincident with prediction locations
    loc.coincide <- apply(v0, 2, function(x, min.dist){any(x < min.dist)},
                          min.dist=krige$dist.epsilon)
    if(any(loc.coincide)) loc.coincide <- (1:ni)[loc.coincide]
    else loc.coincide <- NULL
    if(!is.null(loc.coincide)){
      temp.f <- function(x, data, dist.eps){return(data[x < dist.eps])}
      data.coincide <- apply(v0[, loc.coincide, drop=FALSE], 2, temp.f, data=data, dist.eps=krige$dist.epsilon)
    }
    else data.coincide <- NULL
  }
  else remove("locations")
  ## using nugget interpreted as microscale variation or measurement error
  nug.factor <- ifelse(signal, tausq.rel.micro, tausq.rel)
  ## covariances between data and prediction locations
  v0 <- ifelse(v0 < krige$dist.epsilon, 1+nug.factor,
               cov.spatial(obj = v0, cov.model = cov.model, 
                           kappa = kappa, cov.pars = cpars))
  tv0ivv0 <- diagquadraticformXAX(X = as.vector(v0),
                                  lowerA = invcov$lower.inverse,
                                  diagA = invcov$diag.inverse)
  b <- bilinearformXAY(X = as.vector(cbind(data,trend.d)),
                       lowerA = as.vector(invcov$lower.inverse),
                       diagA = as.vector(invcov$diag.inverse), 
                       Y = as.vector(v0))
  if(n.predictive == 0) remove("v0","invcov")
  tv0ivdata <- drop(b[1,])
  b <- t(trend.l) -  b[-1,, drop=FALSE]
  if(beta.prior == "deg") {
    kc$predict <- tv0ivdata + drop(crossprod(b,beta))
    kc$krige.var <- sill.partial * drop(1+nug.factor - tv0ivv0)
    beta.est <- "Simple kriging performed (beta provided by user)"
  }
  if(beta.prior == "flat"){
    kc$predict <- tv0ivdata + drop(crossprod(b,beta.flat))
    if(beta.size == 1)
      bitb <- drop(b^2) * drop(ittivtt)
    else
      bitb <- diagquadraticformXAX(X = b,
                                   lowerA = (ittivtt[lower.tri(ittivtt)]),
                                   diagA = diag(ittivtt))
    kc$krige.var <- sill.partial * drop(1+nug.factor - tv0ivv0 + bitb)
    kc$beta.est <- beta.flat
    names(kc$beta.est) <- beta.names
    remove("bitb")
  }
  remove("tv0ivv0", "tv0ivdata")
  if(n.predictive == 0) remove("b")
  kc$distribution <- "normal"
  if(any(round(kc$krige.var, dig=12) < 0))
    cat("krige.conv: negative kriging variance found! Investigate why this is happening.\n")
  ##
  ## ########### Sampling from the resulting distribution ###
  ##
  if(n.predictive > 0) {
    seed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
    if(messages.screen)
      cat("krige.conv: sampling from the predictive distribution (conditional simulations)\n")
    if(length(cov.pars) > 2){
      reduce.var <-
        bilinearformXAY(X = as.vector(v0),
                        lowerA = as.vector(invcov$lower.inverse),
                        diagA = as.vector(invcov$diag.inverse), 
                        Y = as.vector(v0))
      remove("v0")
      attr(reduce.var, "dim") <- c(ni, ni)
      if(beta.prior == "flat"){
        if(beta.size == 1)
          ok.add.var <- outer(as.vector(b),as.vector(b)) * as.vector(ittivtt)
        else
          ok.add.var <-
            bilinearformXAY(X = as.vector(b),
                            lowerA = as.vector(ittivtt[lower.tri(ittivtt)]),
                            diagA = as.vector(diag(ittivtt)), 
                            Y = as.vector(b))
        reduce.var <- reduce.var + ok.add.var
      }
      varcov <- (varcov.spatial(coords = locations,
                                cov.model = cov.model,
                                cov.pars = cov.pars,
                                kappa = kappa, nugget = nugget)$varcov) -
                                  reduce.var
      remove("reduce.var")
      if(is.R()) gc(verbose=FALSE)
      kc$simulations <-  kc$predict +
        crossprod(chol(varcov), matrix(rnorm(ni * n.predictive),
                                       ncol=n.predictive))
    }
    else{
      coincide.cond <- (((round(1e12 * nugget) == 0) | !signal) & (!is.null(loc.coincide)))
      if(coincide.cond){
        nloc <- ni - length(loc.coincide)
        ind.not.coincide <- -(loc.coincide) 
        v0 <- v0[,ind.not.coincide, drop=FALSE]
        b <- b[,ind.not.coincide, drop=FALSE]
      }
      else{
        nloc <- ni
        ind.not.coincide <- TRUE
      }
      Dval <- 1.0 + nug.factor
      if(beta.prior == "deg")
        vbetai <- matrix(0, ncol = beta.size, nrow = beta.size)
      else
        vbetai <- matrix(ittivtt, ncol = beta.size, nrow = beta.size)
      df.model <- ifelse(beta.prior == "deg", n, n-beta.size)
      kc$simulations <- matrix(NA, nrow=ni, ncol=n.predictive)
      if(nloc > 0)
        kc$simulations[ind.not.coincide,] <- 
          cond.sim(env.loc = base.env, env.iter = base.env,
                   loc.coincide = loc.coincide,
                   coincide.cond = coincide.cond, 
                   tmean = kc$predict[ind.not.coincide],
                   Rinv = invcov,
                   mod = list(beta.size = beta.size, nloc = nloc,
                     Nsims = n.predictive, n = n, Dval = Dval,
                     df.model = df.model, s2 = sill.partial,
                     cov.model.number = cor.number(cov.model),
                     phi = phi, kappa = kappa),
                   vbetai = vbetai,
                   fixed.sigmasq = TRUE)
      remove("v0", "b", "locations", "invcov")
      if(is.R()) gc(verbose = FALSE)
      if(coincide.cond)
        kc$simulations[loc.coincide,] <- rep(data.coincide, n.predictive)
    }
    ##
    ## Backtransforming simulations
    ##
    if(abs(lambda - 1) > 0.001){
      if(messages.screen)
        cat("krige.conv: back-transforming the simulated values\n")
      if(any(kc$simulations < -1/lambda))
        warning("Truncation in the back-transformation: there are simulated values less than (- 1/lambda) in the normal scale.")
      kc$simulations <-
        BCtransform(kc$simulations, lambda, inv=TRUE)$data
    }
    ##
    ## mean/quantiles/probabilities estimators from simulations
    ##
    if(!is.null(mean.estimator) | !is.null(quantile.estimator) |
       !is.null(probability.estimator)){
      kc <- c(kc, statistics.predictive(simuls= kc$simulations,
                                        mean.var = mean.estimator,
                                        quantile = quantile.estimator,
                                        threshold = probability.estimator))
    }
    kc$.Random.seed <- seed
  }
  ##
  ## Backtransforming moments of the prediction distribution
  ## NOTE: this must be placed here, AFTER the simulations
  ##
  if(abs(lambda-1) > 0.001){
    if(messages.screen){
      cat("krige.conv: back-transforming the predicted mean and variance\n")
      if((abs(lambda) > 0.001) & (abs(lambda-0.5) > 0.001))
        cat("krige.conv: back-transforming by simulating from the predictive.\n           (run the function a few times and check stability of the results.\n")
    }
    kc[c("predict", "krige.var")] <-
      backtransform.moments(lambda = lambda,
                            mean = kc$predict,
                            variance = kc$krige.var,
                            distribution = "normal",
                            n.simul = n.back.moments)[c("mean", "variance")]
  }
  ##
  message <- "krige.conv: Kriging performed using global neighbourhood"
  if(messages.screen) cat(paste(message, "\n"))
  ##
  kc$message <-  message
  kc$call <- call.fc
  ##
  ## Setting classes and attributes 
  ##
  attr(kc, 'sp.dim') <- ifelse(krige1d, "1d", "2d")
  attr(kc, "prediction.locations") <- call.fc$locations
  if(!is.null(call.fc$borders))
    attr(kc, "borders") <- call.fc$borders
  class(kc) <- "kriging"
  return(kc)
}

"krige.control" <-
  function(type.krige = "ok",
           trend.d = "cte", trend.l = "cte",
           obj.model = NULL,
           beta, cov.model, cov.pars, kappa,
           nugget, micro.scale = 0, dist.epsilon = 1e-10, 
           aniso.pars, lambda)
{
  if(type.krige != "ok" & type.krige != "OK" & type.krige != "o.k." & type.krige != "O.K." & type.krige != "sk" & type.krige != "SK" & type.krige != "s.k." & type.krige != "S.K.")
    stop("krige.conv: wrong option in the argument type.krige. It should be \"sk\" or \"ok\"(if ordinary or simple kriging is to be performed)")
  if(type.krige=="OK" | type.krige=="O.K." |type.krige=="o.k.")
    type.krige <- "ok"
  if(type.krige=="SK" | type.krige=="S.K." |type.krige=="s.k.")
    type.krige <- "sk"
  ##
  if(!is.null(obj.model)){
    if(missing(beta)) beta <- obj.model$beta
    if(missing(cov.model)) cov.model <- obj.model$cov.model
    if(missing(cov.pars)) cov.pars <- obj.model$cov.pars
    if(missing(kappa)) kappa <- obj.model$kappa
    if(missing(nugget)) nugget <- obj.model$nugget
    if(missing(lambda)) lambda <- obj.model$lambda
    if(missing(aniso.pars)) aniso.pars <- obj.model$aniso.pars
  }
  else{
    if(missing(beta)) beta <- NULL
    if(missing(cov.model)) cov.model <- "matern"
    if(missing(cov.pars))
      stop("covariance parameters (sigmasq and phi) should be provided")
    if(missing(kappa)) kappa <- 0.5
    if(missing(nugget)) nugget <- 0
    if(missing(lambda)) lambda <- 1
    if(missing(aniso.pars)) aniso.pars <- NULL
  }
  ##
  if(type.krige == "sk")
    if(is.null(beta) | !is.numeric(beta))
      stop("\nkrige.conv: argument beta must be provided in order to perform simple kriging")
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic",
                           "wave", "linear", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(micro.scale > nugget)
    stop("krige.control: micro.scale must be in the interval [0, nugget]")
  ##
  if(!is.null(aniso.pars))
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("krige.control: anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
  ##
  if(inherits(trend.d, "formula") | inherits(trend.l, "formula")){
    if((inherits(trend.d, "formula") == FALSE) | (inherits(trend.l, "formula") == FALSE))
      stop("krige.control: trend.d and trend.l must have similar specification")
  }
  else{
##    if((!is.null(class(trend.d)) && class(trend.d) == "trend.spatial") &
##       (!is.null(class(trend.l)) && class(trend.l) == "trend.spatial")){
    if((length(class(trend.d)) > 0 && class(trend.d) == "trend.spatial") &
       (length(class(trend.l)) > 0 && class(trend.l) == "trend.spatial")){
      if(ncol(trend.d) != ncol(trend.l))
        stop("krige.bayes: trend.d and trend.l do not have the same number of columns")
    }
    else
      if(trend.d != trend.l)
        stop("krige.control: trend.l is different from trend.d")
  }
  ##
  res <- list(type.krige = type.krige,
              trend = trend.d, trend.d = trend.d, trend.l = trend.l, 
              beta = beta,
              cov.model = cov.model, 
              cov.pars = cov.pars, kappa = kappa,
              nugget = nugget,
              micro.scale = micro.scale, dist.epsilon = dist.epsilon, 
              aniso.pars = aniso.pars, lambda = lambda)
  class(res) <- "krige.geoR"
  return(res)
}

"prepare.graph.kriging" <-
  function (locations, borders, values, xlim, ylim) 
{
  if (!is.null(borders)) {
    borders <- as.matrix(as.data.frame(borders))
    require(splancs)
    inout.vec <- as.vector(inout(pts = locations, poly = borders))
    if (length(inout.vec) != length(values)) 
      stop("image.kriging: length of the argument values is incompatible with number of elements inside the borders.")
    temp <- rep(NA, nrow(locations))
    temp[inout.vec] <- values[inout.vec]
    values <- temp
    remove("temp")
  }
  locations <- locations[order(locations[, 2], locations[, 1]), ]
  x <- as.numeric(levels(as.factor(round(locations[, 1], dig = 8))))
  nx <- length(x)
  y <- as.numeric(levels(as.factor(round(locations[, 2], dig = 8))))
  ny <- length(y)
  if (missing(xlim))  xlim <- NULL
  if (missing(ylim))  ylim <- NULL
  coords.lims <- set.coords.lims(coords = locations, xlim = xlim, 
                                 ylim = ylim)
  coords.lims[, 1] <- coords.lims[, 1] + c(-0.025, 0.025) * 
    diff(coords.lims[, 1])
  coords.lims[, 2] <- coords.lims[, 2] + c(-0.025, 0.025) * 
    diff(coords.lims[, 2])
  return(list(x = x, y = y, values = matrix(values, ncol = ny), 
              coords.lims = coords.lims))
}


"image.kriging" <-
  function (x, locations, borders, 
            values = x$predict, coords.data,
            xlim, ylim, x.leg, y.leg, ...) 
{
  pty.prev <- par()$pty
  ldots <- list(...)
  if(missing(x)) x <- NULL
  attach(x)
  on.exit(detach(x))
  if(missing(locations)) locations <-  eval(attr(x, "prediction.locations"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(missing(borders)){
    if(!is.null(attr(x, "borders"))) borders <- eval(attr(x, "borders"))
    else borders <- NULL
  }
  if(missing(coords.data)) coords.data <- NULL
  if(missing(xlim)) xlim <- NULL
  if(missing(ylim)) ylim <- NULL
  if(missing(x.leg)) x.leg <- NULL
  if(missing(y.leg)) y.leg <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D')
    plot.1d(values, xlim=xlim, ylim = ylim,
            x1vals = unique(round(locations[,1], dig=12)), ...)
  else{
    locations <- prepare.graph.kriging(locations=locations,
                                       borders=borders, values=values,
                                       xlim = xlim, ylim = ylim) 
    par(pty = "s")
    image(x=locations$x, y=locations$y, z=locations$values,
          xlim = locations$coords.lims[,1],
          ylim = locations$coords.lims[,2], ...)
    ##
    ## adding points at data locations
    ##
    if(!is.null(coords.data)) points(coords.data, pch=20)
    ##
    ## adding borders
    ##
    if(!is.null(borders)) polygon(borders, lwd=2)
    ##
    ## adding the legend
    ##
    if(!is.null(x.leg) & !is.null(y.leg)){
      if(is.null(ldots$col)) ldots$col <- heat.colors(12)
      legend.krige(x.leg=x.leg, y.leg=y.leg,
                   values=locations$values[!is.na(locations$values)], ...)
    }
  }
  par(pty = pty.prev)
  return(invisible())
}

"persp.kriging" <-
  function(x, locations, borders, values = x$predict, ...)
{
  if(missing(x)) x <- NULL
  attach(x)
  on.exit(detach(x))
  if(missing(locations)) locations <-  eval(attr(x, "prediction.locations"))
  if(is.null(locations)) stop("prediction locations must be provided")
  if(ncol(locations) != 2)
    stop("locations must be a matrix or data-frame with two columns")
  if(missing(borders)) borders <- NULL
  ##
  ## Plotting 1D or 2D
  ##
  if(!is.null(attr(x, 'sp.dim')) && attr(x, 'sp.dim') == '1D')
    plot.1d(values, xlim=xlim, ylim = ylim,
            x1vals = unique(round(locations[,1], dig=12)), ...)
  else{
    locations <- prepare.graph.kriging(locations=locations,
                                       borders=borders, values=values) 
    persp(locations$x, locations$y, locations$values, ...)
  }
  return(invisible())
}

"legend.krige" <-
  function(x.leg, y.leg, values, scale.vals, vertical = FALSE, offset.leg = 1, ...)
{
  values <- values[!is.na(values)]
  if(length(x.leg) != 2 | length(y.leg) != 2)
    stop("x.leg and y.leg require a vector with 2 elements")
  v.r <- range(values[is.finite(values)], na.rm = TRUE)
  lags.x <- function(xs, nl){
    xs.r <- 0.5 * diff(xs/(nl-1))
    return(seq(xs[1]+xs.r, xs[2]-xs.r, l=nl))
  }
  leg.l <- list(...)
  if(is.null(leg.l$br))
    nc <- ifelse(is.null(leg.l$col), 12, length(leg.l$col))
  else
    nc <- length(leg.l$breaks) - 1
  if(is.null(leg.l$col)) leg.l$col <- heat.colors(nc)
  if(is.null(leg.l$zl)) leg.l$zlim <- c(v.r[1], v.r[2])
  if(vertical){
    xy <- list(x=x.leg, y=lags.x(xs=y.leg, nl=nc))
    if(is.null(leg.l$br))
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), nrow=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col)
    else
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), nrow=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col, breaks = leg.l$br)
  }
  else{
    xy <- list(x=lags.x(xs=x.leg, nl=nc), y=y.leg)
    if(is.null(leg.l$br))
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), ncol=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col)
    else
      image(x=xy$x, y=xy$y,
            z=matrix(seq(leg.l$zlim[1], leg.l$zlim[2], l=nc), ncol=1),
            add=TRUE, xaxs = "i", yaxs = "i", xlab="", ylab="",
            zlim = leg.l$zlim, col=leg.l$col, breaks = leg.l$br)
  }
  leg.poly <- rbind(c(x.leg[1], y.leg[1]), c(x.leg[2], y.leg[1]),
                    c(x.leg[2], y.leg[2]), c(x.leg[1], y.leg[2]),
                    c(x.leg[1], y.leg[1]))
  polygon(leg.poly)
#  if(is.null(leg.l$cex)) leg.l$cex <- par()$cex
  if(is.null(leg.l$cex)) leg.l$cex <- 0.8
  if(missing(scale.vals))
    scale.vals <- pretty(c(values,leg.l$zlim), n=5, min.n=4)
  scale.vals <- scale.vals[scale.vals > leg.l$zlim[1] &
                           scale.vals < leg.l$zlim[2]]
  if(vertical){
    y.r <- range(lags.x(xs=y.leg,nl=nc))
    y.text <- y.r[1] + ((scale.vals - leg.l$zlim[1]) * diff(y.r))/diff(leg.l$zlim)
    text((max(x.leg)+ offset.leg * diff(x.leg)), y.text,
         lab=scale.vals, col=1, cex=leg.l$cex)
  }
  else{
    x.r <- range(lags.x(xs=x.leg,nl=nc))
    x.text <- x.r[1] + ((scale.vals - leg.l$zlim[1]) * diff(x.r))/diff(leg.l$zlim)
    text(x.text, (max(y.leg)+ offset.leg * (diff(y.leg)/2)), lab=scale.vals, col=1, cex=leg.l$cex)
  }
  return(invisible())
}

