"krige.bayes" <- 
  function(geodata, coords=geodata$coords, data=geodata$data, locations = "no",
           model = model.control(
             trend.d = "cte", trend.l = "cte",
             cov.model = "matern",
             kappa = 0.5, aniso.pars = NULL, lambda = 1), 
           prior = prior.control(
             beta.prior = c("flat", "normal", "fixed"),
             beta = NULL, beta.var = NULL,
             sill.prior = c("reciprocal", "fixed"), sill = NULL, 
             range.prior = c("uniform", "exponential", "fixed",
               "squared.reciprocal","reciprocal"), exponential.prior.par = 1,
             range = NULL, range.discrete = NULL, 
             nugget.prior = c("fixed", "uniform"), nugget = 0,
             nugget.discrete = NULL),
           output = output.control(
             n.posterior = 1000, n.predictive = NULL, moments = TRUE,
             simulations.predictive = TRUE, keep.simulations = TRUE,
             mean.estimator = TRUE, quantile.estimator = NULL,
             probability.estimator = NULL, signal = TRUE, 
             messages.screen = TRUE)
           )
{
##           smaller.locations = NULL, trend.smaller.l = model$trend.l
  ##
######################### PART 1 ##############################
  ## Input check
  ##
  if(is.R()) require(mva)
  call.fc <- match.call()
  seed <- .Random.seed
  kb.results <- list(posterior = list(), predictive=list())
  ##
  ## reading input
  ##
  cov.model <- model$cov.model
  cov.model.number <- cor.number(cov.model)
  kappa <- model$kappa
  if(cov.model == "powered.exponential" & (kappa <= 0 | kappa > 2))
    stop("for power exponential correlation model the parameter kappa must be in the interval \(0,2\]")
  lambda <- model$lambda
  ##
  beta <- prior$beta
  beta.var <- prior$beta.var
  sill <- prior$sill
  ##range <- prior$range
  nugget <- prior$nugget
  exponential.prior.par <- prior$exponential.prior.par  
  ##
  n.posterior <- output$n.posterior
  n.predictive <- output$n.predictive
  moments <- output$moments    
  messages.screen <- output$messages.screen    
  ##
  ## checking data configuration
  ##
  if(is.vector(coords)){
    coords <- cbind(coords, 0)
    warning("vector of coordinates: one spatial dimention assumed")
  }
  coords <- as.matrix(coords)
  data.dist<- as.vector(dist(coords))
  data.dist.range <- range(data.dist)
  data.dist.min <- data.dist.range[1]
  data.dist.max <- data.dist.range[2]
  ## check == here
  if(round(1e12*data.dist.min) == 0)
    stop("This function does not allow two data at same location")
  if(all(locations == "no")) {
    if(prior$beta.prior != "fixed" & prior$sill.prior != "fixed"  & prior$range.prior != "fixed")
      if(messages.screen){
        cat("krige.bayes: no prediction locations provided.")
        cat("             Only Bayesian estimates of model parameters will be returned. \n ")
      }
      else
        stop("Locations to be predicted not provided. Bayesian parameter estimation allowed only if priors are provided for beta, sill and range")
  }
  ##  if(!is.null(smaller.locations)){
  ##    if(if(is.R()) require(akima) == FALSE)
  ##      stop("package \"akima\" is needed if a non-null object is provided in the argument \"smaller.grid\"")
  ##  }
  ##
  ## Cheking priors input
  ##
  if(prior$beta.prior == "fixed" & prior$sill.prior != "fixed")
    stop("option for fixed beta and random sill not implemented\n")
  if(prior$beta.prior == "fixed" & prior$range.prior != "fixed")
    stop("option for fixed beta and random range not implemented\n")
  if(prior$sill.prior == "fixed" & prior$range.prior != "fixed")
    stop("option for fixed sill and random range not implemented\n")
  if(prior$beta.prior != "flat" & prior$sill.prior != "fixed")
    stop("selected prior for beta not allowed for selected prior for sill\n")
  if(prior$beta.prior == "fixed" & is.null(beta))
    stop("if beta is fixed, its value must be provided in the argument beta\n")
  if(prior$beta.prior == "normal" & (is.null(beta) | is.null(beta.var)))
    stop("if the prior for beta is normal, the prior mean(s) and prior (co)variance must be provided using the argument beta and priorr$beta.var\n")
  if(prior$sill.prior == "fixed" & is.null(sill))
    stop("if sill is fixed, its value must be provided in the argument sill\n")
  if(prior$range.prior == "fixed" & is.null(prior$range))
    stop("if range is fixed, its value must be provided in the argument range\n")
  if(prior$range.prior != "fixed"){
    if (is.null(prior$range.discrete))
      stop("to include the range as random in the Bayesian analysis the argument range.discrete must be provided\n")
    discrete.diff <- diff(prior$range.discrete)
    if(round(max(1e08 * discrete.diff)) != round(min(1e08 * discrete.diff)))
      stop("the current implementation requires equally spaced values in the argument \"range.discrete\"\n")
    discrete.diff <- NULL
    if(is.R()) gc(verbose=FALSE)
  }    
  if(prior$nugget.prior != "fixed"){
    if(is.null(prior$nugget.discrete))
      stop("to include the nugget as random in the Bayesian analysis the argument nugget.discrete must be provided\n")
    if(prior$nugget.prior != "fixed")
      if(any(prior$nugget.discrete > 1))
        warning("when nugget is considered as random in the Bayesian analysis the values in the argument nugget.discrete must be relative nugget: tausq/sigmasq\n")
    discrete.diff <- diff(prior$nugget.discrete)
    if(round(max(1e08 * discrete.diff)) != round(min(1e08 * discrete.diff)))
      stop("the current implementation requires equally spaced values in the argument \"nugget.discrete\"\n")
    discrete.diff <- NULL
    if(is.R()) gc(verbose=FALSE)
  }
  ##
  ## checking output options
  ##
  if(!is.null(output$quantile.estimator)){
    if(is.numeric(output$quantile.estimator))
      if(any(output$quantile.estimator) < 0 | any(output$quantile.estimator) > 1)
        stop("quantiles indicators must be numbers in the interval [0,1]\n"
             )
    if(output$quantile.estimator == TRUE)
      output$quantile.estimator <- c(0.025000000000000001, 0.5, 0.97499999999999998)
  }
  if(!is.null(output$probability.estimator)){
    if(!is.numeric(output$probability.estimator))
      stop("probability.estimator must be a numeric value (or vector) of cut-off value(s)\n")
    if(length(output$probability.estimator)>1 & length(output$probability.estimator) != nrow(locations))
      stop("probability.estimator must either have length 1, or have length = nrow(locations)\n")
  }
  if(lambda != 1 & lambda != 0 & moments) {
    moments <- FALSE
    cat(paste("WARNING: moments cannot be computed for lambda =", 
              lambda, ". Argument moments was set to FALSE\n"))
  }
  if(lambda < 0 & output$mean.estimator) {
    output$mean.estimator <- FALSE
    cat("krige.bayes: mean.predictor set to F. The resulting distribution does not has expectation for $lambda < 0$\n"
        )
  }
  ##
  ## Box-Cox transformation of the data
  ##
  if(lambda != 1) {
    if(prior$range.prior == "fixed")
      stop("transformation option available only for the full Bayesian approach"
           )
    if(messages.screen)
      cat(paste("Box-Cox's transformation performed for $lambda=", lambda, "\n"))
    if(lambda == 0)
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  ##
  ## Checking the dimension of points to be predict (1 or 2)	
  ##
  if(all(locations != "no")) {
    if(is.vector(locations)) {
      if(length(locations) == 2) {
        locations <- t(as.matrix(locations))
        warning("THE FUNCTION IS CONSIDERING THAT YOU HAVE ENTERED WITH 1 POINT TO BE PREDICTED IN A TWO DIMENSION REGION\n")
      }
      else locations <- as.matrix(cbind(locations, 0))
    }
    else locations <- as.matrix(locations)
  }
  ##
  ## Building trend matrices:	
  ##
  if(inherits(model$trend.d, "formula") | inherits(model$trend.l, "formula")){
    if(any(locations != "no"))
      if((inherits(model$trend.d, "formula") == FALSE) | (inherits(model$trend.l, "formula") == FALSE))
        stop("model$trend.d and model$trend.l must have similar specification\n")
    if(messages.screen)
      cat("krige.bayes: Kriging with external trend to be performed using covariates provided by the user\n")
  }
  else{
    if (any(locations != "no") & (model$trend.d != model$trend.l)){
      stop("model$trend.l is different from model$trend.d")
    }
    if(messages.screen){
      if(model$trend.d == "cte")
        cat("krige.bayes: analysis assuming a constant mean\n")
      if(model$trend.d == "1st")
        cat("krige.bayes: analysis assuming a 1st degree polinomial trend\n")
      if(model$trend.d == "2nd") 
        cat("krige.bayes: analysis assuming a 2nd degree polinomial trend\n")
    }
  }
  trend.data <- trend.spatial(trend=model$trend.d, coords=coords)
  dimnames(coords) <- list(NULL, NULL)
  dimnames(trend.data) <- list(NULL, NULL)
  if(all(locations != "no")) {
    trend.l <- trend.spatial(trend=model$trend.l, coords=locations)
    dimnames(locations) <- list(NULL, NULL)
    dimnames(trend.l) <- list(NULL, NULL)
    if(is.matrix(trend.l))
      ni <- dim(trend.l)[1]
    else ni <- length(trend.l)
    ##    if(!is.null(smaller.locations)){
    ##      trend.smaller.l <- trend.spatial(trend=model$trend.l, coords=smaller.locations)
    ##    }
  }
  beta.size <- ncol(trend.data)
  ##
  ## Checking dimensions
  ##
  if(nrow(coords) != length(data)) stop(
           "number of data is different of number of data locations (coordinates)"
           )
  if(all(locations != "no")) {
    if(nrow(locations) != nrow(trend.l))
      stop("number of points to be estimated is different of the number of trend points"
           )
    dimnames(locations) <- list(NULL, NULL)
  }
  dimnames(coords) <- list(NULL, NULL)
  ##
  ## Anisotropy correction (must be placed here after trend matrices be defined)
  ##
  if(!is.null(model$aniso.pars)) {
    if(length(model$aniso.pars) != 2 | !is.numeric(model$aniso.pars))
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle and the anisotropy ratio\n")
    if(messages.screen)
      cat("krige.bayes: anisotropy parameters provided and assumed to be constants\n")
    coords <- coords.aniso(coords = coords, aniso.pars = model$aniso.pars)
    if(all(locations != "no"))
      locations <- coords.aniso(coords = locations, 
				aniso.pars = model$aniso.pars)
    data.dist <- as.vector(dist(coords))
  }
  ##
  n <- length(data)
  tausq <- nugget
  if(prior$sill.prior != "fixed") sigmasq <- 1
  else sigmasq <- sill
  ##
######################### PART 2 ##############################
  ## Prediction with fixed phi
  ##
  if(prior$range.prior == "fixed"){
    phi <- prior$range
    ## covariance matrix for data points
    invcov <- varcov.spatial(dists.lowertri=data.dist, cov.model = cov.model,
                             kappa = kappa, nugget = nugget,
                             cov.pars = c(sigmasq, phi), inv = TRUE,
                             scaled=TRUE, only.inv.lower.diag = TRUE)
    remove("data.dist")
    if(is.R()) gc(verbose=FALSE)
    ## Some matrices computations
###    ttiv <- crossprod(trend.data, invcov)
###    ttivtt <- ttiv %*% trend.data
###    ittivtt <- solve(ttivtt)
    ttivtt <- as.double(rep(0, beta.size * beta.size))
    ttivtt <- .C("bilinearform_XAY",
                 as.double(invcov$lower.inverse),
                 as.double(invcov$diag.inverse),
                 as.double(as.vector(trend.data)),
                 as.double(as.vector(trend.data)),
                 as.integer(beta.size),
                 as.integer(beta.size),
                 as.integer(n),
                 res=ttivtt)$res
    attr(ttivtt, "dim") <- c(beta.size, beta.size)
    ittivtt <- solve(ttivtt)
###    ttivz <- ttiv %*% data
    ttivz <- as.double(rep(0, beta.size))
    ttivz <- .C("bilinearform_XAY",
                as.double(invcov$lower.inverse),
                as.double(invcov$diag.inverse),
                as.double(as.vector(trend.data)),
                as.double(as.vector(data)),
                as.integer(beta.size),
                as.integer(1),
                as.integer(n),
                res=ttivz)$res
    ## covariance vector between data points and prediction locations
    v0mat <- as.double(rep(0, n*ni))
    .C("loccoords",
       as.double(as.vector(locations[,1])),
       as.double(as.vector(locations[,2])),
       as.double(as.vector(coords[,1])),
       as.double(as.vector(coords[,2])),
       as.integer(ni),
       as.integer(n),
       v0mat, DUP=FALSE)
    attr(v0mat, "dim") <- c(n, ni)
    v0mat <- cov.spatial(obj = v0mat,
                         cov.model = cov.model,
                         kappa = kappa, cov.pars = c(1, phi))
###    tv0iv <- t(apply(v0mat, 2, crossprod, y = invcov))
###    tv0ivz <- tv0iv %*% data
    tv0ivz <- as.double(rep(0,ni))
    tv0ivz <- .C("bilinearform_XAY",
                 as.double(invcov$lower.inverse),
                 as.double(invcov$diag.inverse),
                 as.double(as.vector(v0mat)),
                 as.double(data),
                 as.integer(ni),
                 as.integer(1),
                 as.integer(n),
                 res = tv0ivz)$res
###    tv0ivx <- apply(tv0iv, 1, crossprod, y = trend.data)
    tv0ivx <- as.double(rep(0, ni*beta.size))
    tv0ivx <- .C("bilinearform_XAY",
                 as.double(invcov$lower.inverse),
                 as.double(invcov$diag.inverse),
                 as.double(as.vector(v0mat)),
                 as.double(as.vector(trend.data)),
                 as.integer(ni),
                 as.integer(beta.size),
                 as.integer(n),
                 res=tv0ivx)$res
    attr(tv0ivx, "dim") <- c(ni, beta.size)
    tb <- trend.l - tv0ivx
                                        #    b <- t(tb)
###    tv0ivv0 <- diag(tv0iv %*% v0mat)
    tv0ivv0 <- as.double(rep(0,ni))
    tv0ivv0 <- .C("diag_quadraticform_XAX",
                  as.double(invcov$lower.inverse),
                  as.double(invcov$diag.inverse),
                  as.double(as.vector(v0mat)),
                  as.integer(ni),
                  as.integer(n),
                  res = tv0ivv0)$res
    v0mat <- NULL
    if(is.R()) gc(verbose=FALSE)
    ##
########################## PART 2.1 ##############################
    ## Simple kriging (fixed beta and sigmasq)
    ##
    if(prior$beta.prior == "fixed") {
      beta.mean <- paste("SK : beta provided by user: ",
                         beta)
      kb.results$predictive$mean <- as.vector(tv0ivz + as.vector(tb %*% beta))
      if(output$signal)
        kb.results$predictive$variance <- as.vector(sigmasq * (1 - tv0ivv0))
      else kb.results$predictive$variance <- as.vector(nugget + sigmasq * (1 - tv0ivv0))
      priors.messa <- c("no priors for parameters", 
                        "results corresponds to simple kriging")
      sill.mean <- paste("sill provided by user: ", sill)
    }
    ##
######################### PART 2.2 ##################################
    ## uncertainty in beta and sigmasq
    ##
    if(prior$beta.prior == "flat") {
      beta.mean <- ittivtt %*% ttivz
###      predict <- tv0ivz + b %*% beta.mean
      kb.results$predictive$mean <- as.vector(tv0ivz + as.vector(tb %*% beta.mean))
      ##      bi <- apply(tb, 2, crossprod, y = ittivtt)
      ##      if(ncol(trend.data) > 1)
      ##        bi <- t(bi)
      ##      bitb <- diag(bi %*% tb)
      if(beta.size == 1)
        bitb <- as.vector(tb^2) * as.vector(ittivtt)
      else{
        bitb <- as.double(rep(0,ni))
        bitb <- .C("diag_quadraticform_XAX",
                   as.double(ittivtt[lower.tri(ittivtt)]),
                   as.double(diag(ittivtt)),
                   as.double(as.vector(t(tb))),
                   as.integer(ni),
                   as.integer(beta.size),
                   res = bitb)$res
      }
      if(prior$sill.prior == "fixed") {
        ##
        ## 2.2.1 Ordinary kriging: uncertainty only in beta, flat prior
        ##
        if(output$signal)
          kb.results$predictive$variance <- as.vector(sigmasq * (1 - tv0ivv0 + bitb))
        else kb.results$predictive$variance <- as.vector(nugget + sigmasq * (1 - tv0ivv0 + bitb))
        priors.messa <- c(
                          "uninformative prior for beta", 
                          "no prior for sill", 
                          "no prior for range", 
                          "kriging estimates and variances equivalent to ordinary kriging"
                          )
        sill.mean <- paste("sill provided by user: ",
                           sill)
      }
      else {
        ##
        ## 2.2.2 using a conjugate normal prior
        ##
        residu <- as.vector(data - trend.data %*% beta.mean)
###        sill.mean <- (residu %*% invcov %*% residu)/
###          (n - length(beta.mean) - 2)
        resivres <- as.double(0)
        resivres <- .C("diag_quadraticform_XAX",
                       as.double(invcov$lower.inverse),
                       as.double(invcov$diag.inverse),
                       as.double(residu),
                       as.integer(1),
                       as.integer(n),
                       res = resivres)$res
        sill.mean <- resivres / (n - length(beta.mean) - 2)
        if(output$signal)
          kb.results$predictive$variance <- as.vector(sill.mean * (1 - tv0ivv0 + bitb))
        else kb.results$predictive$variance <- as.vector(nugget + sill.mean * (1 - tv0ivv0 + bitb))
        priors.messa <- c(
                          "uninformative prior for beta", 
                          "uninformative prior for sill", 
                          "kriging estimates (but not variances) equivalent to ordinary kriging"
                          )
      }
    }
    ##
    ## 2.2.3 Uncertainty only in the mean parameter (with normal prior)
    ##
    if(prior$beta.prior == "normal") {
      beta.var <- beta.var/sigmasq
      ibeta.var <- solve(beta.var)
      varbgz <- solve(ibeta.var + ttivtt)
      beta.mean <- varbgz %*% (ibeta.var %*% beta + ttivz)
      kb.results$predictive$mean <- as.vector(tv0ivz + as.vector(tb %*% beta.mean))
###      bi <- apply(tb, 2, crossprod, y = varbgz)
###      bitb <- diag(bi %*% tb)
      if(beta.size == 1)
        bitb <- as.vector(tb^2) * as.vector(ittivtt)
      else{
        bitb <- as.double(rep(0,ni))
        bitb <- .C("diag_quadraticform_XAX",
                   as.double(ittivtt[lower.tri(ittivtt)]),
                   as.double(diag(ittivtt)),
                   as.double(as.vector(t(tb))),
                   as.integer(ni),
                   as.integer(beta.size),
                   res = bitb)$res
      }
      if(output$signal)
        kb.results$predictive$variance <- as.vector(sigmasq * (1 - tv0ivv0) + bitb)
      else kb.results$predictive$variance <- as.vector(nugget + sigmasq * (1 - tv0ivv0) + bitb)
      priors.messa <- c("Normal prior for beta", 
                        "no prior for sill", "no prior for range")
      sill.mean <- paste("provided by user: ", sill)
    }
    ##
######################### PART 2.3 ##################################
    ## Preparing output
    ##
###    kb.results$predictive$mean <- as.vector(predict)
###    kb.results$predictive$variance <- as.vector(krige.var)
    kb.results$posterior$beta.mean <- beta.mean
    kb.results$posterior$sigmasq.mean <- sill.mean
    kb.results$posterior$phi.mean <- paste("provided by user: ", prior$range)
    kb.results$type.prediction <- priors.messa
  }
  else {
    ##	
######################### PART 3 ##############################
    ## considering the uncertainty in the parameter phi
    ##
    if(messages.screen)
      cat("krige.bayes: computing the discrete approximation of the posterior of phi/tausq.rel\n"
          )
    ##	
######################### PART 3.1 ##############################
    ## Computing the discrete posterior distribution for phi
    ##
    df.model <- as.vector(n - beta.size)
    phidist <- list()
    if(is.null(prior$range.discrete)){
      phidist$phi <- seq(0, max(data.dist), l = 51)
      warning("argument range.discrete not provided. Default values assumed\n")
    }
    else
      phidist$phi <- prior$range.discrete
    if(is.null(prior$nugget.discrete)){
      nugget.discrete <- nugget
      if(messages.screen)
        cat(paste("krige.bayes: nugget assumed to be fixed and equals to", nugget, "\n"))
    }
    else
      nugget.discrete <- prior$nugget.discrete
    n.range <- length(phidist$phi)
    n.nugget <- length(nugget.discrete)
    phi.val <- phidist$phi
    phidist$phi <- as.matrix(expand.grid(phidist$phi, nugget.discrete))
    dimnames(phidist$phi) <- list(NULL, NULL)
    krige.bayes.aux1 <- function(phinug){
      phi <- phinug[1]
      nugget <- phinug[2]
      if(round(1e12 * phi) == 0)
        covphi <- list(inverse = diag((1/(1 + nugget)), n), 
                       log.det.to.half = ((n/2) * log(1 + nugget)))
      else 
        covphi <- varcov.spatial(dists.lowertri = data.dist, cov.model = 
                                 cov.model, kappa = kappa, nugget = nugget,
                                 cov.pars = c(1, phi), inv = TRUE, det = TRUE)
      ttiv <- crossprod(trend.data, covphi$inverse)
      ttivtt <- ttiv %*% trend.data
      vbetaphi <- solve(ttivtt)
      eival.x <- eigen(ttivtt, symmetric = TRUE, only.values = TRUE)
      ttivttdet <- prod(eival.x$values)
      betaphi <- as.vector(vbetaphi %*% ttiv %*% data)
      res <- data - trend.data %*% betaphi
      sqres <- crossprod(res, covphi$inverse) %*% res
      df.model <- length(data) - length(betaphi)
      s2phi <- (1/df.model) * sqres
      logprobphi <- (-0.5) * log(ttivttdet) - (covphi$log.det.to.half) - (
                                                                          df.model/2) * log(s2phi)
      if(prior$range.prior == "reciprocal" & round(1e+08 * phi) != 0)
        logprobphi <- logprobphi - log(phi)
      if(prior$range.prior == "squared.reciprocal" & round(1e+08 * phi) != 0)
        logprobphi <- logprobphi - log(phi^2)
      if(prior$range.prior == "exponential")
        logprobphi <- logprobphi - exponential.prior.par * phi
      inv.lower <- covphi$inverse[lower.tri(covphi$inverse)]
      inv.diag <- diag(covphi$inverse)
      return(list(vbetaphi = vbetaphi, betaphi = betaphi, s2phi = s2phi,
                  logprobphi = logprobphi, inv.lower = inv.lower, inv.diag = inv.diag))
    }
    temp.res <- apply(phidist$phi, 1, krige.bayes.aux1)
    extract.f <- function(obj){return(obj$s2phi)}
    phidist$s2 <- as.vector(sapply(temp.res, extract.f))
    attr(phidist$s2, "dim") <- c(n.range, n.nugget)
    extract.f <- function(obj){return(obj$logprobphi)}
    phidist$logprobphi <- as.vector(sapply(temp.res, extract.f))
    attr(phidist$logprobphi, "dim") <- c(n.range, n.nugget)
    extract.f <- function(obj){return(obj$inv.lower)}
    inv.lower <-  as.vector(t(sapply(temp.res, extract.f)))
    attr(inv.lower, "dim") <- c(n.range, n.nugget, (n * (n - 1))/2)
    extract.f <- function(obj){return(obj$inv.diag)}
    inv.diag <-  as.vector(t(sapply(temp.res, extract.f)))
    attr(inv.diag, "dim") <- c(n.range, n.nugget, n)
    nframe.del <- sys.nframe()
    extract.f <- function(obj){return(obj$vbetaphi)}
    if(beta.size == 1){
      phidist$vbeta.vector <- as.vector(sapply(temp.res,extract.f))
      attr(phidist$vbeta.vector, "dim") <- c(n.range, n.nugget, 1)
    }
    else{
      phidist$vbeta.vector <- as.vector(t(sapply(temp.res, extract.f)))
##      dimnames(phidist$vbeta.vector) <- list(NULL, NULL)
      attr(phidist$vbeta.vector, "dim") <- c(n.range, n.nugget, (beta.size^2))
    }
    extract.f <- function(obj){return(obj$betaphi)}
    if(beta.size == 1){
      phidist$beta <- as.vector(sapply(temp.res, extract.f))
##      dimnames(phidist$beta) <- list(NULL, NULL)
      attr(phidist$beta, "dim") <- c(n.range, n.nugget, 1)
    }
    else{
      phidist$beta <- as.vector(t(sapply(temp.res, extract.f)))
      attr(phidist$beta, "dim") <- c(n.range, n.nugget, beta.size)
    }
    temp.res <- NULL
    if(is.R()) gc(verbose=FALSE)
    phidist$logprobphi <- phidist$logprobphi + abs(min(phidist$
                                                       logprobphi))
    phidist$probphi <- exp(phidist$logprobphi)
    sumprob <- sum(phidist$probphi)
    phidist$probphi <- phidist$probphi/sumprob
    ##
######################### PART 3.2 ##############################
    ## Sampling from the (parameters) posterior distribution
    ##
    if(messages.screen) cat(
         "krige.bayes: sampling from multivariate distribution of the parameters\n"
         )
    ind <- sample((1:(n.range * n.nugget)), n.posterior, replace = 
                  T, prob = as.vector(phidist$probphi))
    ind.unique <- sort(unique(ind))
    ind.length <- length(ind.unique)
    ind.table <- table(ind)
    phi.unique <- phidist$phi[ind.unique,  ]
    if(messages.screen) {
      cat("krige.bayes: samples and their frequencies from the distribution of  phi and tau.rel\n")
      print(rbind(phi = phi.unique[, 1], nugget.rel = 
                  phi.unique[, 2], frequency = ind.table))
      cat("\n")
    }
    phi.sam <- phidist$phi[ind,  ]
    vecpars.back.order <- order(ind)
    vec.s2 <- rep(as.vector(phidist$s2)[ind.unique], ind.table)
    samples.chisq <- rchisq(n.posterior, df = df.model)
    sigmasq <- (df.model * vec.s2)/samples.chisq
    if(beta.size == 1) {
      samples.beta <- rnorm(n.posterior, mean = 0, sd = 1)
      vec.beta <- rep(as.vector(phidist$beta)[ind.unique],
                      ind.table)
      vec.vbeta <- rep(as.vector(phidist$vbeta.vector)[
                                                       ind.unique], ind.table)
      beta <- vec.beta + sqrt(sigmasq * vec.vbeta) * 
        samples.beta
    }
    else {
      ind.beta <- matrix(phidist$beta, ncol = beta.size)[
                                         ind.unique,  ]
      ind.beta <- ind.beta[rep(1:ind.length, ind.table),
                           ]
      ind.vbeta <- matrix(phidist$vbeta.vector, ncol = 
                          beta.size^2)[ind.unique,  ]
      ind.vbeta <- ind.vbeta[rep(1:ind.length, ind.table),
                             ] * sigmasq
      ##      print("2.4: try to speed up this bit!")
      temp.res <- apply(ind.vbeta, 1, krige.bayes.aux3, 
                        beta.size = beta.size)
      beta <- ind.beta + t(temp.res)
      temp.res <- NULL
      if(is.R()) gc(verbose=FALSE)
    }
    if(beta.size == 1) {
      trend.mean <- mean(beta)
      trend.median <- median(beta)
    }
    else {
      trend.mean <- apply(beta, 2, mean)
      trend.median <- apply(beta, 2, median)
    }
    sill.mean <- mean(sigmasq)
    sill.median <- median(sigmasq)
    range.marg <- apply(phidist$probphi, 1, sum)
    range.marg <- range.marg/(sum(range.marg))
    range.mean <- phi.val %*% range.marg
    range.median <- median(phi.sam[, 1])
    ## check == here
    range.mode <- phi.val[range.marg == max(range.marg)]
    nugget.marg <- apply(phidist$probphi, 2, sum)
    nugget.marg <- nugget.marg/(sum(nugget.marg))
    nugget.mean <- nugget.discrete %*% nugget.marg
    nugget.median <- median(phi.sam[, 2])
    ## check == here
    nugget.mode <- nugget.discrete[nugget.marg == max(nugget.marg)]
    ##
    ## Computing the conditional (on phi and tausq.rel,
    ## the later if the case) modes for beta and sigmasq
    ##
    invcov.mode <- varcov.spatial(dists.lowertri=data.dist,
                                  cov.model = cov.model,
                                  kappa = kappa, nugget = nugget.mode,
                                  cov.pars = c(1, range.mode), inv = TRUE,
                                  only.inv.lower.diag = TRUE)
    remove("data.dist")
    if(is.R()) gc(verbose=FALSE)
    ttivtt.mode <- as.double(rep(0, beta.size * beta.size))
    ttivtt.mode <- .C("bilinearform_XAY",
                      as.double(invcov.mode$lower.inverse),
                      as.double(invcov.mode$diag.inverse),
                      as.double(as.vector(trend.data)),
                      as.double(as.vector(trend.data)),
                      as.integer(beta.size),
                      as.integer(beta.size),
                      as.integer(n),
                      res=ttivtt.mode)$res
    attr(ttivtt.mode, "dim") <- c(beta.size, beta.size)
    ittivtt.mode <- solve(ttivtt.mode)
    ttivz.mode <- as.double(rep(0, beta.size))
    ttivz.mode <- .C("bilinearform_XAY",
                     as.double(invcov.mode$lower.inverse),
                     as.double(invcov.mode$diag.inverse),
                     as.double(as.vector(trend.data)),
                     as.double(as.vector(data)),
                     as.integer(beta.size),
                     as.integer(1),
                     as.integer(n),
                     res=ttivz.mode)$res
    beta.mode.cond <-  ittivtt.mode %*% ttivz.mode
    resid.mode <- as.vector(data - trend.data %*% beta.mode.cond)
    sill.mode.cond <- (.C("diag_quadraticform_XAX",
                          as.double(invcov.mode$lower.inverse),
                          as.double(invcov.mode$diag.inverse),
                          as.double(as.vector(resid.mode)),
                          as.integer(1),
                          as.integer(n),
                          res = as.double(0.0))$res)/
                            (n - length(beta.mode.cond) + 2)
    if(beta.size == 1)
      kb.results$posterior$beta.summary <- c(mean = trend.mean, median = 
                                             trend.median, mode.cond = beta.mode.cond)
    else kb.results$posterior$beta.summary <-
      cbind(mean = trend.mean, median = trend.median,
            mode.cond = beta.mode.cond)
    kb.results$posterior$sigmasq.summary <-
      c(mean = sill.mean, median = sill.median, mode.cond = sill.mode.cond)
    kb.results$posterior$phi.summary <- c(mean = range.mean, median = 
                                          range.median, mode = range.mode)
    if(prior$nugget.prior != "fixed")
      kb.results$posterior$tausq.summary <- c(mean = nugget.mean,
                                              median = nugget.median, mode = 
                                              nugget.mode)
    else
      kb.results$posterior$tausq.summary <- paste("fixed tausq with value =", nugget)
    kb.results$posterior$beta.samples <- as.matrix(beta)[vecpars.back.order,  ]
    beta <- NULL
    if(is.R()) gc(verbose=FALSE)
    if(beta.size == 1)
      kb.results$posterior$beta.samples <- as.vector(kb.results$posterior$
                                                     beta.samples)
    kb.results$posterior$sigmasq.samples <- sigmasq[vecpars.back.order]
    sigmasq <- NULL
    if(is.R()) gc(verbose=FALSE)
    kb.results$posterior$phi.samples <- phi.sam[vecpars.back.order,1]
    if(prior$nugget.prior != "fixed")
      kb.results$posterior$tausq.samples <- phi.sam[vecpars.back.order,2]
    else
      kb.results$posterior$tausq.samples <- paste("fixed tausq with value =", nugget)
    phi.lev <- unique(phidist$phi[, 1])
    kb.results$posterior$phi.marginal <-
      data.frame(phi = phi.lev, expected = apply(phidist$probphi, 1, sum),
                 sampled = as.vector(table(factor(phi.sam[, 1],
                   levels = phi.lev)))/n.posterior)
    nug.lev <- unique(phidist$phi[, 2])
    if(prior$nugget.prior != "fixed")
      kb.results$posterior$nugget.marginal <-
        data.frame(nugget = nug.lev, expected = apply(phidist$probphi, 2, sum),
                   sampled = as.vector(table(factor(phi.sam[, 2],
                     levels = nug.lev)))/n.posterior)
    else
      kb.results$posterior$nugget.marginal <- paste("fixed tausq with value =", nugget)
    ##
######################### PART 3.3 ##############################
    ## Predictive distribution: sampling and computation
    ##
    if(all(locations != "no")) {
      if(messages.screen)
        cat("krige.bayes: prediction at locations provided\n")
      if(output$mean.estimator != FALSE | !is.null(output$probability.estimator) |
         !is.null(output$quantile.estimator))
        prediction.samples <- TRUE     
      message.prediction <- character()
      ni <- dim(trend.l)[1]
      if(is.null(n.predictive)) {
        include.it <- FALSE
        n.predictive <- n.posterior
        phi.sam <- phidist$phi[ind,  ]
        message.prediction <- c(message.prediction, "phi/tausq.rel samples for the predictive are same as for the posterior"
                                )
        if(messages.screen)
          cat("krige.bayes:", message.prediction, "\n")
      }
      else {
        include.it <- TRUE
        ind <- sample((1:(dim(phidist$phi)[1])), n.predictive,
                      replace = TRUE, prob = as.vector(phidist$probphi))
        ind.unique <- sort(unique(ind))
        ind.length <- length(ind.unique)
        ind.table <- table(ind)
        phi.unique <- phidist$phi[ind.unique,  ]
        message.prediction <- c(message.prediction, 
                                "phi/tausq.rel samples for the predictive are NOT the same as for the posterior ")
        if(messages.screen) {
          cat("krige.bayes:", message.prediction, "\n")
          cat("krige.bayes: samples and their frequencies from the distribution of  phi and tau.rel when drawing from the predictive distribution\n")
          print(rbind(phi = phi.unique[, 1], nugget.rel
                      = phi.unique[, 2], frequency = ind.table))
        }
        phi.sam <- phidist$phi[ind,  ]
        vecpars.back.order <- order(ind)
      }
###      d0mat <- apply(locations, 1, d0.krige)
      d0mat <- as.double(rep(0, ni*n))
      .C("loccoords",
         as.double(as.vector(locations[,1])),
         as.double(as.vector(locations[,2])),
         as.double(as.vector(coords[,1])),
         as.double(as.vector(coords[,2])),
         as.integer(ni),
         as.integer(n),
         d0mat, DUP=FALSE)
      attr(d0mat, "dim") <- c(n, ni)
      loc.coincide <- apply(d0mat, 2, function(x){any(x < 1e-10)})
      if(any(loc.coincide))
        loc.coincide <- (1:ni)[loc.coincide]
      else
        loc.coincide <- NULL
      if(!is.null(loc.coincide)){
        temp.f <- function(x, data){return(data[x < 1e-10])}
        data.coincide <- apply(d0mat[,loc.coincide, drop=FALSE],2,temp.f, data=data)
      }
      else
        data.coincide <- NULL
      if(is.R()) gc(verbose=FALSE)
      ##
      ## 3.3.1 Estimating the moments
      ##
      if(moments) {
        krige.bayes.aux10 <- function(phinug){
          counter <- get(".tempM.krige.bayes", pos=1)
          if(messages.screen)          
          krige.bayes.messages(moments = TRUE, n.disc = n.disc,
                                 .temp.ap = counter, ind.length=ind.length)
          phinug <- as.vector(phinug)
          phi.ind <- order(phi.val)[round(100000000. * phi.val) ==
                                    round(100000000. * phinug[1])]
          nug.ind <- order(nugget.discrete)[round(100000000. * nugget.discrete) ==
                                            round(100000000. * phinug[2])]
          v0 <- cov.spatial(obj = d0mat, cov.model = cov.model, kappa
                            = kappa, cov.pars = c(1, phinug[1]))
          tb <- as.double(rep(0, ni*beta.size))
          tb <- .C("bilinearform_XAY",
                   as.double(as.vector(inv.lower[phi.ind, nug.ind,])),
                   as.double(as.vector(inv.diag[phi.ind, nug.ind,])),
                   as.double(as.vector(v0)),
                   as.double(as.vector(trend.data)),
                   as.integer(ni),
                   as.integer(beta.size),
                   as.integer(n),
                   res=tb)$res
          attr(tb, "dim") <- c(ni, beta.size)
          tb <- trend.l - tb
          tv0ivdata <- as.double(rep(0,ni))
          tv0ivdata <- .C("bilinearform_XAY",
                          as.double(as.vector(inv.lower[phi.ind, nug.ind,])),
                          as.double(as.vector(inv.diag[phi.ind, nug.ind,])),
                          as.double(as.vector(v0)),
                          as.double(data),
                          as.integer(ni),
                          as.integer(1),
                          as.integer(n),
                          res = tv0ivdata)$res
          tmean <- tv0ivdata + tb %*% as.vector(phidist$
                                                beta[phi.ind, nug.ind,  ])
          tv0ivdata <- NULL
          if(is.R()) gc(verbose=FALSE)
          if(((round(1e12 * phinug[2]) == 0)) & (!is.null(loc.coincide)))
            tmean[loc.coincide] <- data.coincide
          tv0ivv0 <- as.double(rep(0,ni))
          tv0ivv0 <- .C("diag_quadraticform_XAX",
                        as.double(as.vector(inv.lower[phi.ind, nug.ind,])),
                        as.double(as.vector(inv.diag[phi.ind, nug.ind,])),
                        as.double(as.vector(v0)),
                        as.integer(ni),
                        as.integer(n),
                        res = tv0ivv0)$res
          v0 <- NULL
          if(is.R()) gc(verbose=FALSE)
          s2i <- phidist$s2[phi.ind, nug.ind]        
          vbetai <- matrix(as.vector(phidist$vbeta.vector[phi.ind, nug.ind,  ]),
                           ncol = beta.size, nrow = beta.size)
          if(beta.size == 1)
            tbivbb <- ((as.vector(tb))^2) * vbetai
          else{
            tbivbb <- as.double(rep(0,ni))
            tbivbb <- .C("diag_quadraticform_XAX",
                         as.double(vbetai[lower.tri(vbetai)]),
                         as.double(diag(vbetai)),
                         as.double(as.vector(t(tb))),
                         as.integer(ni),
                         as.integer(beta.size),
                         res = tbivbb)$res          
          }
          if(output$signal)
            tvar <- s2i * (1 - tv0ivv0 + tbivbb)
          else tvar <- s2i * ((1 + phinug[2]) - tv0ivv0 + tbivbb)
          tb  <- tbivbb <- tv0ivv0 <- NULL
          if(is.R()) gc(verbose=FALSE)
          if(((round(1e12 * phinug[2]) == 0) | output$signal) & (!is.null(loc.coincide)))
            tvar[loc.coincide] <- 0
          tvar[tvar < 1e-16] <- 0
          ## take care here, re-using object!
          tvar <- (df.model/(df.model - 2)) * tvar
          tvar <- tvar + (tmean)^2
          ##  here tvar is expec.y0.2 !!!
          assign(".tempM.krige.bayes", (.tempM.krige.bayes + 1), pos=1)
          return(cbind(tmean, tvar))
        }
        if(messages.screen)
          cat("krige.bayes: computing moments of the predictive distributions\n")
        n.disc <- dim(phidist$phi)[1]
        assign(".tempM.krige.bayes", 1, pos=1)
        temp.res <- apply(phidist$phi, 1, krige.bayes.aux10)
        rm(".tempM.krige.bayes", pos=1)
        temp.res <- temp.res %*% as.vector(phidist$probphi)
        kb.results$predictive$moments <-
          data.frame(expect.y0 = as.vector(temp.res[1:ni]),
                     expect.y0.2 = as.vector(temp.res[(ni + 1):(2 * ni)]))
        temp.res <- NULL
        if(is.R()) gc(verbose=FALSE)
        kb.results$predictive$moments$var.y0 <-
          kb.results$predictive$moments$expect.y0.2 -
           ((kb.results$predictive$moments$expect.y0)^2)
        if(lambda == 0) {
          temp <- kb.results$predictive$moments$expect.y0
          kb.results$predictive$moments$expect.y0 <-
            exp(temp + 0.5 * (kb.results$predictive$moments$var.y0))
          kb.results$predictive$moments$var.y0 <-
            (exp(2 * temp + kb.results$predictive$moments$var.y0)) *
              (exp(kb.results$predictive$moments$var.y0) - 1)
          temp <- NULL
          if(is.R()) gc(verbose=FALSE)
          kb.results$predictive$moments$expect.y0.2 <-
            kb.results$predictive$moments$var.y0 +
              (kb.results$predictive$moments$expect.y0^2)
        }
      }
      ##
      ## 3.3.2 Sampling from the predictive
      ##
      if(output$simulations.predictive){
        names(ind.table) <- NULL
        krige.bayes.aux20 <- function(phinug){
          counter <- get(".tempS.krige.bayes", pos=1)
          if(messages.screen)
            krige.bayes.messages(moments = FALSE, n.disc = n.disc,
                                 .temp.ap = counter, ind.length=ind.length)
          phinug <- as.vector(phinug)
          phi.ind <- order(phi.val)[round(100000000. * phi.val) ==
                                    round(100000000. * phinug[1])]
          nug.ind <- order(nugget.discrete)[round(100000000. * nugget.discrete) ==
                                            round(100000000. * phinug[2])]
          v0 <- cov.spatial(obj = d0mat, cov.model = cov.model, kappa
                            = kappa, cov.pars = c(1, phinug[1]))
          tb <- as.double(rep(0, ni*beta.size))
          tb <- .C("bilinearform_XAY",
                   as.double(as.vector(inv.lower[phi.ind, nug.ind,])),
                   as.double(as.vector(inv.diag[phi.ind, nug.ind,])),
                   as.double(as.vector(v0)),
                   as.double(as.vector(trend.data)),
                   as.integer(ni),
                   as.integer(beta.size),
                   as.integer(n),
                   res=tb)$res
          attr(tb, "dim") <- c(ni, beta.size)
          tb <- trend.l - tb
          tv0ivdata <- as.double(rep(0,ni))
          tv0ivdata <- .C("bilinearform_XAY",
                          as.double(as.vector(inv.lower[phi.ind, nug.ind,])),
                          as.double(as.vector(inv.diag[phi.ind, nug.ind,])),
                          as.double(as.vector(v0)),
                          as.double(data),
                          as.integer(ni),
                          as.integer(1),
                          as.integer(n),
                          res = tv0ivdata)$res
          tmean <- tv0ivdata + tb %*% as.vector(phidist$
                                                beta[phi.ind, nug.ind,  ])
          tv0ivdata <- NULL
          if(is.R()) gc(verbose=FALSE)
          if(((round(1e12 * phinug[2]) == 0)) & (!is.null(loc.coincide)))
            tmean[loc.coincide] <- data.coincide
          s2i <- phidist$s2[phi.ind, nug.ind]        
          vbetai <- matrix(as.vector(phidist$vbeta.vector[phi.ind, nug.ind,  ]),
                           ncol = beta.size, nrow = beta.size)
          if(((round(1e12 * phinug[2]) == 0) | output$signal) & (!is.null(loc.coincide))){
            v0 <- v0[,-(loc.coincide)]
            nloc <- ni - length(loc.coincide)
            tmean.coincide <- tmean[loc.coincide]
            tmean <- tmean[-(loc.coincide)]
            tb <- tb[-(loc.coincide),]
          }
          else nloc <- ni
          Nsims <- ind.table[.tempS.krige.bayes]
          normalsc <- rnorm(nloc*Nsims)
          chisc <- rchisq(Nsims, df=df.model)
          sqglchi <- sqrt(df.model/chisc)
          if (output$signal) Dval <- 1.0 else Dval <-  1.0 + phinug[2]
          if(beta.size == 1){
            Blower <- 0
            Bdiag <- vbetai
          }
          else{
            Blower <- vbetai[lower.tri(vbetai)]
            Bdiag <- diag(vbetai)
          }
          R0 <- as.double(rep(0.0, (nloc*(nloc+1))/2))
          if(cov.model.number > 10){
            if(((round(1e12 * phinug[2]) == 0) | output$signal) & (!is.null(loc.coincide))){
              .C("distdiag",
                 as.double(locations[-(loc.coincide),1]),
                 as.double(locations[-(loc.coincide),2]),
                 as.integer(nloc),
                 R0, DUP = FALSE)
            }
            else
              .C("distdiag",
                 as.double(locations[,1]),
                 as.double(locations[,2]),
                 as.integer(ni),
                 R0, DUP = FALSE)
            R0 <- cov.spatial(R0, cov.pars=c(1, phinug[1]), cov.model=cov.model, kappa=kappa)
          }
          else{
            if(((round(1e12 * phinug[2]) == 0) | output$signal) & (!is.null(loc.coincide))){
              .C("cor_diag",
                 as.double(locations[-(loc.coincide),1]),
                 as.double(locations[-(loc.coincide),2]),
                 as.integer(nloc),
                 as.integer(cov.model.number),
                 as.double(phinug[1]),
                 as.double(kappa),
                 R0, DUP = FALSE)
            }
            else
              .C("cor_diag",
                 as.double(locations[,1]),
                 as.double(locations[,2]),
                 as.integer(ni),
                 as.integer(cov.model.number),
                 as.double(phinug[1]),
                 as.double(kappa),
                 R0, DUP = FALSE)
          }  
          normalsc <- .C("kb_sim",
                         as.double(tmean),
                         out = as.double(normalsc),
                         as.double(as.vector(inv.lower[phi.ind, nug.ind,])),
                         as.double(as.vector(inv.diag[phi.ind, nug.ind,])),
                         as.double(as.vector(v0)),
                         as.integer(nloc),
                         as.integer(n),
                         as.double(Dval),
                         as.integer(Nsims),
                         as.double(sqglchi),                      
                         as.double(s2i),                      
                         as.double(Blower),
                         as.double(Bdiag),
                         as.double(as.vector(t(tb))),
                         as.integer(beta.size),
                         as.double(R0))$out
          attr(normalsc, "dim") <- c(nloc, Nsims)
          v0 <- R0 <- tb <- NULL
          if(is.R()) gc(verbose=FALSE)
          ##
          assign(".tempS.krige.bayes", (.tempS.krige.bayes + 1), pos=1)
          if(((round(1e12 * phinug[2]) == 0) | output$signal) & (!is.null(loc.coincide))){
            result <- matrix(0, nrow=ni, ncol=Nsims)
            result[-(loc.coincide),] <- normalsc
            result[loc.coincide,] <- rep(tmean.coincide, Nsims)
            return(result)
          }
          else
            return(normalsc)
        }
        assign(".tempS.krige.bayes", 1, pos=1)
        kb.results$predictive$simulations <-
          matrix(unlist(apply(phi.unique, 1, krige.bayes.aux20)),ncol=n.predictive)
        remove("inv.lower", envir = sys.frame(nframe.del))
        remove("inv.diag", envir = sys.frame(nframe.del))
        remove(".tempS.krige.bayes", pos=1)
        if(is.R()) gc(verbose=FALSE)
        if(messages.screen)
          cat("krige.bayes: preparing output of predictive distribution\n")
        ##
        ## Back transforming
        ##
        if(lambda != 1) {
          cat("krige.bayes: Data transformation (Box-Cox) performed. Results returned in the original scale\n"
              )
          if(lambda == 0)
            kb.results$predictive$simulations <-
              exp(kb.results$predictive$simulations)
          else {
            if(lambda > 0)
              kb.results$predictive$simulations[kb.results$predictive$simulations < (-1/lambda)] <- -1/lambda
            if(lambda < 0)
              kb.results$predictive$simulations[kb.results$predictive$simulations > (-1/lambda)] <- -1/lambda
            kb.results$predictive$simulations <-
              ((kb.results$predictive$simulations * lambda) + 1)^(1/lambda)
          }
        }
        ##
        ## 3.3.3 mean estimators
        ##
        if(output$mean.estimator) {
          kb.results$predictive$mean.simulations <- as.vector(apply(kb.results$predictive$simulations, 1, mean))
          kb.results$predictive$variance.simulations <- as.vector(apply(kb.results$predictive$simulations, 1, var))
          message.prediction <- c(message.prediction, "mean at each location in $predictive$mean and predicted variances in $predictive$variances")
        }
        ##
        ## 3.3.4 quantile estimators
        ##
        if(!is.null(output$quantile.estimator)) {
          kb.results$predictive$quantiles <-
            apply(kb.results$predictive$simulations, 1, quantile, probs
                  = output$quantile.estimator)
          if(length(output$quantile.estimator) > 1) {
            kb.results$predictive$quantiles <- as.data.frame(t(kb.results$predictive$quantiles))
            qname <- rep(0, length(output$quantile.estimator))
            for(i in 1:length(output$quantile.estimator))
              qname[i] <- paste("q", 100 * output$quantile.estimator[
                                                              i], sep = "")
            names(kb.results$predictive$quantiles) <- qname
          }
          else {
            kb.results$predictive$quantiles <- as.vector(kb.results$predictive$
                                                         quantiles)
          }
          message.prediction <- c(message.prediction, 
                                  "Predicted quantile(s) at each location returned in $predictive$quantiles")
        }
        ##
        ## 3.3.5 probability estimators
        ##
        if(!is.null(output$probability.estimator)) {
          kb.results$predictive$probability <-
            apply(kb.results$predictive$simulations, 1, krige.bayes.aux2, cutoff = output$probability.estimator)
          if(length(output$probability.estimator) > 1){
            kb.results$predictive$probability <- as.data.frame(t(kb.results$predictive$probability))
          pname <- rep(0, length(output$probability.estimator))
          for(i in 1:length(output$probability.estimator))
            pname[i] <- paste("cutoff", output$probability.estimator[i], sep = "")
            names(kb.results$predictive$probability) <- pname
          }
          else
          kb.results$predictive$probability <- as.vector(kb.results$predictive$
                                                         probability)  
          message.prediction <- c(message.prediction,
                                  "Estimated probabilities of being less than the provided cutoff(s), at each location, returned in $predictive$probability")
        }
        ##
        ## samples from  predictive
        ##
        if(output$keep.simulations)
          kb.results$predictive$simulations <-
            kb.results$predictive$simulations[, vecpars.back.order]
        else kb.results$predictive$simulations <- NULL
        if(is.R()) gc(verbose=FALSE)
        if(include.it){
          phi.lev <- unique(phidist$phi[, 1])
          kb.results$predictive$phi.marginal <-
            data.frame(phi = phi.lev, expected = apply(phidist$probphi, 1, sum),
                       sampled = as.vector(table(factor(phi.sam[, 1],
                         levels = phi.lev)))/n.predictive)
          nug.lev <- unique(phidist$phi[, 2])
          if(prior$nugget.prior != "fixed")
            data.frame(nugget = nug.lev, expected = apply(phidist$probphi, 2, sum),
                       sampled = as.vector(table(factor(phi.sam[, 2],
                         levels = nug.lev)))/n.predictive)
          else
            kb.results$predictive$nugget.marginal <- paste("fixed tausq with value =", nugget)
          kb.results$predictive$nugget.marginal <-
            data.frame(nugget = nug.lev, expected = apply(phidist$probphi, 2, sum),
                       sampled = as.vector(table(factor(phi.sam[, 2],
                         levels = nug.lev)))/n.predictive)
        }
      }
      kb.results$predictive$type.prediction <-
        "Spatial prediction performed taking into account uncertainty in trend (mean), sill and range parameters"
      kb.results$message.prediction <- message.prediction
    }
    else {
      kb.results$predictive <- "no locations to perform prediction were provided"
      kb.results$message.prediction <- 
        "Only Bayesian estimation of model parameters were computed"
      if(messages.screen)
        print(kb.results$message.prediction)
    }
  }
  if(all(nugget != 0) & prior$nugget.prior == "fixed")
    kb.results$nugget.fixed <- nugget
  kb.results$.Random.seed <- seed
#  if(info.for.prediction == TRUE & prior$range.prior != "fixed")
#    kb.results$info.for.prediction <-
#      list(coords = coords, cov.model = 
#           cov.model, kappa = kappa, trend.d = trend.d,
#           data = data, beta.size = beta.size, n = n,
#           phidist = phidist, phi.val = phi.val, nugget = 
#           nugget, distribution = distribution, 
#           n.posterior = n.posterior, ind.length = 
#           ind.length, vecpars.back.order = 
#           vecpars.back.order, phi.unique = phi.unique,
#           ind.table = ind.table, df.model = df.model, output$mean.estimator
#           = output$mean.estimator, output$quantile.estimator = 
#           output$quantile.estimator, output$probability.estimator = 
#           output$probability.estimator, ind = ind, lambda = 
#           lambda, moments = moments)
  kb.results$call <- call.fc
  class(kb.results) <- c("krige.bayes", "kriging")
  if(messages.screen)
    cat("krige.bayes: done!\n")
  return(kb.results)
}

"krige.bayes.aux2" <- 
function(x, cutoff)
{
	# auxiliary function to perform calculation for the function krige.bayes
	ncut <- length(cutoff)
	lx <- length(x)
	if(ncut > 1) {
		result <- rep(0, ncut)
		for(i in 1:ncut)
			result[i] <- sum(x <= cutoff[i])/lx
	}
	else result <- sum(x <= cutoff)/lx
	return(result)
      }

"krige.bayes.aux3" <- 
  ##
  ## This function produces a sample from  a multivariate normal distribution 
  ## mean is 0 and cov.values is a vector of length beta.size^2
  function(cov.values, beta.size)
{
                                        #  dm <- length(mean)
                                        #  dc <- dim(cov.values)[1]
                                        #  if(dm != dc)
                                        #    stop("mean and cov.values must have compatible dimensions")
  cov.values <- matrix(cov.values, ncol = beta.size)
  cov.svd <- svd(cov.values)
  cov.decomp <- cov.svd$u %*% (t(cov.svd$u) * sqrt(cov.svd$d))
  zsim <- as.vector(cov.decomp %*% rnorm(beta.size))
  return(zsim)
}
"lines.krige.bayes" <- 
  function(obj, max.dist, length = 100, summary.posterior = c("mode", "median", "mean"), ...)
{
  spost <- match.arg(summary.posterior)
  if(is.null(obj$call$cov.model))
    cov.model <- "exponential"
  else {
    cov.model <- obj$call$cov.model
    if(obj$call$cov.model == "matern" | obj$call$cov.model == "powered.exponential" | obj$
       call$cov.model == "cauchy" | obj$call$cov.model == "gneiting-matern")
      kappa <- obj$call$kappa
    else kappa <- NULL
  }
  distance <- seq(0, max.dist, length = length)
  if(spost == "mode")
    spost1 <- "mode.cond"
  else spost1 <- spost
  cov.pars <- c(obj$posterior$sigmasq.summary[spost1], obj$posterior$phi.summary[spost])
  names(cov.pars) <- NULL
  if(is.numeric(obj$posterior$tausq.summary))
    nugget <- obj$posterior$tausq.summary[spost] * cov.pars[1]
  else nugget <- 0
  names(nugget) <- NULL
  sill.total <- nugget + cov.pars[1]
  gamma <- (nugget + cov.pars[1]) - cov.spatial(distance, cov.model = cov.model, kappa = kappa, 
                                                cov.pars = cov.pars)
  lines(distance, gamma, ...)
}

"image.krige.bayes" <-
  function (obj, locations,
            values.to.plot = c("moments.mean", "moments.variance",
              "mean.simulations", "variance.simulations",
              "quantiles", "probability", "simulation"),
            number.col, coords.data, ...) 
{
  if(is.null(locations))
    stop("prediction locations must be provided")
  locations <- locations[order(locations[, 2], locations[,1]), ]
  x <- as.numeric(levels(as.factor(locations[, 1])))
  nx <- length(x)
  y <- as.numeric(levels(as.factor(locations[, 2])))
  ny <- length(y)
  coords.lims <- apply(locations,2,range)
  coords.diff <- diff(coords.lims)
  if(coords.diff[1] != coords.diff[2]){
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[,ind.min] <- coords.lims[,ind.min] + c(-coords.diff.diff, coords.diff.diff)/2
  }
  par(pty = "s")
  if (is.numeric(values.to.plot)){
    image(x, y, matrix(values.to.plot, ncol=ny), xlim= coords.lims[,1], ylim=coords.lims[,2],...)
  }
  else{
    values <- match.arg(values.to.plot)
    switch(values,
           moments.mean =
           {
             values <- matrix(obj$predictive$moments$expect.y0, ncol = ny);
             cat("plotting map the mean of the predictive distribution\n")
           },
           moments.variance = {
             values <- matrix(obj$predictive$moments$var.y0, ncol = ny);
             cat("plotting map the variance of the predictive distribution\n")
           },
           mean.simulations=
           {
             values <- matrix(obj$predictive$mean.simulations, ncol = ny);
             cat("plotting map the mean of the simulations from the predictive distribution\n")
           },
           variance.simulations =
           {
             values <- matrix(obj$predictive$variance.simulations, ncol = ny);
             cat("plotting map the variance of the predictive distribution\n")
           },
           quantile =
           {
             if(!is.vector(obj$predictive$quantiles))
               if(missing(number.col))
                 stop("argument number.col must be provided")
               else
                 values <- matrix(obj$predictive$quantiles[,number.col], ncol = ny)
             else
               values <- matrix(obj$predictive$quantiles, ncol = ny);
             cat("plotting map a quantile of the predictive distribution\n")
           },
           probability =
           {
             if(!is.vector(obj$predictive$probability))
               if(missing(number.col))
                 stop("argument number.col must be provided")
               else
                 values <- matrix(obj$predictive$probability[,number.col], ncol = ny)
             else
               values <- matrix(obj$predictive$probability, ncol = ny);
             cat("plotting map a simulation of the predictive distribution\n")
           },
           simulation =
           {
             values <- matrix(obj$predictive$simulations[,number.col], ncol = ny);
             cat("plotting map the variance of the predictive distribution\n")
           },
           stop("wrong specification for values to plot")
           )
    image(x, y, values, xlim= coords.lims[,1], ylim=coords.lims[,2],...)
  }
  if(!missing(coords.data))
    points(coords.data)
  return(invisible())
}

"persp.krige.bayes" <-
  function (obj, locations,
            values.to.plot = c("moments.mean", "moments.variance",
              "mean.simulations", "variance.simulations",
              "quantiles", "probability", "simulation"), number.col, ...) 
{
  if(is.null(locations))
    stop("prediction locations must be provided")
  locations <- locations[order(locations[, 2], locations[,1]), ]
  x <- as.numeric(levels(as.factor(locations[, 1])))
  nx <- length(x)
  y <- as.numeric(levels(as.factor(locations[, 2])))
  ny <- length(y)
  if (nx == ny) 
    par(pty = "s")
  if (is.numeric(values.to.plot)){
    persp(x,y,matrix(values.to.plot, ncol=ny), ...)
  }
  else{
    values <- match.arg(values.to.plot)
    switch(values,
           moments.mean =
           {
             values <- matrix(obj$predictive$moments$expect.y0, ncol = ny);
             cat("plotting map the mean of the predictive distribution\n")
           },
           moments.variance = {
             values <- matrix(obj$predictive$moments$var.y0, ncol = ny);
             cat("plotting map the variance of the predictive distribution\n")
           },
           mean.simulations=
           {
             values <- matrix(obj$predictive$mean.simulations, ncol = ny);
             cat("plotting map the mean of the simulations from the predictive distribution\n")
           },
           variance.simulations =
           {
             values <- matrix(obj$predictive$variance.simulations, ncol = ny);
             cat("plotting map the variance of the predictive distribution\n")
           },
           quantile =
           {
             if(!is.vector(obj$predictive$quantiles))
               if(missing(number.col))
                 stop("argument number.col must be provided")
               else
                 values <- matrix(obj$predictive$quantiles[,number.col], ncol = ny)
             else
               values <- matrix(obj$predictive$quantiles, ncol = ny);
             cat("plotting map a quantile of the predictive distribution\n")
           },
           probability =
           {
             if(!is.vector(obj$predictive$probability))
               if(missing(number.col))
                 stop("argument number.col must be provided")
               else
                 values <- matrix(obj$predictive$probability[,number.col], ncol = ny)
             else
               values <- matrix(obj$predictive$probability, ncol = ny);
             cat("plotting map a simulation of the predictive distribution\n")
           },
           simulation =
           {
             values <- matrix(obj$predictive$simulations[,number.col], ncol = ny);
             cat("plotting map the variance of the predictive distribution\n")
           },
           stop("wrong specification for values to plot")
           )
    persp(x, y, values, ...)
  }
  return(invisible())
}

"model.control" <-
  function(trend.d = "cte", trend.l = "cte",
           cov.model = "matern",
           kappa=0.5, aniso.pars=NULL, lambda=1) 
{
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  return(list(trend.d = trend.d, trend.l = trend.l,
              cov.model = cov.model,
              kappa=kappa, aniso.pars=aniso.pars, lambda=lambda))
}

"prior.control" <-
  function(beta.prior = c("flat", "normal", "fixed"),
           beta = NULL, beta.var = NULL,
           sill.prior = c("reciprocal", "fixed"), sill = NULL, 
           range.prior = c("uniform", "exponential", "fixed",
             "squared.reciprocal","reciprocal"), exponential.prior.par = 1,
           range = NULL, range.discrete = NULL, 
           nugget.prior = c("fixed", "uniform"), nugget = 0,
           nugget.discrete = NULL)
{
  beta.prior <- match.arg(beta.prior)
  sill.prior <- match.arg(sill.prior)
  range.prior <- match.arg(range.prior)
  nugget.prior <- match.arg(nugget.prior)
  return(list(beta.prior = beta.prior, beta = beta, beta.var = beta.var,
              sill.prior = sill.prior, sill = sill, 
              range.prior = range.prior,
              exponential.prior.par = exponential.prior.par,
              range = range, range.discrete = range.discrete, 
              nugget.prior = nugget.prior, nugget = nugget,
              nugget.discrete = nugget.discrete))
}

"output.control" <-
  function(n.posterior = 1000, n.predictive = NULL,
           simulations.predictive = TRUE, keep.simulations = TRUE,
           mean.estimator = TRUE, quantile.estimator = NULL,
           probability.estimator = NULL, messages.screen = TRUE,
           signal = TRUE, moments = TRUE)
{
  return(list(n.posterior = n.posterior, n.predictive = n.predictive,
              moments = moments,
              simulations.predictive = simulations.predictive,
              keep.simulations = keep.simulations,
              mean.estimator = mean.estimator,
              quantile.estimator = quantile.estimator,
              probability.estimator = probability.estimator,
              signal = signal, messages.screen = messages.screen))
}

"krige.bayes.messages" <- 
function(moments, n.disc, .temp.ap, ind.length)
{
  
  if(moments){
    if(n.disc <= 30)
      cat(paste("computing moments: point",
                .temp.ap, "out of", n.disc, "\n"))
    if(n.disc > 30 & n.disc <= 300)
      if(.temp.ap %% 10 == 1){
        cat(paste("computing moments: point",
                  .temp.ap, "out of", n.disc, "\n"))
        if(.temp.ap == 1) cat("                   (counting for each 10)\n")
      }
    if(n.disc > 300)
      if(.temp.ap %% 100 == 1){
        cat(paste("computing moments: point",
                  .temp.ap, "out of", n.disc, "\n"))
        if(.temp.ap == 1) cat("                   (counting for each 100)\n")
      }
  }
  else{
    if(ind.length <= 30)
      cat(paste("simulating from the predictive for the parameter set", 
                .temp.ap, "out of", ind.length, "\n"))
    if(ind.length > 30 & ind.length <= 300)
      if(.temp.ap %% 10 == 1){
        cat(paste("simulating from the predictive for the parameter set",
                  .temp.ap, "out of", 
                  ind.length, "\n"))
        if(.temp.ap == 1) cat("                   (counting for each 10)\n")
      }
    if(ind.length > 300)
      if(.temp.ap %% 100 == 1){
        cat(paste("simulating from the predictive for the parameter set",
                  .temp.ap, "out of", 
                  ind.length, "\n"))
        if(.temp.ap == 1) cat("                   (counting for each 100)\n")
      }
  }
}




