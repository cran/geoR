"likfit" <-
  function (geodata, coords=geodata$coords, data=geodata$data,
            trend = "cte", ini.cov.pars,
            fix.nugget = FALSE, nugget = 0, 
            fix.kappa = TRUE, kappa = 0.5, 
            fix.lambda = TRUE, lambda = 1, 
            fix.psiA = TRUE, psiA = 0, 
            fix.psiR = TRUE, psiR = 1, 
            cov.model = c("matern", "exponential", "gaussian",
              "spherical", "circular", "cubic", "wave",
              "powered.exponential", "cauchy", "gneiting",
              "gneiting.matern", "pure.nugget"),
            method = "ML",
            components = FALSE, nospatial = TRUE,
            limits = likfit.limits(), messages.screen = TRUE, ...) 
{
  ##
  ## Checking input
  ##
  if(is.R()) require(mva)
  call.fc <- match.call()
  temp.list <- list()
  ##
  cov.model <- match.arg(cov.model)
  temp.list$cov.model <- cov.model
  ##
  if (method == "REML" | method == "reml" | method == "rml") 
    method <- "RML"
  if(method == "ML" | method == "ml")
    method <- "ML"
  temp.list$method <- method
  if(is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)){
    ini.cov.pars <- as.matrix(ini.cov.pars)
    if(nrow(ini.cov.pars) == 1)
      ini.cov.pars <- as.vector(ini.cov.pars)
    else{
      if((cov.model != "pure.nugget") & (ncol(ini.cov.pars) != 2))
        stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq and phi")
    }
  }
  if(is.vector(ini.cov.pars)){
    if((cov.model != "pure.nugget") & (length(ini.cov.pars) != 2))
      stop("\nini.cov.pars must be a vector with 2 components: \ninitial values for sigmasq and phi")
  }
  if(fix.kappa & !is.null(kappa))
    if(cov.model == "matern" & kappa == 0.5)
      cov.model <- "exponential"
  coords <- temp.list$coords <- as.matrix(coords)
  n <- temp.list$n <- length(data)
  if ((2*n) != length(coords))
    stop("\nnumber of locations does not match with number of data")
  temp.list$xmat <- trend.spatial(trend=trend, coords=coords)
  beta.size <- temp.list$beta.size <- dim(temp.list$xmat)[2]
  ##
  ## Checking for multiple initial values for preliminar search of   
  ## best initial value
  ##
  if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1) | (length(lambda) > 1) | (length(psiR) > 1) | (length(psiA) > 1)){
    if(messages.screen)
      cat("likfit: searching for best initial value ...")
    .likGRF.dists.vec <<- as.vector(dist(coords))
    temp.list$z <- as.vector(data)
    temp.list$trend <- trend
    ini.temp <- matrix(ini.cov.pars, ncol=2)
    grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]), phi=unique(ini.temp[,2]), tausq=unique(nugget), kappa=unique(kappa), lambda=unique(lambda), psiR=unique(psiR), psiA=unique(psiA)))
    temp.f <- function(parms, temp.list){
      return(loglik.GRF(coords=temp.list$coords, data=temp.list$z, cov.model=temp.list$cov.model, cov.pars=parms[1:2], nugget=parms["tausq"], kappa=parms["kappa"], lambda=parms["lambda"], psiR=parms["psiR"], psiA=parms["psiA"], trend= temp.list$trend, method=temp.list$method, compute.dists=F))
    }
    grid.lik <- apply(grid.ini, 1, temp.f, temp.list=temp.list)
    ini.temp <- grid.ini[which(grid.lik == max(grid.lik)),]
    if(messages.screen){
      cat(" selected values:\n")
      print(round(as.vector(ini.temp), dig=3))
    }
    names(ini.temp) <- NULL
    ini.cov.pars <- ini.temp[1:2]
    nugget <- ini.temp[3]
    kappa <- ini.temp[4]
    lambda <- ini.temp[5]
    psiR <- ini.temp[6]
    psiA <- ini.temp[7]
    grid.ini <- NULL
    temp.list$trend <- NULL
    if(is.R()) {remove(".likGRF.dists.vec", pos=1); gc(verbose=FALSE)}    
    else remove(".likGRF.dists.vec", where=1)    
  }
  ##
  tausq <- nugget
  ##
  ## Box-Cox transformation for fixed lambda
  ##
  if(fix.lambda) {
    if(round(lambda, dig=4) == 1) {
      temp.list$log.jacobian <- 0
      temp.list$z <- as.vector(data)
    }
    else {
      if(any(data <= 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      Jdata <- data^(lambda - 1)
      if(any(Jdata <= 0))
        temp.list$log.jacobian <- log(prod(Jdata))
      else temp.list$log.jacobian <- sum(log(Jdata))
      Jdata <- NULL
      if(is.R()) gc(verbose=FALSE)
      if(round(lambda, dig=4) == 0)
        temp.list$z <- log(data)
      else temp.list$z <- ((data^lambda) - 1)/lambda
    }
  }
  else{
    temp.list$z <- as.vector(data)
    temp.list$log.jacobian <- NULL
  }
  ##
  ## Coordinates transformation for fixed anisotropy parameters
  ##
  if(fix.psiR & fix.psiA){
    if(psiR != 1 | psiA != 0)
      coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    if(is.R()) assign(".likGRF.dists.vec", as.vector(dist(coords)), pos=1)
    else assign(".likGRF.dists.vec", as.vector(dist(coords)), where=1)
    range.dist <- range(.likGRF.dists.vec)
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
  }
  ##
  ##
  ##
  ini <- ini.cov.pars[2]
  ##  fixed.pars <- NULL
  lower.optim <- c(limits$phi["lower"])
  upper.optim <- c(limits$phi["upper"])
  fixed.values <- list()
  if(fix.nugget) {
    ##    fixed.pars <- c(fixed.pars, 0)
    fixed.values$tausq <- nugget
  }
  else {
    ini <- c(ini, nugget/ini.cov.pars[1])
    lower.optim <- c(lower.optim, limits$tausq.rel["lower"])
    upper.optim <- c(upper.optim, limits$tausq.rel["upper"])
  }
  if(fix.kappa){
##    fixed.kappa <- c(fixed.pars, kappa)
    fixed.values$kappa <- kappa
  }
  else {
    ini <- c(ini, kappa)
    lower.optim <- c(lower.optim, limits$kappa["lower"])
    upper.optim <- c(upper.optim, limits$kappa["upper"])
  }
  if(fix.lambda){
##    fixed.pars <- c(fixed.pars, lambda)
    fixed.values$lambda <- lambda
  }
  else {
    ini <- c(ini, lambda)
    lower.optim <- c(lower.optim, limits$lambda["lower"])
    upper.optim <- c(upper.optim, limits$lambda["upper"])
  }
  if(fix.psiR){
##    fixed.pars <- c(fixed.pars, psiR)
    fixed.values$psiR <- psiR
  }
  else {
    ini <- c(ini, psiR)
    lower.optim <- c(lower.optim, limits$psiR["lower"])
    upper.optim <- c(upper.optim, limits$psiR["upper"])
  }
  if(fix.psiA){
##    fixed.pars <- c(fixed.pars, psiA)
    fixed.values$psiA <- psiA
  }
  else {
    ini <- c(ini, psiA)
    lower.optim <- c(lower.optim, limits$psiA["lower"])
    upper.optim <- c(upper.optim, limits$psiA["upper"])
  }
  ## This must be here, after the previous ones:
  if(fix.nugget & nugget > 0){
    ## Warning: Inverting order here, ini will be now: c(phi,sigmasg)
    ini <- c(ini, ini.cov.pars[1])
    lower.optim <- c(lower.optim, limits$sigmasq["lower"])
    upper.optim <- c(upper.optim, limits$sigmasq["upper"])
##    fixed.pars <- c(fixed.pars, ini.cov.pars[1])
##    fixed.values$sigmasq <- 0
  }
  names(ini) <- NULL
  ##
  ip <- list(f.tausq = fix.nugget, f.kappa = fix.kappa, f.lambda = fix.lambda,
                  f.psiR = fix.psiR, f.psiA = fix.psiA)
  ##  
  if (messages.screen == TRUE) {
    cat("-----------------------------------------------------------------------\n")
    cat("likfit: Initialising likelihood maximisation using the function ")
    if(is.R()) cat("optim.\n") else cat("nlminb.\n")
    cat("likfit: Use control() to pass arguments for the maximisation function.")
    cat("\n        For more details see documentation for ")
    if(is.R()) cat("optim.\n") else cat("nlminb.\n")        
    cat("likfit: It is highly advisable to run this function several\n        times with different initial values for the parameters.\n")
    cat("likfit: WARNING: This step can be time demanding!\n")
    cat("-----------------------------------------------------------------------\n")
  }
  npars <- beta.size + 2 + sum(unlist(ip)==FALSE)
  if(is.R()){
    lik.minim <- optim(par = ini, fn = negloglik.GRF, method="L-BFGS-B",
                       lower=lower.optim, upper=upper.optim,
                       fp=fixed.values, ip=ip, temp.list = temp.list, ...)
  }
  else{
    lik.minim <- nlminb(ini, negloglik.GRF,
                        lower=lower.optim, upper=upper.optim,
                        fp=fixed.values, ip=ip, temp.list = temp.list, ...)
  }
  ##
  if(messages.screen == TRUE) 
    cat("likfit: end of numerical maximisation.\n")
  par.est <- lik.minim$par
  if(any(par.est < 0)) par.est <- round(par.est, dig=12)
  phi <- par.est[1]
  ##
  ## Values of the maximised likelihood
  ##
  if(is.R())
    value.min <- lik.minim$value
  else
    value.min <- lik.minim$objective
  if(method == "ML"){
    if(ip$f.tausq & (tausq > 0))
      loglik.max <-  (- value.min) + (n/2)*(-log(2*pi))
    else
      loglik.max <-  (- value.min) + (n/2)*(-log(2*pi) + log(n) -1)
  }
  if(method == "RML"){
    xx.eigen <- eigen(crossprod(temp.list$xmat), symmetric = TRUE, only.values = TRUE)
    if(ip$f.tausq & (tausq > 0))
      loglik.max <- (- value.min) - ((n-beta.size)/2)*(log(2*pi)) +
        0.5 * sum(log(xx.eigen$values))
    else
      loglik.max <- (- value.min) - ((n-beta.size)/2)*(log(2*pi)) +
        ((n-beta.size)/2)*(log(n-beta.size)) - ((n-beta.size)/2) +
          0.5 * sum(log(xx.eigen$values))
  }
  ##
  ## Assigning values for estimated parameters
  ##
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    psiA <- par.est[2]
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    psiR <- par.est[2]
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    psiR <- par.est[2]
    psiA <- par.est[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    lambda  <- par.est[2]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    lambda  <- par.est[2]
    psiA <- par.est[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    lambda  <- par.est[2]
    psiR <- par.est[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    lambda  <- par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    kappa  <-  par.est[2]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    kappa  <-  par.est[2]
    psiA <- par.est[3]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    kappa  <-  par.est[2]
    psiR <- par.est[3]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    kappa  <-  par.est[2]
    psiR <- par.est[3]
    psiA <- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
    psiA <- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
    psiR<- par.est[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    kappa <-  par.est[2]
    lambda <- par.est[3]
    psiR<- par.est[4]
    psiA<- par.est[5]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    psiA<- par.est[3]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    psiR<- par.est[3]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    psiR<- par.est[3]
    psiA<- par.est[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiA <- par.est[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    lambda <- par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    psiA <- par.est[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    psiR <- par.est[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    psiR <- par.est[4]
    psiA <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
    psiA <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
    psiR <- par.est[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- par.est[2]
    kappa <-  par.est[3]
    lambda <- par.est[4]
    psiR <- par.est[5]
    psiA <- par.est[6]
  }
  ##
  if(fix.nugget & nugget > 0){
    sigmasq <- par.est[length(par.est)]
    if(sigmasq > 1e-12) tausq <- nugget/sigmasq
    check.sigmasq <- TRUE
  }
  else check.sigmasq <- FALSE
  ##
  ##
  ## Transforming data acccording to the estimated lambda (Box-Cox) parameter
  ##
  if(!fix.lambda) {
    if(round(lambda, dig=4) == 1) {
      log.jacobian.max <- 0
      temp.list$z <- data
    }
    else {
      if(any(data^(lambda - 1) <= 0))
        log.jacobian.max <- log(prod(data^(lambda - 1)))
      else log.jacobian.max <- sum(log(data^(lambda - 1)))
      temp.list$z <- ((data^lambda)-1)/lambda
    }
  }
  else{
    log.jacobian.max <- temp.list$log.jacobian
  }
  ##
  ## Transforming coords for estimated anisotropy (if the case)
  ##
  if(fix.psiR & fix.psiA){
    if(is.R()) remove(".likGRF.dists.vec", pos=1)
    else{
      remove(".likGRF.dists.vec", where=1)
    }
  }
  else{
    if(psiR != 1 | psiA != 0)
      coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    .likGRF.dists.vec <- as.vector(dist(coords))
    range.dist <- range(.likGRF.dists.vec)
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
    remove(".likGRF.dists.vec")
  }      
  if(is.R()) gc(verbose=FALSE)
  ##
  ## Computing estimated beta and tausq/sigmasq (if the case)
  ##
  if((phi < 1e-12))
    siv <- diag(x=1/sqrt((1+tausq)), n)
  else{
    if(check.sigmasq){
      if(sigmasq < 1e-12){
        if(!fix.nugget)
          siv <- diag(x=1/sqrt((1+tausq)), n)
        else
          siv <- diag(x=1/sqrt((tausq)), n)          
      }
      else
        siv <- varcov.spatial(coords = coords, cov.model = cov.model,
                              kappa = kappa,
                              nugget = tausq, cov.pars = c(1, phi),
                              inv=TRUE, sqrt.inv = TRUE,
                              det = FALSE)$sqrt.inverse
    }
    else
      siv <- varcov.spatial(coords = coords, cov.model = cov.model,
                            kappa = kappa,
                            nugget = tausq, cov.pars = c(1, phi),
                            inv=TRUE, sqrt.inv = TRUE,
                            det = FALSE)$sqrt.inverse
  }
  sivx <- crossprod(siv, temp.list$xmat)
  xivx <- crossprod(sivx)
  sivy <- crossprod(siv, temp.list$z)
  xivy <- crossprod(sivx, sivy)
  betahat <- solve(xivx, xivy)
  res <- as.vector(temp.list$z - temp.list$xmat %*% betahat)
  if(!fix.nugget | (round(1e+12 * nugget) == 0)){
    res <- as.vector(temp.list$z - temp.list$xmat %*% betahat)
    ssres <- as.vector(crossprod(crossprod(siv,res)))
    if(method == "ML")
      sigmasq <- ssres/n
    else
      sigmasq <- ssres/(n - beta.size)
  }
  if(fix.nugget){
    if(nugget > 0)
      tausq <- nugget
  }
  else tausq <- tausq * sigmasq
  betahat.var <- solve(xivx)
  if(sigmasq > 1e-12) betahat.var <- sigmasq * betahat.var
#  if(!fix.nugget & phi < 1e-16){
#    tausq <- sigmasq + tausq
#    sigmasq <- 0
#  }
  n.model.pars <- beta.size+7
  par.su <- data.frame(status=rep(-9,n.model.pars))
  ind.par.su <- c(rep(0, beta.size), ip$f.tausq, 0, 0, ip$f.kappa,
                  ip$f.psiR, ip$f.psiA,ip$f.lambda)
  par.su$status <- ifelse(ind.par.su,"fixed", "estimated")
  par.su$values <- round(c(betahat, tausq, sigmasq, phi, kappa, psiR, psiA, lambda), dig=4)
  if(beta.size == 1) beta.name <- "beta"
  else beta.name <- paste("beta", 0:(beta.size-1), sep="")
  row.names(par.su) <- c(beta.name, "tausq", "sigmasq", "phi", "kappa",
                             "psiR", "psiA", "lambda")
  par.su <- par.su[c((1:(n.model.pars-3)), n.model.pars-1, n.model.pars-2, n.model.pars),] 
  ##
  ## Preparing output
  ##
  lik.results <- list(cov.model = cov.model,
                      nugget = tausq,
                      cov.pars=c(sigmasq, phi),
                      kappa = kappa,
                      beta = as.vector(betahat),
                      beta.var = betahat.var,
                      lambda = lambda,
                      aniso.pars = c(psiA = psiA, psiR = psiR),
                      method = method,
                      loglik = loglik.max,
                      npars = npars,
                      AIC = (loglik.max - npars),
                      BIC = (loglik.max - 0.5 * log(n) * npars),
                      parameters.summary = par.su,
                      info.minimisation.function = lik.minim,
                      max.dist = max.dist,
                      trend.matrix= temp.list$xmat,
                      transform.info = list(fix.lambda = fix.lambda,
                        log.jacobian = log.jacobian.max))
  ##
  ## Likelihood results for the model without spatial correlation
  ##
  if(nospatial){
    if(fix.lambda){
      beta.ns <- solve(crossprod(temp.list$xmat), crossprod(temp.list$xmat, temp.list$z))
      ss.ns <- sum((as.vector(temp.list$z - temp.list$xmat %*% beta.ns))^2)
      if(method == "ML"){
        nugget.ns <- ss.ns/n
        loglik.ns <- (n/2)*((-log(2*pi)) - log(nugget.ns) - 1) + temp.list$log.jacobian
      }
      if(method == "RML"){
        nugget.ns <- ss.ns/(n-beta.size)
        loglik.ns <- ((n-beta.size)/2)*((-log(2*pi)) - log(nugget.ns) -1) + temp.list$log.jacobian
      }
      npars.ns <- beta.size + 1 + fix.lambda
      lambda.ns <- lambda
    }
    else{
      bc.list <- list(n = n, beta.size = beta.size,
                      data = data, xmat = temp.list$xmat,
                      method = method)
      if(is.R())
        lik.lambda.ns <- optim(par=1, fn = boxcox.ns, method="L-BFGS-B",
                               lower=limits$lambda["lower"], upper=limits$lambda["upper"], bc.list=bc.list)
      else
        lik.lambda.ns <- nlminb(par=1, fn = boxcox.ns,
                                lower=limits$lambda["lower"], upper=limits$lambda["upper"], data=data,
                                bc.list = bc.list)
      bc.list <- NULL
      if(is.R()) gc(verbose=FALSE)
      lambda.ns <- lik.lambda.ns$par
      tdata.ns <- ((data^lambda.ns)-1)/lambda.ns
      beta.ns <- solve(crossprod(temp.list$xmat),crossprod(temp.list$xmat,tdata.ns))
      ss.ns <- sum((as.vector(tdata.ns - temp.list$xmat %*% beta.ns))^2)
      if(is.R())
        value.min.ns <- lik.lambda.ns$value
      else
        value.min.ns <- lik.lambda.ns$objective
      if(method == "ML"){
        loglik.ns <- (- value.min.ns)+ (n/2)*((-log(2*pi)) + log(n) - 1)
        nugget.ns <- ss.ns/n
      }
      if(method == "RML"){
        nugget.ns <- ss.ns/(n-beta.size)
        loglik.ns <- (- value.min.ns)+ ((n-beta.size)/2)*((-log(2*pi)) +
                                                          log(n-beta.size) - 1)
      }      
      npars.ns <- beta.size + 1 + fix.lambda
    }
    lik.results$nospatial <- list(beta.ns = beta.ns, variance.ns = nugget.ns,
                                  loglik.ns = loglik.ns, npars.ns = npars.ns,
                                  lambda.ns = lambda.ns)
  }
  ##
  ## Computing residuals and predicted values
  ## (isolated components of the model)
  ##
  if (components) {
    if(!fix.psiR & !fix.psiA)
      if(psiR != 1 | psiA != 0)
        coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    trend.comp <- temp.list$z - res
    invcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                             kappa = kappa, nugget = tausq,
                             cov.pars = c(sigmasq, phi), inv=TRUE)$inverse 
    covmat.signal <- varcov.spatial(coords = coords, cov.model = cov.model, 
                                    kappa = kappa, nugget = 0,
                                    cov.pars = c(sigmasq, phi))$varcov
    spatial.comp <- as.vector(covmat.signal %*% invcov %*% res)
    predict.comp <- trend.comp + spatial.comp
    residual.comp <- as.vector(temp.list$z - predict.comp)
    residual.std <- as.vector(invcov %*% residual.comp)
    residual.trend.std <- as.vector(invcov %*% res)
    s2.random <- (crossprod(res,invcov) %*% res)/(n - beta.size)
    s2 <- (crossprod(residual.comp,invcov) %*% residual.comp)/(n - beta.size)
  }
  if (length(lik.results$beta.var) == 1)
    lik.results$beta.var <- as.vector(lik.results$beta.var)
  if (length(lik.results$beta) > 1){
    if(inherits(trend, "formula"))
      beta.names <- c("intercept", paste("covar", 1:(ncol(temp.list$xmat)-1), sep = ""))
    else
      if (trend == "1st")
        beta.names <- c("intercept", "x", "y")
      else
        if (trend == "2nd")
          beta.names <- c("intercept", "x", "y", "x2", "xy", "y2")
    names(lik.results$beta) <- beta.names
  }
  if (components) {
    lik.results$model.components <- data.frame(trend = trend.comp, spatial = spatial.comp, residuals = residual.comp)
    lik.results$s2 <- s2
    lik.results$s2.random <- s2.random
  }
  lik.results$call <- call.fc
  ##
  ## Assigning classes
  ##
  if(is.R())
    class(lik.results) <- c("likGRF", "variomodel")
  else{
    if(version$major <= 4)
      class(lik.results) <- c("likGRF", "variomodel")
    else oldClass(lik.results) <- c("likGRF", "variomodel")
  }
  ##
  ## Some warning messages about particular possible results
  ##
  if(messages.screen){
    if((lik.results$cov.pars[1] < (0.01 * (lik.results$nugget + lik.results$cov.pars[1])))& lik.results$cov.pars[2] > 0)
      cat("\nWARNING: estimated sill is less than 1 hundredth of the total variance. Consider re-examine the model excluding spatial dependence\n" )      
    if((lik.results$cov.pars[2] > (10 * max.dist)) & lik.results$cov.pars[1] > 0 )
      cat("\nWARNING: estimated range is more than 10 times bigger than the biggest distance between two points. Consider re-examine the model:\n 1) excluding spatial dependence if estimated sill is too low and/or \n 2) taking trends (covariates) into account\n" ) 
    if(((lik.results$cov.pars[2] < (0.1 * min.dist)) & (lik.results$cov.pars[1] > 0)) & lik.results$cov.pars[2] > 0)
      cat("\nWARNING: estimated range is less than 1 tenth of the minimum distance between two points. Consider re-examine the model excluding spatial dependence\n" ) 

  }
  ##
  return(lik.results)
}

"negloglik.GRF" <-
  function(pars, fp, ip, temp.list)
### pars : values for the parameters to be estimated
  ## sequence is c(phi, tausq, kappa, lambda, psiR, psiA, sigmasq)
### fixed pars: parameters considered fixed
### ind.pars : list indicating which are fixed and which are to be estimated
  ##
  ## Warning:
  ##  if fix.nugget = TRUE and nugget > 0 ,
  ## sigmasq should be passed and fp$nugget is the value of the nugget
  ## otherwise the RELATIVE nugget should be passed
{
  n <- temp.list$n
  p <- temp.list$beta.size
  log.jacobian <- temp.list$log.jacobian
  ## Obligatory parameter:
  phi <- pars[1]
  ## Others
  if(ip$f.tausq){
    if(fp$tausq > 0){
      npars.min <- length(pars)
      sigmasq <- pars[npars.min]
    }
    else sigmasq <- 1
  }
  else sigmasq <- 1
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[2]
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[2]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[2]
    psiA <- pars[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- fp$psiR
    psiA <- pars[3]
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- pars[3]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- fp$kappa
    lambda <- pars[2]
    psiR <- pars[3]
    psiA <- pars[4]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[3]
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- pars[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- pars[4]
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- fp$psiA
  }
  if(ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- fp$tausq
    kappa <- pars[2]
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- pars[5]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[3]
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- fp$lambda
    psiR <- pars[3]
    psiA <- pars[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- fp$psiR
    psiA <- pars[4]
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- fp$kappa
    lambda <- pars[3]
    psiR <- pars[4]
    psiA <- pars[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- fp$psiR
    psiA <- pars[4]
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- pars[4]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- fp$lambda
    psiR <- pars[4]
    psiA <- pars[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- fp$psiR
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- fp$psiR
    psiA <- pars[5]
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- pars[5]
    psiA <- fp$psiA
  }
  if(!ip$f.tausq & !ip$f.kappa & !ip$f.lambda & !ip$f.psiR & !ip$f.psiA){
    tausq <- pars[2]
    kappa <- pars[3]
    lambda <- pars[4]
    psiR <- pars[5]
    psiA <- pars[6]
  }
  ##
  ## Absurd values
  ##
  if(kappa < 1e-04) return(1e+32)
  if(round(1e+16*(tausq+sigmasq)) == 0) return(1e+32)
  ##
  ## Anisotropy
  ##
  if(!ip$f.psiR | !ip$f.psiA){
    coords.c <- coords.aniso(temp.list$coords, aniso.pars=c(psiA, psiR))
    if(is.R()) assign(".likGRF.dists.vec", as.vector(dist(coords.c)), pos=1)
    else assign(".likGRF.dists.vec", as.vector(dist(coords.c)), where=1)
  }
  ##
  ## Box-Cox transformation
  ##
  if(!ip$f.lambda){
    if(round(lambda, dig=4) == 1) {
      log.jacobian <- 0
    }
    else {
      if(any(temp.list$z <= 0))
        stop("Transformation not allowed for zero or negative data")
      data <- temp.list$z^(lambda - 1)
      if(any(data <= 0)) log.jacobian <- log(prod(data))
      else log.jacobian <- sum(log(data))
      data <- NULL
    }
    if(round(lambda, dig=4) == 0)
      data <- log(temp.list$z)
    else data <- ((temp.list$z^lambda) - 1)/lambda
  }
  else data <- temp.list$z
  ##
  ## Computing likelihood
  ##
  ## NOTE: Likelihood for Independent observations 
  ##       arbitrary criteria used here:
  ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
  ##
  if((phi < 1e-16) | (sigmasq < 1e-16)){
    if(ip$f.tausq)
      iv <- list(sqrt.inverse = diag(x=1/sqrt((tausq+sigmasq)), n),
                 log.det.to.half = (n/2) * log(tausq+sigmasq))
    else
      iv <- list(sqrt.inverse = diag(x=1/sqrt((1+tausq)), n),
                 log.det.to.half = (n/2) * log(1+tausq))
  }
  else{
    iv <- varcov.spatial(dists.lowertri = .likGRF.dists.vec,
                         cov.model = temp.list$cov.model, kappa=kappa,
                         nugget = tausq, cov.pars=c(sigmasq, phi),
                         sqrt.inv = TRUE, det = TRUE)
  }
  sivx <- crossprod(iv$sqrt.inverse, temp.list$xmat)
  xivx <- crossprod(sivx)
  sivy <- crossprod(iv$sqrt.inverse, data)
  xivy <- crossprod(sivx, sivy)  
  betahat <- solve(xivx, xivy)
  res <- data - temp.list$xmat %*% betahat
  ssres <- as.vector(crossprod(crossprod(iv$sqrt.inverse,res)))
  if(temp.list$method == "ML"){
    if(ip$f.tausq & (tausq > 0))
      negloglik <- iv$log.det.to.half +  0.5 * ssres - log.jacobian
    else
      negloglik <- (n/2) * log(ssres) +  iv$log.det.to.half - log.jacobian
  }
  if(temp.list$method == "RML"){
    if(length(as.vector(xivx)) == 1) {
      choldet <- 0.5 * log(xivx)
    }
    else {
      chol.xivx <- chol(xivx)
      choldet <- sum(log(diag(chol.xivx)))
    }
    if(ip$f.tausq & (tausq > 0))
      negloglik <- iv$log.det.to.half +  0.5 * ssres + choldet - log.jacobian
    else
      negloglik <- ((n-p)/2) * log(ssres) +  iv$log.det.to.half +
        choldet - log.jacobian
  }
  if(negloglik > 1e+32) negloglik <- 1e32
  return(negloglik)
}

"boxcox.ns" <- function(lambda, bc.list)
{
  data <- bc.list$data
  n <- bc.list$n
  ##
  if(round(lambda, dig=4) == 1) {
    log.jacobian <- 0
    y <- data
  }
  else {
    if(any(data <= 0))
      stop("Transformation option not allowed when there are zeros or negative data")
    Jdata <- data^(lambda - 1)
    if(any(Jdata <= 0))
      log.jacobian <- log(prod(Jdata))
    else log.jacobian <- sum(log(Jdata))
    Jdata <- NULL
    if(is.R()) gc(verbose=FALSE)
    if(round(lambda, dig=4) == 0)
      y <- log(data)
    else y <- ((data^lambda) - 1)/lambda
  }
  beta.ns <- solve(crossprod(bc.list$xmat), crossprod(bc.list$xmat, y))
  ss.ns <- sum((as.vector(y) - as.vector(bc.list$xmat %*% beta.ns))^2)
  if(bc.list$method == "ML")
    neglik <- (n/2) * log(ss.ns) - log.jacobian
  if(bc.list$method == "RML")
    neglik <- ((n-bc.list$beta.size)/2) * log(ss.ns) - log.jacobian
  ##
  return(as.vector(neglik))
}

"likfit.limits" <-
  function(phi = c(lower=0, upper=+Inf),
           sigmasq = c(lower=0, upper=+Inf),
           tausq.rel = c(lower=0, upper=+Inf),
           kappa = c(lower=0, upper=+Inf),
           lambda = c(lower=-3, upper=3),
           psiR = c(lower=1, upper=+Inf),
           psiA = c(lower=0, upper=2*pi)
           )
{
  if(length(phi) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(length(sigmasq) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(length(tausq.rel) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(length(kappa) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(length(lambda) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi")
  if(length(psiR) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(length(psiA) != 2)
    stop("phi must be a 2 components vector with lower and upper limits for the parameter phi") 
  if(phi[1] >= phi[2])
    stop("parameter phi: lower limit greater or equal upper limit")
  if(sigmasq[1] >= sigmasq[2])
    stop("parameter sigmasq: lower limit greater or equal upper limit")
  if(tausq.rel[1] >= tausq.rel[2])
    stop("parameter tausq.rel: lower limit greater or equal upper limit")
  if(kappa[1] >= kappa[2])
    stop("parameter kappa: lower limit greater or equal upper limit")
  if(lambda[1] >= lambda[2])
    stop("parameter lambda: lower limit greater or equal upper limit")
  if(psiR[1] >= psiR[2])
    stop("parameter psiR: lower limit greater or equal upper limit")
  if(psiA[1] >= psiA[2])
    stop("parameter psiA: lower limit greater or equal upper limit")
  names(phi) <- c("lower", "upper")
  names(sigmasq) <- c("lower", "upper")
  names(tausq.rel) <- c("lower", "upper")
  names(kappa) <- c("lower", "upper")
  names(lambda) <- c("lower", "upper")
  names(psiR) <- c("lower", "upper")
  names(psiA) <- c("lower", "upper")
  return(list(phi = phi, sigmasq = sigmasq,
              tausq.rel = tausq.rel, kappa = kappa,
              lambda = lambda, psiR = psiR, psiA = psiA))
}

"print.likGRF" <-
  function(obj, digits = "default", ...)
{
  if(is.R() & digits == "default") digits <- max(3, getOption("digits") - 3)
  else digits <- options()$digits
  est.pars <- as.vector(obj$parameters.summary[obj$parameters.summary[,1] == "estimated",2])
  names.est.pars <- dimnames(obj$parameters.summary[obj$parameters.summary[,1] == "estimated",])[[1]]
  names(est.pars) <- names.est.pars
  cat("likfit: estimated model parameters:\n")
  print(round(est.pars, digits=digits))
  cat("\nlikfit: maximised log-likelihood = ")
  cat(round(obj$loglik, digits=digits))
  cat("\n")
  return(invisible())
}  

"summary.likGRF" <-
  function(obj, ...)
{
  names.pars <- dimnames(obj$parameters.summary)[[1]]
  summ.lik <- list()
  if(obj$method == "ML")
    summ.lik$method <- "maximum likelihood"
  if(obj$method == "RML")
    summ.lik$method <- "restricted maximum likelihood"
  summ.lik$mean.component <- obj$beta
  names(summ.lik$mean.component) <- names.pars[1:length(obj$beta)]
  summ.lik$cov.model <- obj$cov.model
  summ.lik$spatial.component <- obj$parameters.summary[c("sigmasq", "phi"),]
  summ.lik$spatial.component.extra <- obj$parameters.summary[c("kappa", "psiA", "psiR"),]
  summ.lik$nugget.component <- obj$parameters.summary[c("tausq"),, drop=FALSE]
  summ.lik$transformation  <- obj$parameters.summary[c("lambda"),, drop=FALSE]
  summ.lik$likelihood <- c(log.L = obj$loglik, n.params = as.integer(obj$npars),
                               AIC = obj$AIC, BIC = obj$BIC)
  summ.lik$estimated.pars <- dimnames(obj$parameters.summary[obj$parameters.summary[,1] == "estimated",])[[1]]
  likelihood.info <- c(log.L = obj$loglik, n.params = as.integer(obj$npars),
                       AIC = obj$AIC, BIC = obj$BIC)
  summ.lik$call <- obj$call
  class(summ.lik) <- "summary.likGRF"
  return(summ.lik)
}

"print.summary.likGRF" <-
  function(obj, digits = "default", ...)
{
  if(is.R() & digits == "default") digits <- max(3, getOption("digits") - 3)
  else digits <- options()$digits
  cat("Summary of the parameter estimation\n")
  cat("-----------------------------------\n")
  cat(paste("Estimation method:", obj$method, "\n"))
  cat("\n")
  ##
  ## Estimates of the model components
  ## Model: Y(x) = X\beta + S(x) + e 
  ##
  cat("Parameters of the mean component (trend):")
  cat("\n")
  print(round(obj$mean.component, digits=digits))
  cat("\n")
  ##
  cat("Parameters of the spatial component:")
  cat("\n")
  cat(paste("   correlation function:", obj$cov.model))
  cat(paste("\n      (estimated) variance parameter sigmasq (partial sill) = ", round(obj$spatial.component[1,2], dig=digits)))
  cat(paste("\n      (estimated) cor. fct. parameter phi (range parameter)  = ", round(obj$spatial.component[2,2], dig=digits)))
  if(obj$cov.model == "matern" | obj$cov.model == "powered.exponential" |
     obj$cov.model == "cauchy" | obj$cov.model == "gneiting.matern"){
    kappa <- obj$spatial.component.extra["kappa",2]
    if(obj$spatial.component.extra["kappa",1] == "estimated")
      cat(paste("\n      (estimated) extra parameter kappa =", round(kappa, digits=digits)))
    else{
      cat(paste("\n      (fixed) extra parameter kappa = ", kappa))
      if(obj$cov.model == "matern" & round((1e12 *kappa)  == 0.5))
      cat(" (exponential)")
    }
  }
  cat("\n")
  ##
  aniso <-  obj$spatial.component.extra[c("psiA", "psiR"),]
  psiApsiR <- obj$spatial.component.extra[c("psiA", "psiR"),2]
  cat("   anisotropy parameters:")
  if(aniso["psiA",1] == "estimated")
    cat(paste("\n      (estimated) anisotropy angle =",
              round(psiApsiR[1], digits=digits),
              " (",round((psiApsiR[1]*360)/(2*pi), dig=1), "degrees )"))
  else
    cat(paste("\n      (fixed) anisotropy angle =", psiApsiR[1],
              " (",(psiApsiR[1]*360)/(2*pi), "degrees )"))
  if(aniso["psiR",1] == "estimated")
    cat(paste("\n      (estimated) anisotropy ratio =",
              round(psiApsiR[2], digits=digits)))
  else
    cat(paste("\n      (fixed) anisotropy ratio =", psiApsiR[2]))
  cat("\n")
  cat("\n")  
  cat("Parameter of the error component:")
  if(obj$nugget.component[,1] == "estimated")
    cat(paste("\n      (estimated) nugget = ", round(obj$nugget.component[,2], dig=digits)))
  else
    cat(paste("\n      (fixed) nugget =", obj$nugget.component[,2]))
  cat("\n")
  cat("\n")
  cat("Transformation parameter:")
  cat("\n")
  lambda <- obj$transformation[,2]
  if(obj$transformation[,1] == "estimated")
    cat(paste("      (estimated) Box-Cox parameter =", round(lambda, dig=digits)))
  else{
    cat(paste("      (fixed) Box-Cox parameter =", lambda))
    if(lambda == 1) cat(" (no transformation)")
    if(lambda == 0) cat(" (log-transformation)")
  }
  cat("\n")
  cat("\n")
  cat("Maximised Likelihood:")
  cat("\n")
  print(round(obj$likelihood, digits=digits))
  cat("\n")
  cat("Call:")
  cat("\n")
  print(obj$call)
  cat("\n")
  invisible(obj)
}

"loglik.GRF" <-
  function(geodata, coords=geodata$coords, data=geodata$data, cov.model="exp", cov.pars, nugget=0, kappa=0.5, lambda=1, psiR=1, psiA=0, trend="cte", method="ML", compute.dists=T)
{
  if (method == "REML" | method == "reml" | method == "rml") 
    method <- "RML"
  if(method == "ML" | method == "ml")
    method <- "ML"
  n <- nrow(coords)
  xmat <- trend.spatial(trend=trend, coords=coords)
  beta.size <- ncol(xmat)
  z <- data
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  ##
  ## Absurd values
  ##
  if(kappa < 1e-04) return(1e+32)
  if(round(1e+16*(nugget+sigmasq)) == 0) return(1e+32)
  ##
  ## Anisotropy
  ##
  if(psiR != 1 | psiA != 0){
    coords.c <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    .likGRF.dists.vec <-  as.vector(dist(coords.c))
  }
  else
    if(compute.dists) .likGRF.dists.vec <-  as.vector(dist(coords))
  ##
  ## Box-Cox transformation
  ##
  if(round(lambda, dig=4) == 1) {
    log.jacobian <- 0
  }
  else {
    if(any(z <= 0))
      stop("Transformation not allowed for zero or negative data")
    data <- z^(lambda - 1)
    if(any(data <= 0)) log.jacobian <- log(prod(data))
    else log.jacobian <- sum(log(data))
    data <- NULL
    if(round(lambda, dig=4) == 0)
      data <- log(z)
    else data <- ((z^lambda) - 1)/lambda
  }
  ##
  ## Computing likelihood
  ##
  ## NOTE: Likelihood for Independent observations 
  ##       arbitrary criteria used here:
  ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
  ##
  if((phi < 1e-16) | (sigmasq < 1e-16)){
    iv <- list(sqrt.inverse = diag(x=1/sqrt((nugget+sigmasq)), n),
               log.det.to.half = (n/2) * log(nugget+sigmasq))
  }
  else{
    iv <- varcov.spatial(dists.lowertri = .likGRF.dists.vec,
                         cov.model = cov.model, kappa=kappa,
                         nugget = nugget, cov.pars=c(sigmasq, phi),
                         sqrt.inv = TRUE, det = TRUE)
  }
  sivx <- crossprod(iv$sqrt.inverse, xmat)
  xivx <- crossprod(sivx)
  sivy <- crossprod(iv$sqrt.inverse, data)
  xivy <- crossprod(sivx, sivy)  
  betahat <- solve(xivx, xivy)
  res <- data - xmat %*% betahat
  ssres <- as.vector(crossprod(crossprod(iv$sqrt.inverse,res)))
  if(method == "ML"){
    negloglik <- (n/2)*(log(2*pi)) + iv$log.det.to.half +  0.5 * ssres - log.jacobian
  }
  if(method == "RML"){
    if(length(as.vector(xivx)) == 1) {
      choldet <- 0.5 * log(xivx)
    }
    else {
      chol.xivx <- chol(xivx)
      choldet <- sum(log(diag(chol.xivx)))
    }
    negloglik <- iv$log.det.to.half +  0.5 * ssres + choldet - log.jacobian
    xx.eigen <- eigen(crossprod(xmat), symmetric = TRUE, only.values = TRUE)
    negloglik <- negloglik + ((n-beta.size)/2)*(log(2*pi)) - 0.5 * sum(log(xx.eigen$values))
  }
  if(negloglik > 1e+32) negloglik <- 1e32
  return(-negloglik)
}



