"likfit" <-
  function (geodata, coords=geodata$coords, data=geodata$data,
            trend = "cte", ini.cov.pars,
            fix.nugget = FALSE, nugget = 0, 
            fix.kappa = TRUE, kappa = 0.5, 
            fix.lambda = TRUE, lambda = 1, 
            fix.psiA = TRUE, psiA = 0, 
            fix.psiR = TRUE, psiR = 1, 
            cov.model = "matern", realisations,
            method.lik = "ML",
            components = FALSE, nospatial = TRUE,
            limits = likfit.limits(), 
            print.pars = FALSE, messages.screen = TRUE, ...) 
{
  if(is.R()) require(mva)
  ##
  ## Checking input
  ##
  call.fc <- match.call()
  temp.list <- list()
  temp.list$print.pars <- print.pars
  ##
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(fix.kappa & !is.null(kappa))
    if(cov.model == "matern" & kappa == 0.5)
      cov.model <- "exponential"
  temp.list$cov.model <- cov.model
  if(cov.model == "powered.exponential")
    if(limits$kappa["upper"] > 2) limits$kappa["upper"] <- 2
  ##
  ## Likelihood method
  ##
  if(method.lik == "REML" | method.lik == "reml" | method.lik == "rml") 
    method.lik <- "RML"
  if(method.lik == "ML" | method.lik == "ml")
    method.lik <- "ML"
  if(method.lik == "ML" & cov.model == "power")
    stop("\n\"power\" model can only be used with method.lik=\"RML\".\nBe sure that what you want is not \"powered.exponential\"")
  temp.list$method.lik <- method.lik
  ##
  ## Coodinates, data and covariates
  ##
  coords <- as.matrix(coords)
  data <- as.vector(data)
  n <- length(data)
  if((nrow(coords) != n) | (2*n) != length(coords))
    stop("\nnumber of locations does not match with number of data")
  if(missing(geodata))
    xmat <- trend.spatial(trend=trend, geodata=list(coords = coords, data = data))
  else
    xmat <- trend.spatial(trend=trend, geodata=geodata)
  if(nrow(xmat) != n)
    stop("trend matrix has dimension incompatible with the data")
  test.xmat <- solve.geoR(crossprod(xmat))
  test.xmat <- NULL
  beta.size <- temp.list$beta.size <- dim(xmat)[2]
  ##
  ## setting a factor for indicating different realisations
  ##
  if(missing(realisations))
    realisations <- as.factor(rep(1, n))
  else{
    if(!missing(geodata)){
      real.name <- deparse(substitute(realisations))
      if(!is.null(geodata[[real.name]]))
        realisations <- geodata$realisations
    }
    if(length(realisations) != n)
      stop("realisations must be a vector with the same length of the data")
    realisations <- as.factor(realisations)
  }
  temp.list$realisations <- realisations
  nrep <- temp.list$nrep <- length(levels(realisations))
  ind.rep <- split(1:n, realisations)
  vecdist <- function(x){as.vector(dist(x))}
  ##
  ## Initial values for parameters
  ##
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
  ##
  ## Checking for multiple initial values for preliminar search of   
  ## best initial value
  ##
  if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1) | (length(lambda) > 1) | (length(psiR) > 1) | (length(psiA) > 1)){
    if(messages.screen)
      cat("likfit: searching for best initial value ...")
    ini.temp <- matrix(ini.cov.pars, ncol=2)
    grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]), phi=unique(ini.temp[,2]), tausq=unique(nugget), kappa=unique(kappa), lambda=unique(lambda), psiR=unique(psiR), psiA=unique(psiA)))
    .likGRF.dists.vec <<- lapply(split(as.data.frame(coords), realisations), vecdist)
    temp.f <- function(parms, coords, data, temp.list)
      return(loglik.GRF(coords = coords, data = as.vector(data),
                        cov.model=temp.list$cov.model, cov.pars=parms[1:2],
                        nugget=parms["tausq"], kappa=parms["kappa"],
                        lambda=parms["lambda"], psiR=parms["psiR"],
                        psiA=parms["psiA"], trend= trend,
                        method.lik=temp.list$method.lik,
                        compute.dists=FALSE,
                        realisations = realisations))
    grid.lik <- apply(grid.ini, 1, temp.f, coords = coords,
                      data = data, temp.list = temp.list)
    grid.lik <- grid.lik[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik)] 
    grid.ini <- grid.ini[(grid.lik != Inf) & (grid.lik != -Inf) & !is.na(grid.lik) & !is.nan(grid.lik),, drop=FALSE] 
    ini.temp <- grid.ini[which(grid.lik == max(grid.lik)),, drop=FALSE]
    if(is.R()) rownames(ini.temp) <- "initial.value"
    if(messages.screen){
      cat(" selected values:\n")
      print(rbind(round(ini.temp, dig=2), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa, fix.lambda, fix.psiR, fix.psiA), "fix", "est")))
      cat(paste("likelihood value:", max(grid.lik), "\n"))
    }
    dimnames(ini.temp) <- NULL
    ini.cov.pars <- ini.temp[1:2]
    nugget <- ini.temp[3]
    kappa <- ini.temp[4]
    lambda <- ini.temp[5]
    psiR <- ini.temp[6]
    psiA <- ini.temp[7]
    grid.ini <- NULL
    if(is.R()) {remove(".likGRF.dists.vec", pos=1); gc(verbose=FALSE)}    
    else remove(".likGRF.dists.vec", where=1)    
  }
  ##
  tausq <- nugget
  ##
  ## Box-Cox transformation for fixed lambda
  ##
  if(fix.lambda) {
    if(abs(lambda - 1) < 0.0001) {
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
      if(abs(lambda) < 0.0001)
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
    if(is.R())
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords), realisations), vecdist), pos=1)
    else
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords), realisations), vecdist), where=1)
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
  if(messages.screen == TRUE) {
    cat("-------------------------------------------------------\n")
    cat("likfit: Initialising likelihood maximisation using the function ")
    if(is.R()) cat("optim.\n") else cat("nlminb.\n")
    cat("likfit: Use control() to pass arguments for the maximisation function.")
    cat("\n        For more details see documentation for ")
    if(is.R()) cat("optim.\n") else cat("nlminb.\n")        
    cat("likfit: It is highly advisable to run this function several\n        times with different initial values for the parameters.\n")
    cat("likfit: WARNING: This step can be time demanding!\n")
    cat("-------------------------------------------------------\n")
  }
  ##
  npars <- beta.size + 2 + sum(unlist(ip)==FALSE)
  temp.list$coords <- coords
  temp.list$xmat <- split(as.data.frame(xmat), realisations)
  temp.list$xmat <- lapply(temp.list$xmat, as.matrix)
  temp.list$n <- as.vector(unlist(lapply(temp.list$xmat, nrow)))
  ##
  ## Constant term in the likelihood
  ##
  temp.list$loglik.cte <- rep(0, nrep)
  for(i in 1:nrep){
    if(method.lik == "ML"){
      if(ip$f.tausq & (tausq > 0))
        temp.list$loglik.cte[i] <-  (temp.list$n[i]/2)*(-log(2*pi))
      else
        temp.list$loglik.cte[i] <-  (temp.list$n[i]/2)*(-log(2*pi) +
                                                        log(temp.list$n[i]) -1)
    }
    if(method.lik == "RML"){
      xx.eigen <- eigen(crossprod(temp.list$xmat[[i]]),
                        symmetric = TRUE, only.values = TRUE)
      if(ip$f.tausq & (tausq > 0))
        temp.list$loglik.cte[i] <- - ((temp.list$n[i]-beta.size)/2)*(log(2*pi)) +
          0.5 * sum(log(xx.eigen$values))
      else
        temp.list$loglik.cte[i] <-  - ((temp.list$n[i]-beta.size)/2)*(log(2*pi)) +
          ((temp.list$n[i]-beta.size)/2)*(log(temp.list$n[i]-beta.size)) -
            ((temp.list$n[i]-beta.size)/2) + 0.5 * sum(log(xx.eigen$values))
    }
  }
  ##
  ## Numerical minimization of the -loglikelihood
  ##
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
    loglik.max <- - lik.minim$value
  else
    loglik.max <- - lik.minim$objective
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
    if(abs(lambda - 1) < 0.0001) {
      log.jacobian.max <- 0
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
  data.rep <- split(temp.list$z, realisations)
  coords.rep <- split(as.data.frame(coords), realisations)
  coords.rep <- lapply(coords.rep, as.matrix)
  ##
  ## Transforming coords for estimated anisotropy (if the case)
  ##
  if(fix.psiR & fix.psiA){
    if(is.R()) remove(".likGRF.dists.vec", pos=1)
    else remove(".likGRF.dists.vec", where=1)
  }
  else{
    if(round(psiR, dig=6) != 1 | round(psiA, dig=6) != 0)
      coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    rangevecdist <- function(x){range(as.vector(dist(x)))}
    range.dist <- lapply(split(as.data.frame(coords), realisations), rangevecdist)
    range.dist <- range(as.vector(unlist(range.dist)))
    max.dist <- max(range.dist)
    min.dist <- min(range.dist)
  }      
  if(is.R()) gc(verbose=FALSE)
  ##
  ## Computing estimated beta and tausq/sigmasq (if the case)
  ##
  xivx <- matrix(0, ncol=beta.size, nrow=beta.size)
  xivy <- matrix(0, ncol=1, nrow=beta.size)
  yivy <- 0
  for(i in 1:nrep){
    ni <- temp.list$n[i]
    xmati <- temp.list$xmat[[i]]
    if((phi < 1e-12))
      siv <- diag(x=1/sqrt((1+tausq)), ni)
    else{
      if(check.sigmasq){
        if(sigmasq < 1e-12){
          if(!fix.nugget)
            siv <- diag(x=1/sqrt((1+tausq)), ni)
          else
            siv <- diag(x=1/sqrt((tausq)), ni)          
        }
        else
          siv <- varcov.spatial(coords = coords.rep[[i]], cov.model = cov.model,
                                kappa = kappa,
                                nugget = tausq, cov.pars = c(1, phi),
                                inv=TRUE, sqrt.inv = TRUE,
                                det = FALSE)$sqrt.inverse
      }
      else
        siv <- varcov.spatial(coords = coords.rep[[i]], cov.model = cov.model,
                              kappa = kappa,
                              nugget = tausq, cov.pars = c(1, phi),
                              inv=TRUE, sqrt.inv = TRUE,
                              det = FALSE)$sqrt.inverse
    }
    sivx <- crossprod(siv, temp.list$xmat[[i]])
    xivx <- xivx + crossprod(sivx)
    sivy <- crossprod(siv, data.rep[[i]])
    xivy <- xivy + crossprod(sivx, sivy)
    yivy <- yivy + crossprod(sivy)
  }
  betahat <- solve.geoR(xivx, xivy)
  res <- as.vector(temp.list$z - xmat %*% betahat)
  if(!fix.nugget | (nugget < 1e-12)){
    ssres <- as.vector(yivy - 2*crossprod(betahat,xivy) +
                       crossprod(betahat,xivx) %*% betahat)  
    if(method.lik == "ML")
      sigmasq <- ssres/n
    else
      sigmasq <- ssres/(n - beta.size)
  }
  if(fix.nugget){
    if(nugget > 0)
      tausq <- nugget
  }
  else tausq <- tausq * sigmasq
  betahat.var <- solve.geoR(xivx)
  if(sigmasq > 1e-12) betahat.var <- sigmasq * betahat.var
#  if(!fix.nugget & phi < 1e-16){
#    tausq <- sigmasq + tausq
#    sigmasq <- 0
#  }
  n.model.pars <- beta.size + 7
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
                      sigmasq = sigmasq,
                      phi = phi,
                      kappa = kappa,
                      beta = as.vector(betahat),
                      beta.var = betahat.var,
                      lambda = lambda,
                      aniso.pars = c(psiA = psiA, psiR = psiR),
                      method.lik = method.lik, trend = trend,
                      loglik = loglik.max,
                      npars = npars,
                      AIC = -2 * (loglik.max - npars),
                      BIC = -2 * (loglik.max - 0.5 * log(n) * npars),
                      parameters.summary = par.su,
                      info.minimisation.function = lik.minim,
                      max.dist = max.dist,
                      trend.matrix= xmat,
                      transform.info = list(fix.lambda = fix.lambda,
                        log.jacobian = log.jacobian.max))
  ##
  ## Likelihood results for the model without spatial correlation
  ##
  if(nospatial){
    if(fix.lambda){
      beta.ns <- solve.geoR(crossprod(xmat), crossprod(xmat, temp.list$z))
      ss.ns <- sum((as.vector(temp.list$z - xmat %*% beta.ns))^2)
      if(method.lik == "ML"){
        nugget.ns <- ss.ns/n
        loglik.ns <- (n/2)*((-log(2*pi)) - log(nugget.ns) - 1) + temp.list$log.jacobian
      }
      if(method.lik == "RML"){
        nugget.ns <- ss.ns/(n-beta.size)
        loglik.ns <- ((n-beta.size)/2)*((-log(2*pi)) - log(nugget.ns) -1) +
          temp.list$log.jacobian
      }
      npars.ns <- beta.size + 1 + fix.lambda
      lambda.ns <- lambda
    }
    else{
      if(is.R())
        lik.lambda.ns <- optim(par=1, fn = boxcox.negloglik, method = "L-BFGS-B",
                               lower = limits$lambda["lower"],
                               upper = limits$lambda["upper"],
                               data = data, xmat = xmat, lik.method = method.lik)
      else
        lik.lambda.ns <- nlminb(par=1, fn = boxcox.negloglik,
                                lower=limits$lambda["lower"],
                                upper=limits$lambda["upper"],
                                data = data, xmat = xmat, lik.method = method.lik)
      lambda.ns <- lik.lambda.ns$par
      if(abs(lambda) < 0.0001) tdata.ns <- log(data)
      else tdata.ns <- ((data^lambda.ns)-1)/lambda.ns
      beta.ns <- solve.geoR(crossprod(xmat),crossprod(xmat,tdata.ns))
      ss.ns <- sum((as.vector(tdata.ns - xmat %*% beta.ns))^2)
      if(is.R())
        value.min.ns <- lik.lambda.ns$value
      else
        value.min.ns <- lik.lambda.ns$objective
      if(method.lik == "ML"){
        loglik.ns <- (- value.min.ns)+ (n/2)*((-log(2*pi)) + log(n) - 1)
        nugget.ns <- ss.ns/n
      }
      if(method.lik == "RML"){
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
  ## Assigning names to the components of the mean vector beta
  ##
  if(length(lik.results$beta.var) == 1)
    lik.results$beta.var <- as.vector(lik.results$beta.var)
  if(length(lik.results$beta) > 1){
    if(inherits(trend, "formula"))
      beta.names <- c("intercept", paste("covar", 1:(ncol(xmat)-1), sep = ""))
    else
      if(trend == "1st")
        beta.names <- c("intercept", "x", "y")
      else
        if(trend == "2nd")
          beta.names <- c("intercept", "x", "y", "x2", "xy", "y2")
    names(lik.results$beta) <- beta.names
  }
  ##
  ## Computing residuals and predicted values
  ## (isolated components of the model)
  ##
  if(components) {
    if(!fix.psiR & !fix.psiA)
      if(psiR != 1 | psiA != 0)
        coords <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    coords.rep <- split(as.data.frame(coords), realisations)
    res.rep <- split(res, realisations)
    trend.comp <- temp.list$z - res
    spatial.comp <- list()
    for(i in 1:nrep){
      invcov <- varcov.spatial(coords = coords[ind.rep[[i]],], cov.model = cov.model, 
                               kappa = kappa, nugget = tausq,
                               cov.pars = c(sigmasq, phi), inv=TRUE)$inverse 
      covmat.signal <- varcov.spatial(coords = coords[ind.rep[[i]],],
                                      cov.model = cov.model, 
                                      kappa = kappa, nugget = 0,
                                      cov.pars = c(sigmasq, phi))$varcov
      spatial.comp[[i]] <- as.vector(covmat.signal %*% invcov %*% res[ind.rep[[i]]])
    }
    spatial.comp <- as.vector(unlist(spatial.comp))[as.vector(unlist(ind.rep))]
    predict.comp <- trend.comp + spatial.comp
    residual.comp <- as.vector(temp.list$z - predict.comp)
#    residual.std <- as.vector(invcov %*% residual.comp)
#    residual.trend.std <- as.vector(invcov %*% res)
    lik.results$model.components <-
      data.frame(trend = trend.comp, spatial = spatial.comp, residuals = residual.comp)
#    lik.results$s2.random <- (crossprod(res,invcov) %*% res)/(n - beta.size)
#    lik.results$s2 <- (crossprod(residual.comp,invcov) %*% residual.comp)/(n - beta.size)
  }
  ##
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
  if(temp.list$print.pars){
    running.pars <- c(sigmasq, phi, tausq, kappa, psiA, psiR, lambda)
    names(running.pars) <-
      c("sigmasq", "phi", "tausq", "kappa", "psiA", "psiR", "lambda")
    print(running.pars)
  }
  ##
  ## Absurd values
  ##
  if(kappa < 1e-04) return(.Machine$double.xmax/10000)
  if((tausq+sigmasq) < 1e-16) return(.Machine$double.xmax/10000)
  ##
  ## Anisotropy
  ##
  if(!ip$f.psiR | !ip$f.psiA){
    coords.c <- coords.aniso(temp.list$coords, aniso.pars=c(psiA, psiR))
    vecdist <- function(x){as.vector(dist(x))}
    if(is.R())
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords.c),
                                               temp.list$realisations), vecdist), pos=1)
    else
      assign(".likGRF.dists.vec", lapply(split(as.data.frame(coords.c),
                                               temp.list$realisations), vecdist), where=1)
  }
  ##
  ## Box-Cox transformation
  ##
  if(!ip$f.lambda){
    if(abs(lambda - 1) < 0.0001) {
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
    if(abs(lambda) < 0.0001)
      data <- log(temp.list$z)
    else data <- ((temp.list$z^lambda) - 1)/lambda
  }
  else data <- temp.list$z
  data <- split(data, as.factor(temp.list$realisations))
  ##
  ## Computing likelihood
  ##
  sumnegloglik <- 0
  for(i in 1:temp.list$nrep){
    ## NOTE: Likelihood for Independent observations 
    ##       arbitrary criteria used here:
    ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
    ##
    n <- temp.list$n[i]
    xmat <- temp.list$xmat[[i]]
    z <- data[[i]]
    if((phi < 1e-16) | (sigmasq < 1e-16)){
      if(ip$f.tausq)
        iv <- list(sqrt.inverse = diag(x=1/sqrt((tausq+sigmasq)), n),
                   log.det.to.half = (n/2) * log(tausq+sigmasq))
      else
        iv <- list(sqrt.inverse = diag(x=1/sqrt((1+tausq)), n),
                   log.det.to.half = (n/2) * log(1+tausq))
    }
    else{
      iv <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]],
                           cov.model = temp.list$cov.model, kappa=kappa,
                           nugget = tausq, cov.pars=c(sigmasq, phi),
                           sqrt.inv = TRUE, det = TRUE)
    }
    if(!is.null(iv$crash.parms)) return(.Machine$double.xmax/100)
    sivx <- crossprod(iv$sqrt.inverse, xmat)
    xivx <- crossprod(sivx)
    sivy <- crossprod(iv$sqrt.inverse, z)
    xivy <- crossprod(sivx, sivy)
    betahat <- solve.geoR(xivx, xivy)
    if(inherits(betahat, "try-error")){
      error.now <- options()$show.error.message
      options(show.error.messages = FALSE)
      t.ei <- eigen(xivx, symmetric = TRUE)
      betahat <- try(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec) %*% xivy)
      if(is.null(error.now) || error.now == TRUE)
        options(show.error.messages = TRUE)        
    }
    if(inherits(betahat, "try-error"))
      stop("Covariates have very different orders of magnitude. Try to multiply and/or divide them to bring them to similar orders of magnitude") 
    res <- z - xmat %*% betahat
    ssres <- as.vector(crossprod(crossprod(iv$sqrt.inverse,res)))
    if(temp.list$method.lik == "ML"){
      if(ip$f.tausq & (tausq > 0))
        negloglik <- iv$log.det.to.half +  0.5 * ssres
      else
        negloglik <- (n/2) * log(ssres) +  iv$log.det.to.half
    }
    if(temp.list$method.lik == "RML"){
      if(length(as.vector(xivx)) == 1) {
        choldet <- 0.5 * log(xivx)
      }
      else {
        chol.xivx <- chol(xivx)
        choldet <- sum(log(diag(chol.xivx)))
      }
      if(ip$f.tausq & (tausq > 0))
        negloglik <- iv$log.det.to.half +  0.5 * ssres + choldet
      else
        negloglik <- ((n-p)/2) * log(ssres) +  iv$log.det.to.half + choldet
    }  
    negloglik <- negloglik - temp.list$loglik.cte[i]
    sumnegloglik <- sumnegloglik + negloglik
  }
  sumnegloglik <- sumnegloglik - log.jacobian
  if(sumnegloglik > .Machine$double.xmax/1000 | sumnegloglik == Inf | sumnegloglik == -Inf)
    sumnegloglik <- (.Machine$double.xmax/1000)
  if(temp.list$print.pars)
    cat(paste("log-likelihood = ", -sumnegloglik, "\n"))
  return(sumnegloglik) 
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
  function(x, digits = "default", ...)
{
  if(is.R() & digits == "default") digits <- max(3, getOption("digits") - 3)
  else digits <- options()$digits
  est.pars <- as.vector(x$parameters.summary[x$parameters.summary[,1] == "estimated",2])
  names.est.pars <- dimnames(x$parameters.summary[x$parameters.summary[,1] == "estimated",])[[1]]
  names(est.pars) <- names.est.pars
  cat("likfit: estimated model parameters:\n")
  print(round(est.pars, digits=digits))
  cat("\nlikfit: maximised log-likelihood = ")
  cat(round(x$loglik, digits=digits))
  cat("\n")
  return(invisible())
}  

"summary.likGRF" <-
  function(object, ...)
{
  names.pars <- dimnames(object$parameters.summary)[[1]]
  summ.lik <- list()
  if(object$method.lik == "ML")
    summ.lik$method.lik <- "maximum likelihood"
  if(object$method.lik == "RML")
    summ.lik$method.lik <- "restricted maximum likelihood"
  summ.lik$mean.component <- object$beta
  names(summ.lik$mean.component) <- names.pars[1:length(object$beta)]
  summ.lik$cov.model <- object$cov.model
  summ.lik$spatial.component <- object$parameters.summary[c("sigmasq", "phi"),]
  summ.lik$spatial.component.extra <- object$parameters.summary[c("kappa", "psiA", "psiR"),]
  summ.lik$nugget.component <- object$parameters.summary[c("tausq"),, drop=FALSE]
  summ.lik$transformation  <- object$parameters.summary[c("lambda"),, drop=FALSE]
  summ.lik$likelihood <- c(log.L = object$loglik, n.params = as.integer(object$npars),
                               AIC = object$AIC, BIC = object$BIC)
  summ.lik$estimated.pars <- dimnames(object$parameters.summary[object$parameters.summary[,1] == "estimated",])[[1]]
  likelihood.info <- c(log.L = object$loglik, n.params = as.integer(object$npars),
                       AIC = object$AIC, BIC = object$BIC)
  summ.lik$call <- object$call
  class(summ.lik) <- "summary.likGRF"
  return(summ.lik)
}

"print.summary.likGRF" <-
  function(x, digits = "default", ...)
{
  if(class(x) != "summary.likGRF")
    stop("object is not of the class \"summary.likGRF\"")
  if(is.R() & digits == "default") digits <- max(3, getOption("digits") - 3)
  else digits <- options()$digits
  cat("Summary of the parameter estimation\n")
  cat("-----------------------------------\n")
  cat(paste("Estimation method:", x$method.lik, "\n"))
  cat("\n")
  ##
  ## Estimates of the model components
  ## Model: Y(x) = X\beta + S(x) + e 
  ##
  cat("Parameters of the mean component (trend):")
  cat("\n")
  print(round(x$mean.component, digits=digits))
  cat("\n")
  ##
  cat("Parameters of the spatial component:")
  cat("\n")
  cat(paste("   correlation function:", x$cov.model))
  cat(paste("\n      (estimated) variance parameter sigmasq (partial sill) = ", round(x$spatial.component[1,2], dig=digits)))
  cat(paste("\n      (estimated) cor. fct. parameter phi (range parameter)  = ", round(x$spatial.component[2,2], dig=digits)))
  if(x$cov.model == "matern" | x$cov.model == "powered.exponential" |
     x$cov.model == "cauchy" | x$cov.model == "gneiting.matern"){
    kappa <- x$spatial.component.extra["kappa",2]
    if(x$spatial.component.extra["kappa",1] == "estimated")
      cat(paste("\n      (estimated) extra parameter kappa =", round(kappa, digits=digits)))
    else{
      cat(paste("\n      (fixed) extra parameter kappa = ", kappa))
      if(x$cov.model == "matern" & (round(kappa, digits=digits)  == 0.5))
      cat(" (exponential)")
    }
  }
  cat("\n")
  ##
  aniso <-  x$spatial.component.extra[c("psiA", "psiR"),]
  psiApsiR <- x$spatial.component.extra[c("psiA", "psiR"),2]
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
  if(x$nugget.component[,1] == "estimated")
    cat(paste("\n      (estimated) nugget = ", round(x$nugget.component[,2], dig=digits)))
  else
    cat(paste("\n      (fixed) nugget =", x$nugget.component[,2]))
  cat("\n")
  cat("\n")
  cat("Transformation parameter:")
  cat("\n")
  lambda <- x$transformation[,2]
  if(x$transformation[,1] == "estimated")
    cat(paste("      (estimated) Box-Cox parameter =", round(lambda, dig=digits)))
  else{
    cat(paste("      (fixed) Box-Cox parameter =", lambda))
    if(abs(lambda - 1) <  0.0001) cat(" (no transformation)")
    if(abs(lambda) < 0.0001) cat(" (log-transformation)")
  }
  cat("\n")
  cat("\n")
  cat("Maximised Likelihood:")
  cat("\n")
  print(round(x$likelihood, digits=digits))
  cat("\n")
  cat("Call:")
  cat("\n")
  print(x$call)
  cat("\n")
  invisible(x)
}

"loglik.GRF" <-
  function(geodata, coords=geodata$coords, data=geodata$data,
           cov.model="exp", cov.pars,
           nugget=0, kappa=0.5, lambda=1, psiR=1, psiA=0,
           trend="cte", method.lik="ML",
           compute.dists = TRUE, realisations = NULL)
{
  if(method.lik == "REML" | method.lik == "reml" | method.lik == "rml") 
    method.lik <- "RML"
  if(method.lik == "ML" | method.lik == "ml")
    method.lik <- "ML"
  sigmasq <- cov.pars[1]
  phi <- cov.pars[2]
  if(is.null(realisations))
    realisations <- as.factor(rep(1, length(data)))
  else
    realisations <- as.factor(realisations)
  nrep <- length(levels(realisations))
  ##
  ## Absurd values
  ##
  if(kappa < 1e-04) return(-(.Machine$double.xmax/100))
  if((nugget+sigmasq) < 1e-16) return(-(.Machine$double.xmax/100))
  ##
  ## Trend matrix
  ##
  if(missing(geodata))
    xmat <- trend.spatial(trend=trend, geodata = list(coords = coords, data = data))
  else
    xmat <- trend.spatial(trend=trend, geodata = geodata)
  beta.size <- ncol(xmat)
  xmat <- split(as.data.frame(xmat), realisations)
  xmat <- lapply(xmat, as.matrix)
  ##
  ## Anisotropy
  ##
  vecdist <- function(x){as.vector(dist(x))}
  if(psiR != 1 | psiA != 0){
    coords.c <- coords.aniso(coords, aniso.pars=c(psiA, psiR))
    vecdist <- function(x){as.vector(dist(x))}
    .likGRF.dists.vec <- lapply(split(as.data.frame(coords.c),
                                      as.factor(realisations)), vecdist)
  }
  else
    if(compute.dists)
      .likGRF.dists.vec <- lapply(split(as.data.frame(coords),
                                        as.factor(realisations)), vecdist)
  ##
  ## Box-Cox transformation
  ##
  z <- data
  if(abs(lambda - 1) < 0.0001)
    log.jacobian <- 0
  else {
    if(any(z <= 0))
      stop("Transformation not allowed for zero or negative data")
    data <- z^(lambda - 1)
    if(any(data <= 0)) log.jacobian <- log(prod(data))
    else log.jacobian <- sum(log(data))
    data <- NULL
    if(abs(lambda) < 0.0001)
      data <- log(z)
    else data <- ((z^lambda) - 1)/lambda
  }
  data <- split(data, as.factor(realisations))
  ##
  ## Computing likelihood
  ##
  sumnegloglik <- 0
  for(i in 1:nrep){
    ## NOTE: Likelihood for Independent observations 
    ##       arbitrary criteria used here:
    ##       (phi < 1-e16) or (sigmasq < 1-e16)  ==> independence
    ##
    n <- length(data[[1]])
    if((phi < 1e-16) | (sigmasq < 1e-16)){
      iv <- list(sqrt.inverse = diag(x=1/sqrt((nugget+sigmasq)), n),
                 log.det.to.half = (n/2) * log(nugget+sigmasq))
    }
    else{
      iv <- varcov.spatial(dists.lowertri = .likGRF.dists.vec[[i]],
                           cov.model = cov.model, kappa=kappa,
                           nugget = nugget, cov.pars=c(sigmasq, phi),
                           sqrt.inv = TRUE, det = TRUE)
    }
    if(!is.null(iv$crash.parms)){
      cat("varcov.spatial: improper matrix for following the given parameters:")
      print(iv$crash.parms)
      stop()
    }
    sivx <- crossprod(iv$sqrt.inverse, xmat[[i]])
    xivx <- crossprod(sivx)
    sivy <- crossprod(iv$sqrt.inverse, data[[i]])
    xivy <- crossprod(sivx, sivy)  
    betahat <- solve.geoR(xivx, xivy)
    res <- data[[i]] - xmat[[i]] %*% betahat
    ssres <- as.vector(crossprod(crossprod(iv$sqrt.inverse,res)))
    if(method.lik == "ML"){
      negloglik <- (n/2)*(log(2*pi)) + iv$log.det.to.half +  0.5 * ssres
    }
    if(method.lik == "RML"){
      if(length(as.vector(xivx)) == 1) {
        choldet <- 0.5 * log(xivx)
      }
      else {
        chol.xivx <- chol(xivx)
        choldet <- sum(log(diag(chol.xivx)))
      }
      negloglik <- iv$log.det.to.half +  0.5 * ssres + choldet
      xx.eigen <- eigen(crossprod(xmat[[i]]), symmetric = TRUE, only.values = TRUE)
      negloglik <- negloglik + ((n-beta.size)/2)*(log(2*pi)) - 0.5 * sum(log(xx.eigen$values))
    }
    sumnegloglik <- sumnegloglik + negloglik 
  }
  sumnegloglik <- sumnegloglik - log.jacobian
  if(sumnegloglik > .Machine$double.xmax/1000)
    sumnegloglik <- (.Machine$double.xmax/1000)
  return(as.vector(-sumnegloglik))
}

"solve.geoR" <-
  function(a, b = NULL, ...)
{
  error.now <- options()$show.error.message
  options(show.error.messages = FALSE)
  if(is.null(b))
    res <- try(solve(a, ...))
  else
    res <- try(solve(a, b, ...))
  if(inherits(res, "try-error")){
    t.ei <- eigen(a, symmetric = TRUE)
    if(is.null(b))
      res <- try(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec))
    else
      res <- try(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec) %*% b)
  }  
  if(is.null(error.now) || error.now == TRUE)
    options(show.error.messages = TRUE)        
  if(inherits(res, "try-error"))
    stop("Covariates have very different orders of magnitude. Try to multiply and/or divide them to bring them to similar orders of magnitude")
  return(res)
}

