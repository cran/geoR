"olsfit" <-
  function(...)
{
  stop("this function is now obsolete.\nuse variofit() instead.")
}

"wlsfit" <-
  function(...)
{
  stop("this function is now obsolete.\nuse variofit() instead.")
}

"variofit" <-
  function (vario, ini.cov.pars, cov.model = "matern",
            fix.nugget = FALSE, nugget = 0, 
            fix.kappa = TRUE, kappa = 0.5,
            simul.number = NULL,  max.dist = "all",
            weights = c("npairs", "equal", "cressie"),
            minimisation.function,
            messages.screen = TRUE, ...) 
{
  call.fc <- match.call()
  if(class(vario) != "variogram")
    warning("object vario should preferably  be of the class \"variogram\"")
  weights <- match.arg(weights)
  if(missing(minimisation.function)){
    if(weights == "equal") minimisation.function <- "nls"
    else minimisation.function <- "optim"
  }
  if(minimisation.function == "nls" & weights != "equal"){
    warning("variofit: minimisation function nls can only be used with weights=\"equal\".\n          changed for \"optim\".\n")
    minimisation.function <- "optim"
  }
  if(messages.screen){
    cat(paste("variofit: weights used:", weights, "\n"))
    cat(paste("variofit: minimisation function used:", minimisation.function, "\n"))
  }
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (is.matrix(vario$v) & is.null(simul.number)) 
    stop("object in vario$v is a matrix. This function works for only 1 empirical variogram at once\n")
  if (!is.null(simul.number)) 
    vario$v <- vario$v[, simul.number]
  ##
  ## Setting maximum distance
  ##
  if (max.dist == "all") 
    XY <- data.frame(u = vario$u, v = vario$v, n=vario$n)
  else
    XY <- data.frame(u = vario$u[vario$u <= max.dist],
                     v = vario$v[vario$u <= max.dist],
                     n = vario$v[vario$u <= max.dist])
  ##
  ##  Checking initial values
  ##
  if(is.matrix(ini.cov.pars) | is.data.frame(ini.cov.pars)){
    ini.cov.pars <- as.matrix(ini.cov.pars)
    if(nrow(ini.cov.pars) == 1)
      ini.cov.pars <- as.vector(ini.cov.pars)
    else{
      if((cov.model != "pure.nugget") & (ncol(ini.cov.pars) != 2))
        stop("\nini.cov.pars must be a matrix or data.frame with 2 components: \ninitial values for sigmasq (partial sill) and phi (range parameter)\n")
    }
  }
  if(is.vector(ini.cov.pars)){
    if((cov.model != "pure.nugget") & (length(ini.cov.pars) != 2))
      stop("\nini.cov.pars must be a vector with 2 components: \ninitial values for sigmasq and phi\n")
  }
  ##
  ## Preparing grid of initial values and choosing the best
  ##
  if(is.matrix(ini.cov.pars) | (length(nugget) > 1) | (length(kappa) > 1)) {
    if(messages.screen)
      cat("variofit: searching for best initial value ...")
    ini.temp <- matrix(ini.cov.pars, ncol=2)
    grid.ini <- as.matrix(expand.grid(sigmasq=unique(ini.temp[,1]), phi=unique(ini.temp[,2]), tausq=unique(nugget), kappa=unique(kappa)))
    ##  loss function:
    v.loss <- function(parms, u, v, n, cov.model, weights){
      sigmasq <- parms[1]
      phi <- parms[2]
      tausq <- parms[3]
      kappa <- parms[4]
      v.mod <- (sigmasq + tausq) -
        cov.spatial(u, cov.pars=c(sigmasq, phi), cov.model = cov.model,
                    kappa = kappa)
      if(weights == "equal")
        loss <- sum((v - v.mod)^2)
      if (weights == "npairs") 
        loss <- sum(n * (v - v.mod)^2)
      if (weights == "cressie") 
        loss <- sum((n/(v.mod^2)) * (v - v.mod)^2)
      return(loss)
    }
    grid.loss <- apply(grid.ini, 1, v.loss, u=XY$u, v=XY$v, n=XY$n, cov.model = cov.model, weights = weights)
    ini.temp <- grid.ini[which(grid.loss == min(grid.loss)),, drop=FALSE]
    if(is.R()) rownames(ini.temp) <- "initial.value"
    if(messages.screen){
      cat(" selected values:\n")
      print(rbind(round(ini.temp, dig=2), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa), "fix", "est")))
      cat(paste("loss value:", max(grid.loss), "\n"))
    }
    names(ini.temp) <- NULL
    ini.cov.pars <- ini.temp[1:2]
    nugget <- ini.temp[3]
    kappa <- ini.temp[4]
    grid.ini <- NULL
  }
  ##
  ## transforming kappa for constraint minimisation
  ##
  if(fix.kappa == FALSE){
    if(cov.model == "powered.exponential")
      Tkappa.ini <- log(kappa/(2-kappa))
    else
      Tkappa.ini <- log(kappa)
  }
  ##
  ## minimisation using "nls"
  ##
  if (minimisation.function == "nls") {
    require(nls)
    if(ini.cov.pars[2] == 0) ini.cov.pars <- max(XY$u)/10
    if(kappa == 0) kappa <- 0.5
    Tphi.ini <- log(ini.cov.pars[2])
    XY$cov.model <- cov.model
    ##
    if (fix.nugget) {
      XY$nugget <- nugget
      if(fix.kappa){
        XY$kappa <- kappa
        res <- nls((v-nugget) ~ matrix((1-cov.spatial(u,cov.pars=c(1,exp(Tphi)), cov.model=cov.model, kappa=kappa)), ncol=1), start=list(Tphi=Tphi.ini), data=XY, alg="plinear", ...)
      }
      else{
        if(cov.model == "powered.exponential")
          res <- nls((v-nugget) ~ matrix((1-cov.spatial(u,cov.pars=c(1,exp(Tphi)), cov.model=cov.model, kappa=(2*exp(Tkappa)/(1+exp(Tkappa))))), ncol=1), start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini), data=XY, alg="plinear", ...)
        else
          res <- nls((v-nugget) ~ matrix((1-cov.spatial(u,cov.pars=c(1,exp(Tphi)), cov.model=cov.model, kappa=exp(Tkappa))), ncol=1), start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini), data=XY, alg="plinear", ...)       
        kappa <- exp(coef(res)["Tkappa"])
        names(kappa) <- NULL
      }
      cov.pars <- coef(res)[c(".lin", "Tphi")]
      names(cov.pars) <- NULL
      cov.pars[2] <- exp(cov.pars[2])
    }
    else{
      if(fix.kappa){
        XY$kappa <- kappa
        res <- nls(v ~ cbind(1,(1- cov.spatial(u, cov.pars=c(1,exp(Tphi)), cov.model = cov.model, kappa=kappa))), start=list(Tphi=Tphi.ini), alg="plinear", data=XY, ...)
      }
      else{
        if(cov.model == "powered.exponential")
          res <- nls(v ~ cbind(1, (1-cov.spatial(u, cov.pars=c(1, exp(Tphi)), cov.model = cov.model, kappa=exp(Tkappa)))), start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini), alg="plinear", data=XY, ...)
        else
          res <- nls(v ~ cbind(1, (1-cov.spatial(u, cov.pars=c(1, exp(Tphi)), cov.model = cov.model, kappa=(2*exp(Tkappa)/(1+exp(Tkappa)))))), start=list(Tphi=Tphi.ini, Tkappa = Tkappa.ini), alg="plinear", data=XY, ...)
        kappa <- exp(coef(res)["Tkappa"]);names(kappa) <- NULL
      }
      nugget <- coef(res)[".lin1"];names(nugget) <- NULL
      cov.pars <- coef(res)[c(".lin2", "Tphi")]
      names(cov.pars) <- NULL
      cov.pars[2] <- exp(cov.pars[2]) 
    }
    if(nugget < 0 | cov.pars[1] < 0){
      warning("\nvariofit: negative variance parameter found using the default option \"nls\".\n        Try another minimisation function and/or fix some of the parameters.\n")
      temp <- c(sigmasq=cov.pars[1], phi=cov.pars[2], tausq=nugget, kappa=kappa)
      print(rbind(round(temp, dig=4), status=ifelse(c(FALSE, FALSE, fix.nugget, fix.kappa), "fix", "est")))
      return(invisible())
    }
    value <- sum(resid(res)^2)
    message <- "nls does not provides convergence message"
  }
  ##
  ## minimisation using "optim" or "nlm"
  ##
  if (minimisation.function == "nlm" | minimisation.function == "optim") {
    ##
    ## Preparing lists for the minimiser
    ##
    .global.list <- list(u = XY$u, v = XY$v, n=XY$n, fix.nugget = fix.nugget,
                         nugget = nugget, fix.kappa = fix.kappa, kappa = kappa,
                         cov.model = cov.model, m.f = minimisation.function,
                         weights = weights)
    ##
    ## Preparing initial value
    ##
    ini <- ini.cov.pars
    if(fix.nugget == FALSE) ini <- c(ini, nugget)
    if(fix.kappa == FALSE) ini <- c(ini, Tkappa.ini)
    names(ini) <- NULL
    if(minimisation.function == "nlm"){
      result <- nlm(loss.vario, ini, g.l = .global.list, ...)
      result$par <- result$estimate
      result$value <- result$minimum
      result$convergence <- result$code
      if(is.R()){
        if(!is.null(get(".temp.theta", pos =1)))
          result$par <- get(".temp.theta", pos=1)
      }
      else{
        if(!is.null(get(".temp.theta", where = 1)))
          result$par <- get(".temp.theta", where = 1)
      }
    }
    else{
      if(fix.kappa == FALSE){
        if(fix.nugget) lower <- c(0, 0, -Inf)
        else lower <- c(0, 0, 0, -Inf)
      }
      else lower <- 0
      result <- optim(ini, loss.vario, method = "L-BFGS-B",
                      lower = lower, g.l = .global.list, ...)
    }
    value <- result$value
    message <- paste(minimisation.function, "convergence code:", result$convergence)
    cov.pars <- as.vector(result$par[1:2])
    if(fix.kappa == FALSE){
      if (fix.nugget)
        Tkappa <- result$par[3]
      else{
        nugget <- result$par[3]
        Tkappa <- result$par[4]
      }
      if(.global.list$cov.model == "powered.exponential")
        kappa <- 2*(exp(Tkappa))/(1+exp(Tkappa))
      else kappa <- exp(Tkappa)
    }
    else
      if(fix.nugget == FALSE)
        nugget <- result$par[3]        
  }
  ##
  ## Estimating implicity beta
  ##
  
  ##
  ## Preparing output
  ##
  estimation <- list(nugget = nugget, cov.pars = cov.pars, 
                     cov.model = cov.model, kappa = kappa, value = value, 
                     trend = vario$trend, max.dist = max(vario$u),
                     minimisation.function = minimisation.function,
                     message = message)
  estimation$weights <- weights
  if(weights == "equal") estimation$method <- "OLS"
  else estimation$method <- "WLS"
  estimation$fix.nugget <- fix.nugget
  estimation$fix.kappa <- fix.kappa
  estimation$lambda <- vario$lambda
  estimation$call <- call.fc
  class(estimation) <- "variomodel"
  return(estimation)
}


"loss.vario" <-
  function (theta, g.l) 
{
  ##
  ## Imposing constraints for nlm
  ##
  if(g.l$m.f == "nlm"){
    .temp.theta <<- NULL
    if(g.l$fix.kappa == FALSE){
      if(g.l$fix.nugget){
        theta.minimiser <- theta[1:2]
        Tkappa <- theta[3]
      }
      else{
        theta.minimiser <- theta[1:3]
        Tkappa <- theta[4]
      }
    }
    else
      theta.minimiser <- theta
    penalty <- 10000 * sum(0 - pmin(theta.minimiser, 0))
    theta <- pmax(theta.minimiser, 0)
    if(g.l$fix.kappa == FALSE)
      theta <- c(theta.minimiser, Tkappa)
    if (any(theta.minimiser < 0))
      .temp.theta <<- theta
    else penalty <- 0
  }
  else penalty <- 0
  ##
  ## reading parameters
  ##
  if(g.l$fix.kappa == FALSE){
    if (g.l$fix.nugget){
      tausq <- g.l$nugget
      Tkappa <- theta[3]
    }
    else{
      tausq <- theta[3]
      Tkappa <- theta[4]
    }
    if(g.l$cov.model == "powered.exponential")
      kappa <- 2*(exp(Tkappa))/(1+exp(Tkappa))
    else kappa <- exp(Tkappa)
  }
  else{
    kappa <- g.l$kappa
    if (g.l$fix.nugget)
      tausq <- g.l$nugget
    else
      tausq <- theta[3]
  }
  ##
  sigmasq <- theta[1]
  phi <- theta[2]
  sill.total <- sigmasq + tausq
  ##
  ## Computing values for the theoretical variogram 
  ##
  gamma <- sill.total - cov.spatial(g.l$u, cov.model = g.l$cov.model, 
                                    kappa = kappa, cov.pars = c(sigmasq, phi))
  ##
  ## Computing loss function
  ##
  if(g.l$weight == "equal")
    loss <- sum((g.l$v - gamma)^2)
  if (g.l$weights == "npairs") 
    loss <- sum(g.l$n * (g.l$v - gamma)^2)
  if (g.l$weights == "cressie") 
    loss <- sum((g.l$n/(gamma^2)) * (g.l$v - gamma)^2)
  return(loss + penalty)
}

"print.variomodel" <-
  function(x, digits = "default", ...)
{
  if(is.R() & digits == "default") digits <- max(3, getOption("digits") - 3)
  else digits <- options()$digits
  if(x$fix.nugget){
    est.pars <- c(sigmasq = x$cov.pars[1], phi=x$cov.pars[2])
    if(x$fix.kappa == FALSE)
      est.pars <- c(est.pars, kappa = x$kappa)
  }
  else{
    est.pars <- c(tausq = x$nugget, sigmasq = x$cov.pars[1], phi=x$cov.pars[2])    
    if(x$fix.kappa == FALSE)
      est.pars <- c(est.pars, kappa = x$kappa)
  }
  if(x$weights == "equal")
    cat("variofit: model parameters estimated by OLS (ordinary least squares):\n")
  else
    cat("variofit: model parameters estimated by WLS (weighted least squares):\n")
  cat(paste("covariance model is:", x$cov.model))
  if(x$cov.model == "matern" | x$cov.model == "powered.exponential" |
     x$cov.model == "cauchy" | x$cov.model == "gneiting.matern")
    if(x$fix.kappa) cat(paste(" with fixed kappa =", x$kappa)) 
  if(x$cov.model == "matern" & x$fix.kappa & x$kappa == 0.5)
    cat(" (exponential)")
  cat("\n")
  if(x$fix.nugget)
    cat(paste("fixed value for tausq = ", x$nugget,"\n"))
  cat("parameter estimates:\n")
  print(round(est.pars, digits=digits))
  if(x$weights == "equal")
    cat("\nvariofit: minimised sum of squares = ")
  else
      cat("\nvariofit: minimised weighted sum of squares = ")
  cat(round(x$value, digits=digits))
  cat("\n")
  return(invisible())
}  

"summary.variomodel" <-
  function(object, ...)
{
  summ.lik <- list()
  if(object$weights == "equal")
    summ.lik$pmethod <- "OLS (ordinary least squares)"
  else
    summ.lik$pmethod <- "WLS (weighted least squares)"
  summ.lik$cov.model <- object$cov.model
  summ.lik$spatial.component <- c(sigmasq = object$cov.pars[1], phi=object$cov.pars[2])
  summ.lik$spatial.component.extra <- c(kappa = object$kappa)
  summ.lik$nugget.component <- c(tausq = object$nugget)
  summ.lik$fix.nugget <- object$fix.nugget
  summ.lik$fix.kappa <- object$fix.kappa
  summ.lik$sum.of.squares <- c(value = object$value)
  if(object$fix.nugget){
    summ.lik$estimated.pars <- c(sigmasq = object$cov.pars[1], phi=object$cov.pars[2])
    if(object$fix.kappa == FALSE)
      summ.lik$estimated.pars <- c(summ.lik$estimated.pars, kappa = object$kappa)
  }
  else{
    summ.lik$estimated.pars <- c(tausq = object$nugget, sigmasq = object$cov.pars[1], phi=object$cov.pars[2])
    if(object$fix.kappa == FALSE)
      summ.lik$estimated.pars <- c(summ.lik$estimated.pars, kappa = object$kappa)
  }
  summ.lik$weights <- object$weights
  summ.lik$call <- object$call
  class(summ.lik) <- "summary.variomodel"
  return(summ.lik)
}

"print.summary.variomodel" <-
  function(x, digits = "default", ...)
{
  if(class(x) != "summary.variomodel")
    stop("object is not of the class \"summary.variomodel\"")
  if(is.R() & digits == "default") digits <- max(3, getOption("digits") - 3)
  else digits <- options()$digits
  cat("Summary of the parameter estimation\n")
  cat("-----------------------------------\n")
  cat(paste("Estimation method:", x$pmethod, "\n"))
  cat("\n")
  ##
  ## Estimates of the model components
  ## Model: Y(x) = X\beta + S(x) + e 
  ##
#  cat("Parameters of the mean component (trend):")
#  cat("\n")
#  print(round(x$mean.component, digits=digits))
#  cat("\n")
  ##
  cat("Parameters of the spatial component:")
  cat("\n")
  cat(paste("   correlation function:", x$cov.model))
  if(x$cov.model == "matern" & x$fix.kappa & x$spatial.component.extra == 0.5)
    cat(" (exponential)")
  if(x$cov.model == "matern" | x$cov.model == "powered.exponential" |
     x$cov.model == "cauchy" | x$cov.model == "gneiting.matern"){
    if(x$fix.kappa)
      cat(paste("\n      (fixed) extra parameter kappa = ", round(x$spatial.component.extra, digits=digits)))
    else
      cat(paste("\n      (estimated) extra parameter kappa = ", round(x$spatial.component.extra, digits=digits)))
  }
  cat(paste("\n      (estimated) variance parameter sigmasq (partial sill) = ", round(x$spatial.component[1], dig=digits)))
  cat(paste("\n      (estimated) cor. fct. parameter phi (range parameter)  = ", round(x$spatial.component[2], dig=digits)))
  cat("\n")
  ##
  cat("\n")  
  cat("Parameter of the error component:")
  if(x$fix.nugget)
    cat(paste("\n      (fixed) nugget =", round(x$nugget.component, digits = digits)))
  else
    cat(paste("\n      (estimated) nugget = ", round(x$nugget.component, dig=digits)))
  cat("\n")
  cat("\n")
  names(x$sum.of.squares) <- NULL
  if(x$weights == "equal") cat("Minimised sum of squares: ")
  else cat("Minimised weighted sum of squares: ")
  cat(round(x$sum.of.squares, digits=digits))
  cat("\n")
  cat("\n")
  cat("Call:")
  cat("\n")
  print(x$call)
  cat("\n")
  invisible(x)
}

##"beta.variofit" <-
##  function(geodata, coords = geodata$coords, data=geodata$data,
##           obj.variofit)
##  {
##
##  }
