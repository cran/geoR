"likfit.old" <-
  function (geodata, coords=geodata$coords, data=geodata$data, trend = "cte",
            ini, fix.nugget = FALSE, nugget = 0, 
            cov.model = "matern",
            kappa = 0.5, fix.lambda = TRUE, lambda = 1, method = "ML", 
            predicted = FALSE, residuals = FALSE, 
            minimisation.function = c("optim","nlmP", "nlm"),
            automatic.refit = FALSE, range.limits,
            messages.screen = TRUE, ...) 
{
  call.fc <- match.call()
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (cov.model=="pure.nugget"){
    if(fix.nugget == TRUE) ini <- rep(0,2)
    else
      if(fix.nugget == TRUE) ini <- rep(0,3)
  }
  if(!is.null(kappa))
    if(cov.model == "matern" & kappa == 0.5)
      cov.model <- "exponential"
  minimisation.function <- match.arg(minimisation.function)
  if(is.R()) require(mva)
  ftau <- nugget
  fixtau <- fix.nugget
  coords <- as.matrix(coords)
  dists.vec <- as.vector(dist(coords))
  range.dist <- range(dists.vec)
  max.dist <- max(range.dist)
  min.dist <- min(range.dist)
  if(missing(range.limits)){
    lower.phi <- 0
    upper.phi <- +Inf
  }
  else{
    lower.phi <- range.limits[1]
    upper.phi <- range.limits[2]
  }
  z <- as.vector(data)
  if(fix.lambda) {
    if(lambda == 1) {
      log.jacobian <- 0
    }
    else {
      if(any(z <= 0))
        stop("Transformation option not allowed when there are zeros or negative data"
             )
      if(any(z^(lambda - 1) <= 0))
        log.jacobian <- log(prod(z^(lambda - 1)))
      else log.jacobian <- sum(log(z^(lambda - 1)))
      if(lambda == 0)
        z <- log(z)
      else z <- ((z^lambda) - 1)/lambda
    }
  }
  n <- length(z)
  if ((2*n) != length(coords))
    stop("Number of locations does not match with number of data")
  reduce.pars <- 0
  if (method == "REML" | method == "reml" | method == "rml") 
    method <- "RML"
  if(method == "ML" | method == "ml")
    method <- "ML"
  if(method == "ML" & cov.model == "power")
    stop("\n\"power\" model can only be used with method=\"RML\".\nBe sure that what you want is not \"powered.exponential\"")
  xmat <- trend.spatial(trend=trend, coords=coords)
  fit.ols <- lm(z ~ xmat + 0)
  trend.ols <- list(coefficients = fit.ols$coefficients)
  var.z <- sum((fit.ols$residual)^2)/(n-length(fit.ols$coefficients))
  dimnames(xmat) <- list(NULL, NULL)
  txmat <- t(xmat)
  beta.size <- dim(xmat)[2]  
  if(missing(ini) | ini=="default"){
    cat("likfit: no initial values for the parameters was provided. Default initial values will be used\n")
    if(fix.nugget==FALSE) ini <- c(.2*var.z, 0.8*var.z, max.dist/5)
    else ini <- c(0.8*var.z, max.dist/5)
  }    
  if(all(ini==0)){
    cov.model <- "pure.nugget"
    cat("likfit: all initial values equal to zero. Model without spatial correlation will be fitted\n")
  }
  else{
    if(is.matrix(ini)) {
      inilength <- dim(ini)[2]
      if(fixtau == FALSE & inilength != 3)
        stop("wrong number of columns for ini (must be 3)")
      if(fixtau == TRUE & inilength != 2)
        stop("wrong number of columns for ini (must be 2)")
    }
    else {
      inilength <- length(ini)
      if (fixtau == FALSE & inilength != 3) 
        stop("wrong length for ini (must be 3)")
      if (fixtau == TRUE & inilength != 2) 
        stop("wrong length for ini (must be 2)")
    }
  }
  assign(".temp.list", list(z = z, xmat = xmat,  txmat = txmat, fixtau = fixtau, 
                            ftau = ftau, method = method, kappa = kappa,
                            cov.model = cov.model, beta.size = beta.size,
                            lower.phi = lower.phi, 
                            dists.lowertri = dists.vec, var.z = var.z,
                            fix.lambda = fix.lambda, n = n,
                            minimisation.function=minimisation.function), pos=1)
  if(fix.lambda == TRUE) {
    .temp.list$lambda <<- lambda
    .temp.list$log.jacobian <<- log.jacobian
  }
  if ((cov.model == "pure.nugget") | all(ini==0) ){  
    if(messages.screen == TRUE)
      cat("likfit: fitting model without spatial correlation\n")
    lik.results <- likfit.nospatial(.temp.list, ...)
    if (fix.nugget == FALSE)
      temp.pars <- c(lik.results$tausqhat, 0, 0)
    else
      temp.pars <- c(ftau, (lik.results$tausqhat - ftau), 0)        
    lambda <- lik.results$lambda
  }
  else{
    if(is.matrix(ini) | is.data.frame(ini)) {
      ini <- as.matrix(ini)
      if(messages.screen == TRUE)
        cat("likfit: searching for the best initial value\n")
      ini.search <- ini
      if(fix.nugget == TRUE)
        ini.search <- cbind(nugget, ini.search)
      if(length(lambda) == 1)
        ini.search <- cbind(ini.search, lambda)
      else {
        temp <- ini.search
        for(i in 1:(length(lambda) - 1)) {
          ini.search <- rbind(ini.search, temp)
        }
        ini.search <- cbind(ini.search, rep(lambda, each = dim(
                                                      temp)[1]))
      }
      dimnames(ini.search) <- list(NULL, NULL)
      loglik.ini <- round(100000000. * apply(ini.search, 1, 
                                             loglik.spatial))
      ini.max <- as.vector(ini.search[loglik.ini == max(loglik.ini),
                                      ])
      if(fixtau == TRUE) {
        ini <- as.vector(ini.max[2:3])
        if(minimisation.function == "nlmP" & ini[2] == 0)
          ini[2] <- min(ini.search[(ini.search[,3] != 0),3])
      }
      else {
        ini <- as.vector(ini.max[1:3])
        if(minimisation.function == "nlmP" & ini[3] == 0)
          ini[3] <- min(ini.search[(ini.search[,3] != 0),3])
        
      }
      if(messages.screen == TRUE) {
        cat("likfit: best initial value:\n")
        names(ini.max) <- c("nugget", "sill", "range", "lambda"
                            )
        print(ini.max)
      }
      if(fix.lambda == FALSE)
        lambda <- as.vector(ini.max[4])
    }
    if (messages.screen == TRUE) {
      cat(paste("likfit: Initialising likelihood maximisation using the function", minimisation.function, "\n"))
      cat("------------------------------------------------------------\n")
      cat("likfit: consider providing additional (non-default) arguments for the minimisation function.\n")
      if(minimisation.function == "nlm" | minimisation.function == "nlmP"){
        cat("likfit: some relevant arguments are: iterlim, steptol, stepmax, print.level, ndigit. For more details see documentation for the function nlm.\n")
      }
      if(minimisation.function == "optim"){
      cat("likfit: parameters can be passed to the function optim using the argument control(). For more details see documentation for the function optim.\n")
    }
      cat("likfit: it is highly advisable to run the function several times with different initial values for the parameters (argument ini).\n")
      cat("------------------------------------------------------------\n")
      "nice" <-
        function (x, decimal = 2, fixed = FALSE) 
          {
            ergb <- x
            index <- (x != 0) & is.finite(x)
            if (fixed) 
              n <- 0
            else n <- floor(log(abs(x[index]))/log(10))
            ergb[index] <- trunc(x[index]/10^(n - decimal)) * 10^(n - decimal)
            return(ergb)
          }
      cat(paste("likfit: range of values allowed for the parameter:", nice(lower.phi), "to", nice(upper.phi), "\n"))
      cat("likfit: WARNING: This step can be time demanding!\n")
      cat("\n")
    }
    if (fixtau == FALSE | ftau != 0) {
      if (fixtau == TRUE & ftau != 0) {
                                        #        if (messages.screen == TRUE)
                                        #          print("covariance parameters used in the minimization function are $\sigma^2$ and $\phi$")
        if(minimisation.function == "nlm") assign(".temp.lower", c(0, lower.phi), pos=1)
        if(fix.lambda == TRUE) {
          assign(".temp.lower", c(0, lower.phi), pos=1)
          if(minimisation.function == "nlm"){
            lik.results <- nlm(proflik.ftau, ini, ...)
            if(exists(".temp.sill")){
              lik.results$estimate[1] <- .temp.sill
              remove(".temp.sill", pos=1)
            }
            if(exists(".temp.phi")){
              lik.results$estimate[2] <- .temp.phi
              remove(".temp.phi", pos=1)
            }
            rm(.temp.lower, inherits = TRUE, pos=1)
          }
          if(minimisation.function == "nlmP"){
            assign(".ind.prof.phi", 2, pos=1)
            lik.results <- nlmP(proflik.ftau, ini, lower=c(0, lower.phi), upper=c(10000*var.z, upper.phi), ...)
          }            
          if(minimisation.function == "optim"){
            lik.results <- optim(ini, proflik.ftau, method="L-BFGS-B", lower=c(0, lower.phi), upper=c(10000*var.z, upper.phi), ...)
            lik.results$estimate <- lik.results$par
          }            
        }
        else{
          if(minimisation.function == "nlm"){
            assign(".temp.lower", c(0, lower.phi), pos=1)
            assign(".temp.lower.lambda", -2, pos=1)
            assign(".temp.upper.lambda", 2, pos=1)
            lik.results <- nlm(proflik.ftau, c(ini,lambda), ...)
            if(exists(".temp.sill")){
              lik.results$estimate[1] <- .temp.sill
              remove(".temp.sill", pos=1)
            }
            if(exists(".temp.phi")){
              lik.results$estimate[2] <- .temp.phi
              remove(".temp.phi", pos=1)
            }
            if(exists(".temp.lambda")){
              lambda <- .temp.lambda
              remove(".temp.lambda", pos=1)
            }
            else{
              lambda <- lik.results$estimate[3]
            }
            rm(.temp.lower, .temp.lower.lambda, .temp.upper.lambda, inherits = TRUE, pos=1)
          }
          if(minimisation.function == "nlmP"){
            assign(".ind.prof.phi", 2, pos=1)
            lik.results <- nlmP(proflik.ftau, c(ini,lambda), lower = c(0, lower.phi, -2), upper = c(10000*var.z, upper.phi, 2), ...)
            lambda <- lik.results$estimate[3]
          }
          if(minimisation.function == "optim"){
            lik.results <- optim(c(ini,lambda), proflik.ftau, method="L-BFGS-B", lower = c(0, lower.phi, -2), upper = c(10000*var.z, upper.phi, 2), ...)
            lik.results$estimate <- lik.results$par              
            lambda <- lik.results$estimate[3]
          }
          lik.results$estimate <- as.vector(lik.results$estimate[1:2])
          if(lambda == 0)
            z <- log(as.vector(data))
          else z <- (((as.vector(data))^lambda) - 1)/
            lambda
        }        
        lik.results$estimate <- temp.pars <- as.vector(c(ftau, lik.results$estimate))
      }
      if (fixtau == FALSE) {
                                        #        if (messages.screen == TRUE) 
                                        #          print("parameters used in the minimization function are the ratio (tau^2/sigma^2) and $\phi$")
        ini.m <- c(ini[1]/ini[2], ini[3])
        if(fix.lambda == TRUE) {
          if (minimisation.function=="nlm"){
            assign(".temp.lower", c(0, lower.phi), pos=1)
            lik.results <- nlm(proflik.nug, ini.m, ...) 
            if(exists(".temp.nugget")){
              lik.results$estimate[1] <- .temp.nugget
              remove(".temp.nugget", pos=1)
            }
            if(exists(".temp.phi")){
              lik.results$estimate[2] <- .temp.phi
              remove(".temp.phi", pos=1)
            }
            rm(.temp.lower, inherits = TRUE, pos=1)
          }
          if (minimisation.function=="nlmP"){
            if(ini.m[1] == 0) ini.m[1] <- 0.05
            assign(".ind.prof.phi", 2, pos=1)
            lik.results <- nlmP(proflik.nug, ini.m, lower=c(0, lower.phi), upper=c(100, upper.phi),...) 
          }
          if (minimisation.function=="optim"){
            lik.results <- optim(ini.m, proflik.nug, method="L-BFGS-B", lower=c(0, lower.phi), upper=c(100, upper.phi),...) 
            lik.results$estimate <- lik.results$par              
          }
        }
        else{
          if (minimisation.function=="nlm"){
            assign(".temp.lower", c(0, lower.phi), pos=1)
            assign(".temp.lower.lambda", -2, pos=1)
            assign(".temp.upper.lambda", 2, pos=1)
            lik.results <- nlm(proflik.nug, c(ini.m, lambda), ...)
            if(exists(".temp.nugget")){
              lik.results$estimate[1] <- .temp.nugget
              remove(".temp.nugget", pos=1)
            }
            if(exists(".temp.phi")){
              lik.results$estimate[2] <- .temp.phi
              remove(".temp.phi", pos=1)
            }
            if(exists(".temp.lambda")){
              lambda <- .temp.lambda
              remove(".temp.lambda", pos=1)
            }
            else{
              lambda <- lik.results$estimate[3]
            }
            rm(.temp.lower, .temp.lower.lambda,  .temp.upper.lambda, inherits = TRUE, pos=1)
          }
          if (minimisation.function=="nlmP"){
            assign(".ind.prof.phi", 2, pos=1)
            lik.results <- nlmP(proflik.nug, c(ini.m, lambda), lower=c(0, lower.phi, -2), upper=c(100, upper.phi, 2),...)
            lambda <- lik.results$estimate[3]
          }
          if (minimisation.function=="optim"){
            lik.results <- optim(c(ini.m, lambda), proflik.nug, method="L-BFGS-B", lower=c(0, lower.phi, -2), upper=c(100, upper.phi, 2),...)
            lik.results$estimate <- lik.results$par              
            lambda <- lik.results$estimate[3]
          }            
          lik.results$estimate <- as.vector(lik.results$estimate[1:2])
          if(lambda == 0)
            z <- log(as.vector(data))
          else z <- (((as.vector(data))^lambda) - 1)/
            lambda
        }
        if(messages.screen == TRUE) {
          if(minimisation.function == "nlm" | minimisation.function == "nlmP") 
          if(minimisation.function == "optim") cat(paste("likfit: optim convergence code: ",lik.results$convergence, "\n"))
        }
        if(automatic.refit == TRUE & (lik.results$estimate[1] < 0.01)) {
          if (messages.screen == TRUE)
            cat(paste("likfit: WARNING: ratio of estimates tau^2/sigma^2 < 0.01 (",round(lik.results$estimate[1], dig = 4), ")", sep = ""))
          cat("\n")
          reduce.pars <- 1
          .temp.list$ftau <<- 0
          .temp.list$fixtau <<- T
          if(fix.lambda == TRUE) {
            if (minimisation.function=="nlm"){
              assign(".temp.lower.phi", lower.phi, pos=1)
              lik.results <- nlm(proflik.phi, ini[3],  ...)
              if(exists(".temp.phi")){
                lik.results$estimate <- .temp.phi
                remove(".temp.phi", pos=1)
              }
              rm(.temp.lower, inherits = TRUE, pos=1)
            }
            if (minimisation.function=="nlmP"){
              assign(".ind.prof.phi", 1, pos=1)
              lik.results <- nlmP(proflik.phi, ini[3],  lower=lower.phi, upper=upper.phi,...)
            }
            if (minimisation.function=="optim"){
              lik.results <- optim(ini[3], proflik.phi, method="L-BFGS-B",  lower=lower.phi, upper=upper.phi,...)
              lik.results$estimate <- lik.results$par  
            }
          }
          else {
            if (minimisation.function=="nlm"){
              assign(".temp.lower.phi", lower.phi, pos=1)
              assign(".temp.lower.lambda", -2, pos=1)
              assign(".temp.upper.lambda", 2, pos=1)
              lik.results <- nlm(proflik.phi, c(ini[3], lambda), ...)
              if(exists(".temp.lambda")){
                lambda <- .temp.lambda
                remove(".temp.lambda", pos=1)
              }
              else{
                lambda <- lik.results$estimate[2]
              }
              if(exists(".temp.phi")){
                lik.results$estimate <- .temp.phi
                remove(".temp.phi", pos=1)
              }
              else{
                lik.results$estimate <- as.vector(lik.results$estimate[1])
              }
              rm(.temp.lower.phi, .temp.lower.lambda, .temp.upper.lambda, inherits = TRUE, pos=1)
            }
            if (minimisation.function=="nlmP"){
              lik.results <- nlmP(proflik.phi, c(ini[3], lambda), lower=c(lower.phi, -2), upper=c(upper.phi, 2),...)
              lambda <- lik.results$estimate[2]
              lik.results$estimate <- lik.results$estimate[1]
            }
            if (minimisation.function=="optim"){
              lik.results <- optim(c(ini[3], lambda), proflik.phi, method="L-BFGS-B", lower=c(lower.phi, -2), upper=c(upper.phi, 2),...)
              lik.results$estimate <- lik.results$par  
              lambda <- lik.results$estimate[2]
              lik.results$estimate <- lik.results$estimate[1]
            }
            if(lambda == 0)
              z <- log(as.vector(data))
            else z <- (((as.vector(data))^lambda) -
                       1)/lambda
          }
          if (messages.screen == TRUE)        
            cat("likfit: model re-fitted without nugget effect (tausq = 0)\n")
          lik.results$estimate <- as.vector(c(0, lik.results$estimate))
          if(messages.screen == TRUE) {
            if(minimisation.function == "nlm" | minimisation.function == "nlmP") cat(paste("likfit: nlm optimisation code: ",lik.results$code,"\n"))
            if(minimisation.function == "optim") cat(paste("likfit: optim convergence code: ",lik.results$convergence,"\n"))
          }
        }          
        nugget.rel <- lik.results$estimate[1]
        if (lik.results$estimate[2] < 1e-08)
          icovhat <- diag(n)
        else
          icovhat <- varcov.spatial(coords = coords, cov.model = 
                                    cov.model, kappa = kappa, nugget = nugget.rel,
                                    cov.pars = c(1, lik.results$estimate[
                                      2]), inv = TRUE, det = FALSE)$inverse
        txiv <- crossprod(xmat, icovhat)
        sigmasqhat <- (z %*% (icovhat - crossprod(txiv,solve(txiv %*% xmat)) %*% txiv) %*% z)/n
        if(method == "RML") sigmasqhat <- sigmasqhat * n/(n-beta.size)
        nuggethat <- lik.results$estimate[1] * sigmasqhat
        lik.results$estimate <- temp.pars <- as.vector(c(nuggethat, sigmasqhat, lik.results$estimate[2]))
      }
      lik.results$estimate <- as.vector(lik.results$estimate)
      if((automatic.refit == TRUE & (lik.results$estimate[3] <= lower.phi)) | lik.results$estimate[3] < 1e-12) {
        if (messages.screen == TRUE){
          cat("likfit: WARNING: phi estimate < minimum value allowed\n")
          cat("likfit: model re-fitted without spatial correlation (phi=0)\n")
        }
        reduce.pars <- 2
        lik.results <- likfit.nospatial(.temp.list, ...)
        lambda <- lik.results$lambda
        if(fix.nugget == TRUE) {
          lik.results$parameters <- temp.pars <-
            as.vector(c(ftau, (lik.results$tausqhat - ftau), 0))
        }
        else {
          lik.results$parameters <- temp.pars <-
            as.vector(c(lik.results$tausqhat, 0, 0))
        }
      }
    }
    else {
                                        # case 3: parameters = $(\sigma^2, \phi)$ ; fixed nugget: tau^2= 0$
      ini.m <- ini[2]
      if(fix.lambda == TRUE) {
        if (minimisation.function=="nlm"){
          assign(".temp.lower.phi", lower.phi, pos=1)
          lik.results <- nlm(proflik.phi,ini.m,   ...)
          if(exists(".temp.phi")){
            lik.results$estimate <- .temp.phi
            remove(".temp.phi", pos=1)
          }
          rm(.temp.lower, inherits = TRUE, pos=1)
        }
        if (minimisation.function=="nlmP"){
          assign(".ind.prof.phi", 1, pos=1)
          lik.results <- nlmP(proflik.phi,ini.m, lower=lower.phi, upper=upper.phi,...)
        }
        if (minimisation.function=="optim"){
          lik.results <- optim(ini.m, proflik.phi, method="L-BFGS-B", lower=lower.phi, upper=upper.phi,...)
          lik.results$estimate <- lik.results$par
        }
      }
      else {
        if (minimisation.function=="nlm"){
          assign(".temp.lower.phi", lower.phi, pos=1)
          assign(".temp.lower.lambda", -2, pos=1)
          assign(".temp.upper.lambda", 2, pos=1)
          lik.results <- nlm(proflik.phi, c(ini.m, lambda), ...)
          if(exists(".temp.lambda")){
            lambda <- .temp.lambda
            remove(".temp.lambda", pos=1)
          }
          else{
            lambda <- lik.results$estimate[2]
          }
          if(exists(".temp.phi")){
            lik.results$estimate <- .temp.phi
            remove(".temp.phi", pos=1)
          }
          else{
            lik.results$estimate <- as.vector(lik.results$estimate[1])
          }
          rm(.temp.lower.phi, .temp.lower.lambda, .temp.upper.lambda, inherits = TRUE, pos=1)
        }
        if (minimisation.function=="nlmP"){
          assign(".ind.prof.phi", 1, pos=1)
          lik.results <- nlmP(proflik.phi, c(ini.m, lambda), lower=c(lower.phi, -2), upper=c(upper.phi, 2),...)
          lambda <- as.vector(lik.results$estimate[2])
          lik.results$estimate <- as.vector(lik.results$estimate[1])
        }
        if (minimisation.function=="optim"){
          lik.results <- optim(c(ini.m, lambda), proflik.phi, method="L-BFGS-B", lower=c(lower.phi, -2.5), upper=c(upper.phi, 2.5),...)
          lik.results$estimate <- lik.results$par        
          lambda <- as.vector(lik.results$estimate[2])
          lik.results$estimate <- as.vector(lik.results$estimate[1])
        }
        if(lambda == 0)
          z <- log(as.vector(data))
        else z <- (((as.vector(data))^lambda) - 1)/lambda
      }    
      if(messages.screen == TRUE) {
        if(minimisation.function == "nlm" | minimisation.function == "nlmP")
          cat(paste("likfit: nlm optimisation code: ",lik.results$code, "\n"))
        if(minimisation.function == "optim") cat(paste("likfit: optim convergence code: ",lik.results$convergence, "\n"))
      }      
      if(automatic.refit == TRUE & (lik.results$estimate <= lower.phi)) {
        if (messages.screen == TRUE) {
          cat("likfit: WARNING: phi estimate < minimum value allowed\n")
          cat("likfit: model without spatial correlation was fitted (phi=0 and sigma^2=0)\n")
        }
        reduce.pars <- 1
        lik.results <- likfit.nospatial(.temp.list, ...)
        lambda <- lik.results$lambda
        if(fix.nugget == TRUE) {
          lik.results$parameters <- temp.pars <-
            as.vector(c(ftau, (lik.results$tausqhat - ftau), 0))
        }
        else {
          lik.results$parameters <- temp.pars <-
            as.vector(c(lik.results$tausqhat, 0, 0))
        }
      }
      else {
        if(lik.results$estimate < 1e-08)
          icovhat <- diag(n)
        else
          icovhat <- varcov.spatial(coords = coords, cov.model = 
                                    cov.model, kappa = kappa,
                                    nugget = 0, cov.pars
                                    = c(1, lik.results$estimate),
                                    inv = TRUE, det = FALSE)$inverse
        txiv <- crossprod(xmat, icovhat)
        sigmasqhat <- (z %*% (icovhat - crossprod(txiv, solve(txiv %*% xmat
                                                              ) %*% txiv)) %*% z)/n
        if(method == "RML") sigmasqhat <- sigmasqhat * n/(n-beta.size)
        temp.pars <- as.vector(c(0, sigmasqhat, lik.results$estimate))
        lik.results$estimate <- as.vector(c(0,sigmasqhat, lik.results$estimate))
      }
    }
  }
  if(messages.screen == TRUE) {
    cat("likfit: end of likelihood maximisation\n")
  }
  if(any(temp.pars < 0)){
    temp.pars <- round(temp.pars, dig=14)
    lik.results$estimate <- round(lik.results$estimate, dig=14)
  }
  if(minimisation.function == "optim") lik.results$minimum <- lik.results$value
  loglik <- -lik.results$minimum
  npars <- length(trend.ols$coefficients) + length(ini) - reduce.pars
  if (fix.lambda == FALSE) npars <- npars + 1
  AIC <- loglik - npars
  BIC <- loglik - 0.5 * log(n) * npars
  if (messages.screen == TRUE) 
    cat("likfit: computing the beta estimate\n")
  if(any(temp.pars[2:3]) != 0)
    invcov <- varcov.spatial(coords = coords, cov.model = cov.model, 
                             kappa = kappa, nugget = temp.pars[1],
                             cov.pars = temp.pars[2:3], 
                             inv = TRUE, det = FALSE)$inverse
  else invcov <- diag((1/temp.pars[1]), n)
  txmatinvcov <- crossprod(xmat, invcov)
  beta <- solve(txmatinvcov %*% xmat) %*% txmatinvcov %*% z
  beta.var <- solve(txmatinvcov %*% xmat)
  if (residuals == TRUE | predicted == TRUE) {
    cat("likfit: computing predicted values and residuals\n")
    trend.est <- as.vector(xmat %*% beta)
    residuals.trend <- as.vector(z - trend.est)
    covmat.signal <- varcov.spatial(coords = coords, cov.model = cov.model, 
                                    kappa = kappa, nugget = 0,
                                    cov.pars = temp.pars[2:3])$varcov
    signal.est <- as.vector(covmat.signal %*% invcov %*% 
                            residuals.trend)
    predict.est <- trend.est + signal.est
    residuals.est <- as.vector(z - predict.est)
    residuals.std <- as.vector(invcov %*% residuals.est)
    residuals.trend.std <- as.vector(invcov %*% residuals.trend)
    s2.trend <- (crossprod(residuals.trend,invcov) %*% residuals.trend)/(n - 
                                                                         length(beta))
    s2 <- (crossprod(residuals.est,invcov) %*% residuals.est)/(n - 
                                                               length(beta))
  }
  if (messages.screen == TRUE) 
    cat("likfit: preparing output\n")
  results <- list()
  results$cov.model <- cov.model
  results$nugget <- temp.pars[1]
  results$cov.pars <- as.vector(c(sigmasq = temp.pars[2], phi = temp.pars[3]))
  if (is.null(kappa))
    results$kappa <- "not used"
  else
    results$kappa <- kappa
  results$beta <- as.vector(beta)
  results$beta.var <- beta.var
  if (length(results$beta.var) == 1)
    results$beta.var <- as.vector(results$beta.var)
  if (length(results$beta) > 1){
    if(inherits(trend, "formula"))
      beta.names <- c("intercept", paste("covar", 1:(ncol(xmat)-1), sep = ""))
    else
      if (trend == "1st")
        beta.names <- c("1", "x", "y")
      else
        if (trend == "2nd")
          beta.names <- c("1", "x", "y", "x2", "xy", "y2")
    names(results$beta) <- beta.names
  }
  results$lambda <- lambda
  results$loglik <- loglik
  results$npars <- npars
  results$AIC <- AIC
  results$BIC <- BIC
  results$trend.ols <- as.vector(trend.ols$coefficients)
  names(results$trend.ols) <- names(results$beta)
  if (residuals == TRUE) {
    results$s2 <- s2
    results$s2.trend <- s2.trend
  }
  if (predicted == TRUE) 
    results$predicted <- cbind(predicted = predict.est, trend.est = trend.est, 
                               signal.est = signal.est)
  if (residuals == TRUE) 
    results$residuals <- round(cbind(residuals = residuals.est, 
                                     resid.trend = residuals.trend, resid.std = residuals.std, 
                                     resid.trend.std = residuals.trend.std), dig = 12)
  if(fix.lambda == FALSE) {
    if(lambda == 1) {
      log.jacobian <- 0
    }
    else {
      if(any(data^(lambda - 1) <= 0))
        log.jacobian <- log(prod(data^(lambda - 1)))
      else log.jacobian <- sum(log(data^(lambda - 1)))
    }
  }
  results$info.lambda <- list(fix.lambda = fix.lambda, log.jacobian = 
                              log.jacobian)
  lik.results$estimate <- NULL
  lik.results$aux <- NULL
  lik.results$minimum <- NULL
  results$method <- method
  results$info <- lik.results
  results$max.dist <- max.dist
  results$trend.matrix <- xmat
  results$call <- call.fc
  class(results) <- "variomodel"
  if(messages.screen == TRUE){
    cat("likfit: estimated model parameters are:\n")
    cat(paste("covariance model:", cov.model))
    if(cov.model == "matern" | cov.model == "powered.exponential" | 
       cov.model == "cauchy" | cov.model == "gneiting.matern")
      cat(paste(" with kappa =", kappa))
    if(!is.null(kappa))
      if(cov.model == "matern" & kappa == 0.5)
        cat(" (exponential)")
    cat("\n")
    print(c(nugget=results$nugget, sill=results$cov.pars[1], range=results$cov.pars[2]))
    if (fix.lambda == FALSE)
      cat(paste("Box-Cox transformation parameter:", round(results$lambda, dig=4),"\n"))
    if((results$cov.pars[1] < (0.01 * (results$nugget + results$cov.pars[1])))& results$cov.pars[2] > 0)
      cat("\nWARNING: estimated sill is less than 1 hundredth of the total variance. Consider re-examine the model excluding spatial dependence\n" )      
    if((results$cov.pars[2] > (10 * max.dist)) & results$cov.pars[1] > 0 )
      cat("\nWARNING: estimated range is more than 10 times bigger than the biggest distance between two points. Consider re-examine the model:\n 1) excluding spatial dependence if estimated sill is too low and/or \n 2) taking trends (covariates) into account\n" ) 
    if(((results$cov.pars[2] < (0.1 * min.dist)) & (results$cov.pars[1] > 0)) & results$cov.pars[2] > 0)
      cat("\nWARNING: estimated range is less than 1 tenth of the minimum distance between two points. Consider re-examine the model excluding spatial dependence\n" ) 
  }
  remove(".temp.list", pos=1)
  return(results)
}

"proflik.ftau" <-
  function (theta) 
{
  if (any(is.na(theta)) | any(theta==Inf) | any(is.nan(theta)))
    neglik <- 1e+32
  else{
    if(length(theta) == 3) include.lambda <- TRUE else include.lambda <- FALSE 
    if(.temp.list$minimisation.function == "nlm"){
      if (exists(".temp.phi", w=1)) remove(".temp.phi", pos=1, inherits = TRUE)
      if (exists(".temp.lambda", w=1)) remove(".temp.lambda", pos=1, inherits = TRUE)
      if (exists(".temp.sill", w=1)) remove(".temp.sill", pos=1, inherits = TRUE)
      theta.minimiser <- theta
      penalty <- 10000 * sum(.temp.lower - pmin(theta[1:2], .temp.lower))
      theta[1:2] <- pmax(theta[1:2], .temp.lower)
      if (theta.minimiser[1] <  .temp.lower[1])
        assign(".temp.sill", theta[1], pos=1)
      if (theta.minimiser[2] < 1.001 * .temp.lower[2])
        assign(".temp.phi", theta[2], pos=1)
      if (include.lambda){
        lambda <- theta[3]
        penalty <- penalty + 1000 * (.temp.lower.lambda - min(lambda, .temp.lower.lambda))
        lambda <- max(lambda, .temp.lower.lambda)
        penalty <- penalty + 1000 * (.temp.upper.lambda - max(lambda, .temp.upper.lambda))
        lambda <- min(lambda, .temp.upper.lambda)
        if (round(1000 * theta.minimiser[3]) <= round(1000 * .temp.lower.lambda))
          assign(".temp.lambda", lambda, pos=1)
        if (round(1000 * theta.minimiser[3]) >= round(1000 * .temp.upper.lambda))
          assign(".temp.lambda", lambda, pos=1)
      }
    }
    else{
      if (include.lambda) lambda <- theta[3]
    }
    z <- .temp.list$z
    n <- length(z)
    if (include.lambda){
      if(lambda == 1) {
        .temp.list$log.jacobian <<- 0
      }
      else {
        if(any(z < 0))
          stop("Transformation option not allowed when there are zeros or negative data"
               )
        if(any(z^(lambda - 1) <= 0))
          .temp.list$log.jacobian <<- log(prod(z^(lambda - 1)))
        else .temp.list$log.jacobian <<- sum(log(z^(lambda - 1)))
        if(lambda == 0)
          z <- log(z)
        else z <- ((z^lambda) - 1)/lambda
      }
    }
    beta.size <- .temp.list$beta.size
    kappa <- .temp.list$kappa
    ftau <- .temp.list$ftau
    sigmasq <- theta[1]
    sill.total <- ftau + sigmasq
    phi <- theta[2]
    covinf <- varcov.spatial(dists.lowertri = .temp.list$dists.lowertri,
                             cov.model = .temp.list$cov.model, kappa = kappa, 
                             nugget = ftau, cov.pars = c(sigmasq, phi), 
                             det = TRUE, func.inv = "eigen",
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
    ssresmat <- as.vector(yiy - crossprod(beta.hat,xiy))
    if(.temp.list$method == "ML") {
      neglik <- ((n/2) * log(2 * pi) +
                 covinf$log.det.to.half +
                 0.5 * ssresmat -
                 .temp.list$log.jacobian
                 )
    }
    if(.temp.list$method == "RML") {
      xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE, only.values = TRUE)
      neglik <- (((n - beta.size)/2) * log(2 * pi) +
                 covinf$log.det.to.half +
                 0.5 * ssresmat +
                 choldet -
                 0.5 * sum(log(xx.eigen$values)) -
                 .temp.list$log.jacobian
                 )
    }
  }
  if(.temp.list$minimisation.function == "nlm")
    return(as.vector(neglik + penalty))
  else
    return(as.vector(neglik))
}
"proflik.lambda" <-
function(lambda)
{
  if (any(is.na(lambda)) | any(lambda==Inf) | any(is.nan(lambda)))
    neglik <- 1e+32
  else{
    if(.temp.list$minimisation.function == "nlm"){
      if (exists(".temp.lambda", w=1)) remove(".temp.lambda", pos=1, inherits = TRUE)
      lambda.minimiser <- lambda
      penalty <-  1000 * (.temp.lower.lambda - min(lambda, .temp.lower.lambda))
      lambda <- max(lambda, .temp.lower.lambda)
      penalty <- penalty + 1000 * (.temp.upper.lambda - max(lambda, .temp.upper.lambda))
      lambda <- min(lambda, .temp.upper.lambda)
      if (round(1000 * lambda.minimiser) <= round(1000 * .temp.lower.lambda))
        assign(".temp.lambda", lambda, pos=1)
      if (round(1000 * lambda.minimiser) >= round(1000 * .temp.upper.lambda))
        assign(".temp.lambda", lambda, pos=1)
    }
    z <- .temp.list$z
    n <- .temp.list$n
    if(lambda == 1) {
      .temp.list$log.jacobian <- 0
    }
    else {
      if(any(z < 0))
        stop("Transformation option not allowed when there are zeros or negative data"
             )
      if(any(z^(lambda - 1) <= 0))
        .temp.list$log.jacobian <- log(prod(z^(lambda - 1)))
      else .temp.list$log.jacobian <- sum(log(z^(lambda - 1)))
      if(lambda == 0)
        z <- log(z)
      else z <- ((z^lambda) - 1)/lambda
    }
    beta.size <- .temp.list$beta.size
    kappa <- .temp.list$kappa
    xmat <- .temp.list$xmat
    txmat <- .temp.list$txmat
    ixx <- solve(crossprod(xmat))
    tausqhat <- (z %*% (diag(n) - xmat %*% ixx %*% txmat) %*% z)/n
    if(.temp.list$method == "ML")
      neglik <- ((n/2) * log(2 * pi) +
                 (n/2) * log(tausqhat) +
                 (n/2) -
                 .temp.list$log.jacobian
                 )
    if(.temp.list$method == "RML") {
      eigentrem <- eigen(ixx, symmetric = TRUE, only.values = TRUE)
      neglik <- (((n - beta.size)/2) * log(2 * pi) +
                 ((n - beta.size)/2) * log(tausqhat) +
                 (n/2) -
                 0.5 * sum(log(eigentrem$values)) -
                 .temp.list$log.jacobian
                 )
    }
  }
  if(.temp.list$minimisation.function == "nlm")
    return(as.vector(neglik + penalty))
  else
    return(as.vector(neglik))
}

"proflik.nug" <-
  function (theta) 
{
  if (any(is.na(theta)) | any(theta==Inf) | any(is.nan(theta)))
    neglik <- 1e+32
  else{
    if(length(theta) == 3) include.lambda <- TRUE else include.lambda <- FALSE 
    if(.temp.list$minimisation.function == "nlm"){
      if (exists(".temp.phi", w=1)) remove(".temp.phi", pos=1, inherits = TRUE)
      if (exists(".temp.lambda", w=1)) remove(".temp.lambda", pos=1, inherits = TRUE)
      if (exists(".temp.nugget", w=1)) remove(".temp.nugget", pos=1, inherits = TRUE)
      theta.minimiser <- theta
      penalty <- 10000 * sum(.temp.lower - pmin(theta[1:2], .temp.lower))
      theta[1:2] <- pmax(theta[1:2], .temp.lower)
      if (theta.minimiser[1] <  .temp.lower[1])
        assign(".temp.nugget", theta[1], pos=1)
      if (theta.minimiser[2] < 1.001 * .temp.lower[2])
        assign(".temp.phi", theta[2], pos=1)
      if (include.lambda){
        lambda <- theta[3]
        penalty <- penalty + 1000 * (.temp.lower.lambda - min(lambda, .temp.lower.lambda))
        lambda <- max(lambda, .temp.lower.lambda)
        penalty <- penalty + 1000 * (.temp.upper.lambda - max(lambda, .temp.upper.lambda))
        lambda <- min(lambda, .temp.upper.lambda)
        if (round(1000 * theta.minimiser[3]) <= round(1000 * .temp.lower.lambda))
          assign(".temp.lambda", lambda, pos=1)
        if (round(1000 * theta.minimiser[3]) >= round(1000 * .temp.upper.lambda))
          assign(".temp.lambda", lambda, pos=1)
      }
    }
    else{
      if(include.lambda) lambda <- theta[3]
    }
    z <- .temp.list$z
    n <- .temp.list$n
    if(include.lambda){
      if(lambda == 1) {
        .temp.list$log.jacobian <<- 0
      }
      else {
        if(any(z < 0))
          stop("Transformation option not allowed when there are zeros or negative data"
               )
        if(any(z^(lambda - 1) <= 0))
          .temp.list$log.jacobian <<- log(prod(z^(lambda - 1)))
        else .temp.list$log.jacobian <<- sum(log(z^(lambda - 1)))
        if(lambda == 0)
          z <- log(z)
        else z <- ((z^lambda) - 1)/lambda
      }
    }
    beta.size <- .temp.list$beta.size
    kappa <- .temp.list$kappa
    tausq.rel <- theta[1]
    phi <- theta[2]
    covinf <- varcov.spatial(dists.lowertri = .temp.list$dists.lowertri,
                             cov.model = .temp.list$cov.model, kappa = kappa,
                             nugget = tausq.rel, cov.pars = c(1, phi),
                             det = TRUE, func.inv = "eigen",
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
    ssresmat <- as.vector(yiy - crossprod(beta.hat,xiy))
    if(.temp.list$method == "ML") {
      neglik <- ((n/2) * log(2 * pi) +
                 covinf$log.det.to.half +
                 (n/2) * log(ssresmat/n) +
                 (n/2) -
                 .temp.list$log.jacobian
                 )
    }
    if(.temp.list$method == "RML") {
      xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE, only.values = TRUE)
      neglik <- (((n - beta.size)/2) * log(2 * pi) +
                 covinf$log.det.to.half +
                 ((n - beta.size)/2) * log(ssresmat/(n-beta.size)) +
                 (n/2) +
                 choldet -
                 0.5 * sum(log(xx.eigen$values)) -
                 .temp.list$log.jacobian
                 )
    }
  }
  if(.temp.list$minimisation.function == "nlm")
    return(as.vector(neglik + penalty))
  else
    return(as.vector(neglik))
}

"proflik.phi" <-
  function (theta) 
{
  if (any(is.na(theta)) | any(theta==Inf) | any(is.nan(theta)))
    neglik <- 1e+32
  else{
    if(length(theta) == 2) include.lambda <- TRUE else include.lambda <- FALSE 
    if(.temp.list$minimisation.function == "nlm"){
      if (exists(".temp.phi", w=1)) remove(".temp.phi", pos=1, inherits = TRUE)
      if (exists(".temp.lambda", w=1)) remove(".temp.lambda", pos=1, inherits = TRUE)
      phi <- phi.minimiser <- theta[1]
      penalty <-  100000 * (.temp.lower.phi - min(phi, .temp.lower.phi))
      phi <- max(phi, .temp.lower.phi)
      if (phi.minimiser < 1.001 * .temp.lower.phi)
        assign(".temp.phi", phi, pos=1)
      if(include.lambda){
        lambda <- lambda.minimiser <- phi.lambda[2]
        penalty <-  penalty + 1000 * (.temp.lower.lambda - min(lambda, .temp.lower.lambda))
        lambda <- max(lambda, .temp.lower.lambda)
        penalty <- penalty + 1000 * (.temp.upper.lambda - max(lambda, .temp.upper.lambda))
        lambda <- min(lambda, .temp.upper.lambda)
        if (round(1000 * lambda.minimiser) <= round(1000 * .temp.lower.lambda))
          assign(".temp.lambda", lambda, pos=1)
        if (round(1000 * lambda.minimiser) >= round(1000 * .temp.upper.lambda))
          assign(".temp.lambda", lambda, pos=1)
      }
    }
    else{
      phi <- theta[1]
      if(include.lambda) lambda <- theta[2]
    }
    z <- .temp.list$z
    n <- .temp.list$n
    if(include.lambda){
      if(lambda == 1) {
        .temp.list$log.jacobian <<- 0
      }
      else {
        if(any(z <= 0))
          stop("Transformation option not allowed when there are zeros or negative data"
               )
        if(any(z^(lambda - 1) <= 0))
          .temp.list$log.jacobian <<- log(prod(z^(lambda - 1)))
        else .temp.list$log.jacobian <<- sum(log(z^(lambda - 1)))
        if(lambda == 0)
          z <- log(z)
        else z <- ((z^lambda) - 1)/lambda
      }
    }
    beta.size <- .temp.list$beta.size
    kappa <- .temp.list$kappa
    covinf <- varcov.spatial(dists.lowertri = .temp.list$
                             dists.lowertri,
                             cov.model = .temp.list$cov.model,
                             kappa = kappa, nugget = 0,
                             cov.pars = c(1, phi),
                             det = TRUE, func.inv = "eigen",
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
    ssresmat <- as.vector(yiy - crossprod(beta.hat,xiy))
    if(.temp.list$method == "ML") {
      neglik <- ((n/2) * log(2 * pi) + covinf$log.det.to.half +
                 (n/2) * log(ssresmat/n) + (n/2)) - .temp.list$
      log.jacobian
    }
    if(.temp.list$method == "RML") {
      xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE, only.values = TRUE)
      neglik <- (((n - beta.size)/2) * log(2 * pi) +
                 covinf$log.det.to.half +
                 ((n - beta.size)/2) * log(ssresmat/(n-beta.size)) +
                 (n/2) +
                 choldet -
                 0.5 * sum(log(xx.eigen$values)) -
                 .temp.list$log.jacobian
                 )
    }
  }
  if(.temp.list$minimisation.function == "nlm")
    return(as.vector(neglik + penalty))
  else
    return(as.vector(neglik))
}








