"olsfit" <-
  function (vario, ini.cov.pars, cov.model = c("exponential", "matern", "gaussian",
                          "spherical", "circular", "cubic", "wave",
                          "powered.exponential", "cauchy", "gneiting",
                          "gneiting.matern", "pure.nugget"),
            fix.nugget = FALSE, nugget = 0, 
            kappa = NULL, simul.number = NULL,  max.dist = "all",
            minimisation.function = c("optim", "nlm", "nls"),
            lower = 0, messages.screen=TRUE, ...) 
{
  call.fc <- match.call()
  ini <- ini.cov.pars
  ini.cov.pars <- NULL
  minimisation.function <- match.arg(minimisation.function)
  cov.model <- match.arg(cov.model)
  if (is.matrix(vario$v) & is.null(simul.number)) 
    stop("object in vario$v is a matrix. This function works for only 1 empirical variogram at once")
  if (!is.null(simul.number)) 
    vario$v <- vario$v[, simul.number]
##  if(is.matrix(ini)) {
##    inilength <- dim(ini)[2]
##    if(fix.nugget == FALSE & inilength != 3)
##      stop("wrong number of columns for ini (must be 3)")
##    if(fix.nugget == TRUE & inilength != 2)
##      stop("wrong number of columns for ini (must be 2)")
##  }
##  else {
##    inilength <- length(ini)
##    if (fix.nugget == FALSE & inilength != 3) 
##      stop("wrong length for ini (must be 3)")
##    if (fix.nugget == TRUE & inilength != 2) 
##      stop("wrong length for ini (must be 2)")
## }
  ##
  ##  Checking and/or preparing initial values
  ##
  if(is.matrix(ini)) {
    inilength <- dim(ini)[2]
    if(inilength != 2)
      stop("wrong number of columns for ini (must be always 2)")
    if(fix.nugget == FALSE)
      ini <- as.matrix(nugget, expand.grid(as.numeric(names(table(ini[,1]))),
                                           as.numeric(names(table(ini[,1])))))    
  }
  else {
    inilength <- length(ini)
    if(inilength != 2)
      stop("wrong length for ini (must be 2)")
    if(fix.nugget == FALSE){
      if(length(nugget) == 1)
        ini <- c(ini, nugget)
      else
      ini <- as.matrix(expand.grid(nugget, ini[1], ini[2]))    
    }
  }
  ##
  ## Setting maximum distance
  ##
  if (max.dist == "all") 
    data <- data.frame(u = vario$u, v = vario$v)
  else
    data <- data.frame(u = vario$u[vario$u <= max.dist],
                       v = vario$v[vario$u <= max.dist])
  ##
  ## Preparing lists for the minimiser
  ##
  .global.list <- list(u = data$u, v = data$v, fix.nugget = fix.nugget,
                       nugget = nugget, kappa = kappa,
                       cov.model = cov.model, m.f = minimisation.function)
  if(minimisation.function == "nlm")
    assign(".global.lower",lower,where=1)
  ##
  ## Choosing the best initial value
  ##
  if(is.matrix(ini) | is.data.frame(ini)) {
    ini <- as.matrix(ini)
    if(messages.screen) cat("olsfit: searching for the best initial value\n")
    ini.search <- ini
    dimnames(ini.search) <- list(NULL, NULL)
    ols.ini <- round(100000000. * apply(ini.search, 1, loss.olsvario, .global.list = .global.list))
    ini <- as.vector(ini.search[ols.ini == min(ols.ini),])
    if(messages.screen){
      cat("olsfit: best initial value:\n")
      if(fix.nugget == TRUE){
        names(ini) <- c("sill", "range")
        print(ini)
      }
      else{
        names(ini) <- c("nugget", "sill", "range")
        print(ini)
      }
    }
    remove(".global.lower", where=1)
  }
  ##
  ## minimisation using "nls"
  ##
  if (minimisation.function == "nls") {
    require(nls)
    data$kappa <- kappa
    data$cov.model <- cov.model
    if (fix.nugget == FALSE) {
      result <- nls(v ~ ((sigmasq + tausq) -
                         cov.spatial(u, cov.model = cov.model, kappa = kappa,
                                     cov.pars = c(sigmasq, phi))),
                    start = list(tausq = ini[1], sigmasq = ini[2], phi = ini[3]),
                    data = data)
      nugget <- result$parameter[1]
      cov.pars <- as.vector(result$parameter[2:3])
      value <- data$v - ((nugget + cov.pars[1]) -
                          cov.spatial(data$u, cov.model = cov.model, kappa = kappa,
                                      cov.pars = cov.pars))
    }
    else {
      data$ftau <- nugget
      result <- nls(v ~ ((sigmasq + ftau) -
                         cov.spatial(u, cov.model = cov.model, kappa = kappa,
                                     cov.pars = c(sigmasq, phi))),
                    start = list(sigmasq = ini[1], phi = ini[2]), data = data)
      cov.pars <- as.vector(result$parameter[1:2])
      value <- data$v - ((nugget + cov.pars[1]) -
                          cov.spatial(data$u, cov.model = cov.model, kappa = kappa,
                                      cov.pars = cov.pars))
    }
    value <- sum(value * value)
    message <- "nls does not provides convergence message"
  }
  ##
  ## minimisation using "optim" or "nlm"
  ##
  if (minimisation.function == "nlm" | minimisation.function == "optim") {
    if(minimisation.function == "nlm"){
      result <- nlm(loss.olsvario, ini, .global.list= .global.list, ...)
      result$par <- result$estimate
      result$value <- result$minimum
      result$convergence <- result$code
      if (exists(".temp.theta", where =1)){
        result$par <- .temp.theta
        remove(".temp.theta", pos=1)
      }
      remove(".global.lower", pos=1, inherits = TRUE)
    }
    else{
      result <- optim(ini, loss.olsvario, method = "L-BFGS-B",
                      lower = lower, .global.list = .global.list, ...)
    }
    value <- result$value
    message <- paste(minimisation.function, "convergence code:", result$convergence)
    if (fix.nugget == FALSE) {
      nugget <- result$par[1]
      cov.pars <- as.vector(result$par[2:3])
    }
    else {
      cov.pars <- as.vector(result$par[1:2])
    }
  }
  ##
  ## Preparing output
  ##
  dmax <- max(vario$u)
  estimation <- list(nugget = nugget, cov.pars = cov.pars, 
                     cov.model = cov.model, kappa = kappa, value = value, 
                     max.dist = dmax, minimisation.function = minimisation.function,
                     message = message)
  estimation$method <- "OLS"
  estimation$fix.nugget <- fix.nugget
  estimation$call <- call.fc
  class(estimation) <- "variomodel"
  if(messages.screen){
    cat("olsfit: estimated model parameters are:\n")
    cat(paste("covariance model:", cov.model))
    if(cov.model == "matern" | cov.model == "powered.exponential" | 
       cov.model == "cauchy" | cov.model == "gneiting.matern")
      cat(paste(" with kappa =", kappa))
    if(!is.null(kappa))
      if(cov.model == "matern" & kappa == 0.5)
        cat(" (exponential)")
    cat("\n")
    print(c(nugget=estimation$nugget, sill=estimation$cov.pars[1], range=estimation$cov.pars[2]))
  }
  return(estimation)
}


"wlsfit" <-
function (vario, ini.cov.pars, cov.model = c("exponential", "matern", "gaussian", "spherical",
                        "circular", "cubic", "wave", "powered.exponential", "cauchy",
                        "gneiting", "gneiting.matern", "pure.nugget"),
          fix.nugget = FALSE, nugget = 0,
          kappa = NULL, simul.number = NULL, max.dist = "all",
          minimisation.function = c("optim", "nlm"), lower = 0,
          weight = c("npairs", "cressie"), messages.screen=TRUE, ...) 
{
  call.fc <- match.call()
  ini <- ini.cov.pars
  ini.cov.pars <- NULL
  weight <- match.arg(weight)
  minimisation.function <- match.arg(minimisation.function)
  cov.model <- match.arg(cov.model)
  if (is.matrix(vario$v) & is.null(simul.number)) 
    stop("object in vario$v is a matrix. This function works for only 1 empirical variogram")
  if (!is.null(simul.number)) 
    vario$v <- vario$v[, simul.number]
##  if(is.matrix(ini)) {
##    inilength <- dim(ini)[2]
##    if(fix.nugget == FALSE & inilength != 3)
##      stop("wrong number of columns for ini (must be 3)")
##    if(fix.nugget == TRUE & inilength != 2)
##      stop("wrong number of columns for ini (must be 2)")
##  }
##  else {
##    inilength <- length(ini)
##    if (fix.nugget == FALSE & inilength != 3) 
##      stop("wrong length for ini (must be 3)")
##    if (fix.nugget == TRUE & inilength != 2) 
##      stop("wrong length for ini (must be 2)")
##  }
  ##
  ##  Checking and/or preparing initial values
  ##
  if(is.matrix(ini)) {
    inilength <- dim(ini)[2]
    if(inilength != 2)
      stop("wrong number of columns for ini (must be always 2)")
    if(fix.nugget == FALSE)
      ini <- as.matrix(expand.grid(nugget, as.numeric(names(table(ini[,1]))),
                                   as.numeric(names(table(ini[,1])))))    
  }
  else {
    inilength <- length(ini)
    if(inilength != 2)
      stop("wrong length for ini (must be 2)")
    if(fix.nugget == FALSE){
      if(length(nugget) == 1)
        ini <- c(ini, nugget)
      else
      ini <- as.matrix(expand.grid(nugget, ini[1], ini[2]))    
    }
  }
  ##
  ## Setting maximum distance
  ##
  if (max.dist == "all") 
    data <- list(u = vario$u, v = vario$v, n = vario$n)
  else data <- list(u = vario$u[vario$u <= max.dist],
                    v = vario$v[vario$u <= max.dist],
                    n = vario$n[vario$u <= max.dist])
  holdit <- data$v
  holdit[is.na(holdit)] <- 0
  data$v <- holdit
  ##
  ## Preparing lists for the minimiser
  ##
  .global.list <- list(u = data$u, v = data$v, n = data$n,
                       fix.nugget = fix.nugget,
                       nugget = nugget, 
                       kappa = kappa, cov.model = cov.model,
                       weight = weight, m.f = minimisation.function)
  if(minimisation.function == "nlm")
    assign(".global.lower",lower,where=1)
  ##
  ## Searchin for best initial value 
  ##
  if(is.matrix(ini) | is.data.frame(ini)) {
    ini <- as.matrix(ini)
    if(messages.screen) cat("wlsfit: searching for the best initial value\n")
    ini.search <- ini
    dimnames(ini.search) <- list(NULL, NULL)
    wls.ini <- round(100000000. * apply(ini.search, 1, loss.wlsvario,
                                        .global.list = .global.list))
    ini <- as.vector(ini.search[wls.ini == min(wls.ini),])
    if(messages.screen){
      cat("wlsfit: best initial value:\n")
      if(fix.nugget == TRUE){
        names(ini) <- c("sill", "range")
        print(ini)
      }
      else{
        names(ini) <- c("nugget", "sill", "range")
        print(ini)
      }
    }
    remove(".global.lower", where=1)
  }
  ##
  ## Minimisation isung "optim" or "nlm" 
  ##
  if(minimisation.function == "nlm"){
    result <- nlm(loss.wlsvario, ini, .global.list = .global.list, ...)
    result$par <- result$estimate
    result$value <- result$minimum
    result$convergence <- result$code
    if (exists(".temp.theta", w=1)){
      result$estimate <- .temp.theta
      remove(".temp.theta", pos=1, inherits = TRUE)
    }
    remove(".global.lower", pos=1, inherits = TRUE)
  }
  else{
    result <- optim(ini, loss.wlsvario, method = "L-BFGS-B",
                    lower = lower, .global.list = .global.list, ...)
  }
  ##
  ## Preparing output 
  ##
  value <- result$value
  message <- paste(minimisation.function, "convergence code:", result$convergence)
  if (fix.nugget == FALSE) {
    nugget <- result$par[1]
    cov.pars <- as.vector(result$par[2:3])
  }
  else {
    cov.pars <- as.vector(result$par[1:2])
  }
  estimation <- list(nugget = nugget, cov.pars = cov.pars, 
                     cov.model = cov.model, kappa = kappa, value = value, 
                     max.dist = max(data$u),
                     minimisation.function = minimisation.function,
                     message = message)
  estimation$method <- "WLS"
  estimation$fix.nugget <- fix.nugget
  estimation$call <- call.fc
  class(estimation) <- "variomodel"
  if(messages.screen){
    cat("wlsfit: estimated model parameters are:\n")
    cat(paste("covariance model:", cov.model))
    if(cov.model == "matern" | cov.model == "powered.exponential" | 
       cov.model == "cauchy" | cov.model == "gneiting.matern")
      cat(paste(" with kappa =", kappa))
    if(!is.null(kappa))
      if(cov.model == "matern" & kappa == 0.5)
        cat(" (exponential)")
    cat("\n")
    print(c(nugget=estimation$nugget, sill=estimation$cov.pars[1], range=estimation$cov.pars[2]))
  }
  return(estimation)
}


"loss.olsvario" <-
  function (theta, .global.list) 
{
  if(.global.list$m.f == "nlm"){
    if (exists(".temp.theta", w=1))
      remove(".temp.theta", pos=1, inherits = TRUE)
    theta.minimiser <- theta
    penalty <- 10000 * sum(.global.lower - pmin(theta, .global.lower))
    theta <- pmax(theta, .global.lower)
    if (any(theta.minimiser < .global.lower))
      assign(".temp.theta", theta, pos=1)
  }
  else penalty <- 0
  u <- .global.list$u
  v <- .global.list$v
  kappa <- .global.list$kappa
  if (.global.list$fix.nugget == FALSE) {
    tausq <- theta[1]
    sigmasq <- theta[2]
    phi <- theta[3]
    tau1 <- tausq
  }
  else {
    tausq <- .global.list$nugget
    sigmasq <- theta[1]
    phi <- theta[2]
    tau1 <- tausq + 1
  }
  sill.total <- sigmasq + tausq
  if (.global.list$fix.nugget == FALSE) {
    gamma <- sill.total - cov.spatial(u, cov.model = .global.list$cov.model, 
                                      kappa = kappa, cov.pars = c(theta[2], theta[3]))
  }
  else {
    gamma <- sill.total - cov.spatial(u, cov.model = .global.list$cov.model, 
                                      kappa = kappa, cov.pars = c(theta[1], theta[2]))
  }
  loss <- sum((v - gamma)^2)
  return(loss + penalty)
}

"loss.wlsvario" <-
  function (theta, .global.list) 
{
  if(.global.list$m.f == "nlm"){
    if (exists(".temp.theta", w=1))
      remove(".temp.theta", pos=1, inherits = TRUE)
    theta.minimiser <- theta
    penalty <- 10000 * sum(.global.lower - pmin(theta, .global.lower))
    theta <- pmax(theta, .global.lower)
    if (any(theta.minimiser < .global.lower))
      assign(".temp.theta", theta, pos=1)
  }
  else penalty <- 0
  u <- .global.list$u
  v <- .global.list$v
  n <- .global.list$n
  kappa <- .global.list$kappa
  if (.global.list$fix.nugget == FALSE) {
    tausq <- theta[1]
    sigmasq <- theta[2]
    phi <- theta[3]
    tau1 <- tausq
  }
  else {
    tausq <- .global.list$nugget
    sigmasq <- theta[1]
    phi <- theta[2]
    tau1 <- tausq + 1
  }
  sill.total <- sigmasq + tausq
  if (.global.list$fix.nugget == FALSE) {
    gamma <- sill.total - cov.spatial(u, cov.model = .global.list$cov.model, 
                                      kappa = kappa, cov.pars = c(theta[2], theta[3]))
  }
  else {
    gamma <- sill.total - cov.spatial(u, cov.model = .global.list$cov.model, 
                                      kappa = kappa, cov.pars = c(theta[1], theta[2]))
  }
  if (.global.list$weight == "npairs") 
    loss <- sum(n * (v - gamma)^2)
  if (.global.list$weight == "cressie") 
    loss <- sum((n/(gamma^2)) * (v - gamma)^2)
  return(loss + penalty)
}

