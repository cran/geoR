"xvalid" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            model, reestimate = FALSE, variog.obj = NULL,
            output.reestimate = FALSE, locations.xvalid = "all",
            data.xvalid = NULL, messages.screen = TRUE, ...) 
{
  n <- nrow(coords)
  data <- as.vector(data)
  if (length(data) != n) 
    stop("coords and data have incompatible dimentions")
  xmat <- trend.spatial(trend = model$trend, coords = coords)
  ##
  ## Locations to be used in the cross-validation
  ##
  if(locations.xvalid == "all" | is.vector(locations.xvalid)){
    if(locations.xvalid == "all")
      locations.xvalid <- 1:n
    else
      if(any(locations.xvalid > n) | !is.numeric(locations.xvalid))
        stop("\nxvalid: vector indicating locations to be validated is not a numeric vector and/or has element(s) with value greater than the number of data loccations")
    crossvalid <- TRUE
  }
  else{
    if(is.matrix(locations.xvalid) | is.data.frame(locations.xvalid))
      if(dim(locations.xvalid)[2] <= 2){
        if(dim(locations.xvalid)[2] == 1){
          locations.xvalid <- is.vector(locations.xvalid)
          crossvalid <- TRUE
        }
        else{
          if(messages.screen)
            cat("xvalid: cross-validation to be performed on locations provided by the user\n")
          if(is.null(data.xvalid))
            stop("the argument \"data.xvalid\" must be provided in order to perform validation on a set of locations different from the original data")
          crossvalid <- FALSE
        }
      }
      else
        if(!is.vector(locations.xvalid) | !is.numeric(locations.xvalid))
          stop("\nargument locations.xvalid must be either:\n a numeric vector with numbers indicating the locations to be cross-validated\n a matrix with coordinates for the locations to be cross-validated.")
    if(any(locations.xvalid) > n | length(locations.xvalid) > n)
      stop("incorrect value to the argument locations.xvalid.\nThis must be a numeric vector with numbers indicating the locations to be cross-validated")
  }
  n.pt.xv <- length(locations.xvalid)
  if(crossvalid == FALSE) n.pt.xv <- dim(locations.xvalid)[[1]]
  if(messages.screen){
    cat(paste("xvalid: number of data locations       =", n))
    cat("\n")
    cat(paste("xvalid: number of validation locations =", n.pt.xv))
    cat("\n")
    if(crossvalid) cat("xvalid: performing cross-validation at location ... ")
    else  cat("xvalid: performing validation at the locations provided")
    }
  ##
  ## Defining a function to predict at one point
  ##
  if(crossvalid){
    cv.f <- function(ndata, ...) {
      if(messages.screen) cat(paste(ndata, ", ", sep=""))
      ## excluding data point
      coords.out <- coords[ndata, , drop = FALSE]
      data.out <- data[ndata]
      xmat.out <- xmat[ndata, , drop = FALSE]
      cv.coords <- coords[-ndata, ]
      cv.data <- as.vector(data)[-ndata]
      cv.xmat <- xmat[-ndata, , drop = FALSE]
      ## re-estimating the model
      if (reestimate) {
        if(model$method == "ML" | model$method == "REML" | model$method == "RML"){
          fix.pars <- (model$parameters.summary[c("tausq", "kappa", "psiA",
                                                  "psiR", "lambda"), 1] == "fixed")
          val.pars <- model$parameters.summary[c("tausq", "kappa", 
                                                 "psiA", "psiR", "lambda"), 2]
          names(fix.pars) <- c("tausq", "kappa", "psiA", "psiR", 
                               "lambda")
          names(val.pars) <- c("tausq", "kappa", "psiA", "psiR", 
                               "lambda")
          CVmod <- likfit(coords = cv.coords, data = cv.data, 
                          ini = model$cov.pars, fix.nugget = fix.pars["tausq"], 
                          nugget = val.pars["tausq"], fix.kappa = fix.pars["kappa"], 
                          kappa = val.pars["kappa"], fix.psiR = fix.pars["psiR"], 
                          psiR = val.pars["psiR"], fix.psiA = fix.pars["psiA"], 
                          psiA = val.pars["psiA"], fix.lambda = fix.pars["lambda"], 
                          lambda = val.pars["lambda"], cov.model = model$cov.model, 
                          trend = ~cv.xmat + 0, method = model$method, 
                          messages.screen = F, ...)
          if(output.reestimate){
            CVpars <- (CVmod$parameters.summary[c("tausq", "kappa", "psiA", "psiR", "lambda"), 2])
            CVpars <- c(CVmod$cov.pars, CVpars[fix.pars == FALSE])
          } 
        }
        if(model$method == "OLS" | model$method == "WLS"){
          if(is.null(variog.obj))
            stop("xvalid: when argument reestimate = TRUE an object with the fitted variogram model must be provided in the argument variog.obj ")
          CVvar <- variog(coords = cv.coords, data = cv.data, uvec = variog.obj$uvec,
                          trend = variog.obj$trend, lambda = variog.obj$lambda,
                          option = variog.obj$output.type,
                          estimator.type = variog.obj$estimator.type,
                          nugget.tolerance = variog.obj$nugget.tolerance,
                          max.dist = max(variog.obj$u), pairs.min = 2,
                          bin.cloud = FALSE, direction = variog.obj$direction,
                          tolerance = variog.obj$tolerance,
                          unit.angle = "radians",
                          messages.screen = FALSE, ...)
          CVmod <- variofit(vario = CVvar, ini=model$cov.pars, cov.model = model$cov.model,
                            fix.nugget = model$fix.nugget, nugget = model$nugget,
                            fix.kappa = model$fix.kappa, kappa = model$kappa, max.dist = model$max.dist,
                            minimisation.function = model$minimisation.function,
                            weights = model$weights, messages.screen = FALSE, ...)
          if(output.reestimate){
            CVpars <- CVmod$cov.pars
            if(CVmod$fix.nugget == FALSE) CVpars <- c(CVpars, CVmod$nugget)
            if(CVmod$fix.kappa == FALSE) CVpars <- c(CVpars, CVmod$kappa)
          }
        }
      }
      else CVmod <- model
      if(model$method == "ML" | model$method == "REML" | model$method == "RML"){
        fix.pars <- (CVmod$parameters.summary[c("tausq", "kappa", 
                                                "psiA", "psiR", "lambda"), 1] == "fixed")
        val.pars <- CVmod$parameters.summary[c("tausq", "kappa", 
                                               "psiA", "psiR", "lambda"), 2]
      }
      if(model$method == "OLS" | model$method == "WLS"){
        fix.pars <- c(CVmod$fix.nugget, CVmod$fix.kappa,T,T,T)
        if(is.null(CVmod$kappa)) CVmod$kappa <- 0.5
        val.pars <- c(CVmod$nugget, CVmod$kappa, 0, 1, CVmod$lambda)
      }
      names(fix.pars) <- c("tausq", "kappa", "psiA", "psiR", "lambda")
      names(val.pars) <- c("tausq", "kappa", "psiA", "psiR", "lambda")
      kr <- krige.conv(coords = cv.coords, data = cv.data, loc = coords.out,
                       krige = krige.control(trend.d = ~cv.xmat + 0,
                         trend.l = ~xmat.out + 0, cov.model = CVmod$cov.model, 
                         cov.pars = CVmod$cov.pars, nugget = CVmod$nugget, 
                         kappa = val.pars["kappa"], lambda = val.pars["lambda"], 
                         aniso.pars = val.pars[c("psiA", "psiR")]), mess = FALSE)
      res <- c(data.out, kr$pred, kr$krige.var)
      if(output.reestimate) res <- c(res, CVpars)
      ##, err = (data.out - kr$pred), e.rel = (data.out - kr$pred)/sqrt(kr$krige.var), 
      ##pval = pnorm(data.out, mean = kr$pred, sd = sqrt(kr$krige.var)))
      return(res)
    }
    res <- as.data.frame(t(apply(matrix(locations.xvalid), 1, cv.f)))
  }
  else{
    xmat.val.loc <- trend.spatial(trend = model$trend, coords = locations.xvalid)
    if(model$method == "ML" | model$method == "REML" | model$method == "RML"){
      fix.pars <- (model$parameters.summary[c("tausq", "kappa", 
                                              "psiA", "psiR", "lambda"), 1] == "fixed")
      val.pars <- model$parameters.summary[c("tausq", "kappa", 
                                             "psiA", "psiR", "lambda"), 2]
    }
    if(model$method == "OLS" | model$method == "WLS"){
      fix.pars <- c(model$fix.nugget, model$fix.kappa,T,T,T)
      if(is.null(model$kappa)) model$kappa <- 0.5
      val.pars <- c(model$nugget, model$kappa, 0, 1, model$lambda)
    }
    names(fix.pars) <- c("tausq", "kappa", "psiA", "psiR", "lambda")
    names(val.pars) <- c("tausq", "kappa", "psiA", "psiR", "lambda")
    res <- krige.conv(coords = coords, data = data, loc = locations.xvalid,
                     krige = krige.control(trend.d = ~xmat + 0,
                       trend.l = ~xmat.val.loc + 0, cov.model = model$cov.model, 
                       cov.pars = model$cov.pars, nugget = model$nugget, 
                       kappa = val.pars["kappa"], lambda = val.pars["lambda"], 
                       aniso.pars = val.pars[c("psiA", "psiR")]), mess = FALSE)[1:2]
    res <- data.frame(data.xvalid, res$pred, res$krige.var)
  } 
  if(messages.screen) cat("\nxvalid: end of cross-validation\n")
  if(output.reestimate){
    pars.names <- c("sigmasq", "phi")
    if(model$method == "ML" | model$method == "REML" | model$method == "RML"){
      fix.pars <- (model$parameters.summary[c("tausq", "kappa", 
                                              "psiA", "psiR", "lambda"), 1] == "fixed")
      pars.names <- c(pars.names,(c("tausq", "kappa", "psiA", "psiR", "lambda"))[fix.pars == FALSE])
    }
    if(model$method == "OLS" | model$method == "WLS"){
      if(model$fix.nugget == FALSE) pars.names <- c(pars.names, "tausq")
      if(model$fix.kappa == FALSE) pars.names <- c(pars.names, "kappa")
    }
      names(res) <- c(c("data", "predicted", "krige.var"), pars.names)
  }
  else names(res) <- c("data", "predicted", "krige.var")
  res$error <- res$data - res$pred
  res$std.error <- res$err/sqrt(res$krige.var)
  res$prob <- pnorm(res$data, mean = res$pred, sd = sqrt(res$krige.var))
  if(output.reestimate){
    np <- length(pars.names)
    res <- res[,c((1:3), ((3+np+1):(6+np)),(4:(3+np)))] 
  }
  attr(res,"row.names") <- NULL
  attr(res, "class") <- "xvalid"
  return(res)
}

"plot.xvalid" <-
  function (obj, valid.obj, coords=valid.obj$coords, 
            error.plots = TRUE, std.error.plots = TRUE, ask = TRUE)
{
  ##
  ## Saving original par() parameters
  ##
  if (is.R()) 
    par.ori <- par(no.readonly = TRUE)
  else par.ori <- par()
  on.exit(par(par.ori))
  ##
  ## auxiliary computations for plots
  ##
  n <- length(obj$data)
  xylim <- range(c(obj$data, obj$pred))
  prelim <- range(obj$pred)
  datlim <- range(obj$data)
  errlim <- max(abs(range(obj$error)))
  errlim <- c(-errlim, errlim)
  err.std <- sqrt(var(obj$error))
  if(n > 90){
    seqerr <- seq(-3.5*err.std, 3.5*err.std, l=15)
    seqstd <- seq(-3.5, 3.5, l=15)
  }
  else{
    seqerr <- seq(-4*err.std, 4*err.std, l=9)
    seqstd <- seq(-4, 4, l=9)
  }
  stdlim <- max(c(3, abs(range(obj$std.error))))
  stdlim <- c(-stdlim, stdlim)
  # indicator for negative and positive errors
  error.cut <- cut(obj$error, breaks=c(errlim[1], 0, errlim[2]), include.l=TRUE, labels=FALSE)
  ##
  ## Data vs predicted
  ##
  par(pty = "s")
  plot(obj$data, obj$pred, type = "n", xlim = xylim, ylim = xylim,
       xlab = "data", ylab = "predicted")
  points(obj$data, obj$pred, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
  abline(0,1)
  ##
  ##
  ##
  par(ask = ask)
  ##
  if(!error.plots | !std.error.plots){
    ##
    ## P-P plot
    ##
    par(pty = "s")  
    plot(ppoints(n), obj$prob[order(obj$prob)], xlim=c(0,1), ylim=c(0,1), xlab="theoretical prob", ylab="observed prob")
    abline(0,1)
  }
  if(error.plots){
    ##
    ## Plotting errors
    ##
    ## sizes proportional to errors values
    err.abs <- abs(obj$error)
    coords.order <- coords[order(err.abs), ]
    err.order <- err.abs[order(err.abs)]
    cut.order <- error.cut[order(err.abs)]
    r.y <- range(err.order)
    err.size <- 0.7 + ((err.order - r.y[1]) * (2 - 0.7))/(r.y[2] - r.y[1])
    ## equal scale for plot
    coords.lims <- apply(coords, 2, range)
    coords.diff <- diff(coords.lims)
    if (coords.diff[1] != coords.diff[2]) {
      coords.diff.diff <- abs(diff(as.vector(coords.diff)))
      ind.min <- which(coords.diff == min(coords.diff))
      coords.lims[, ind.min] <- coords.lims[, ind.min] + c(-coords.diff.diff, 
                                                           coords.diff.diff)/2
    }
    par(pty = "s")
    ##
    plot(coords, xlab = "Coord X", ylab = "Coord Y", type = "n", 
         xlim = coords.lims[, 1], ylim = coords.lims[, 2])
    if (is.R()) {
      points(coords.order, pch = (c("x", "+"))[cut.order], col=(c("red", "blue"))[cut.order], cex = err.size)
    }
    ##
    ## errors histogram
    ##
    par(pty = "m")
    if(min(obj$error) < min(seqerr)) seqerr <- c(min(obj$error), seqerr)
    if(max(obj$error) > max(seqerr)) seqerr <- c(seqerr, max(obj$error))
    hist(obj$error, prob=T, main="", breaks=seqerr, xlab="data - predicted")
    ##
    ## errors vs predicted
    ##
    par(pty = "m")
    plot(obj$pred, obj$error, type = "n", xlim = prelim, ylim = errlim,
         xlab = "predicted", ylab = "data - predicted")
    points(obj$pred, obj$error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
    abline(h=0)
    ##
    ## errors vs data
    ##
    par(pty = "m")
    plot(obj$data, obj$error, type = "n", xlim = datlim, ylim = errlim,
         xlab = "data", ylab = "data - predicted")
    points(obj$data, obj$error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
    abline(h=0)
    ##
  }
  if(error.plots & std.error.plots){
    ##
    ## P-P plot
    ##
    par(pty = "s")  
    plot(ppoints(n), obj$prob[order(obj$prob)], xlim=c(0,1), ylim=c(0,1), xlab="theoretical prob", ylab="observed prob")
    abline(0,1)
  }
  if(std.error.plots){
    ##
    ## Plotting std errors
    ##
    ## sizes proportional to errors values
    err.abs <- abs(obj$std.error)
    coords.order <- coords[order(err.abs), ]
    err.order <- err.abs[order(err.abs)]
    cut.order <- error.cut[order(err.abs)]
    r.y <- range(err.order)
    err.size <- 0.7 + ((err.order - r.y[1]) * (2 - 0.7))/(r.y[2] - r.y[1])
    ## equal scale for plot
    coords.lims <- apply(coords, 2, range)
    coords.diff <- diff(coords.lims)
    if (coords.diff[1] != coords.diff[2]) {
      coords.diff.diff <- abs(diff(as.vector(coords.diff)))
      ind.min <- which(coords.diff == min(coords.diff))
      coords.lims[, ind.min] <- coords.lims[, ind.min] + c(-coords.diff.diff, 
                                                           coords.diff.diff)/2
    }
    par(pty = "s")
    ##
    plot(coords, xlab = "Coord X", ylab = "Coord Y", type = "n", 
         xlim = coords.lims[, 1], ylim = coords.lims[, 2])
    if (is.R()) {
      points(coords.order, pch = (c("x", "+"))[cut.order], col=(c("red", "blue"))[cut.order], cex = err.size)
    }
    ##
    ## std. errors histogram
    ##
    par(pty = "m")
    if(min(obj$std.error) < min(seqstd)) seqstd <- c(min(obj$std.error), seqstd)
    if(max(obj$std.error) > max(seqstd)) seqstd <- c(seqstd, max(obj$std.error))
    hist(obj$std.error, prob=T, main="", breaks = seqstd, xlab="std error")
    ##
    ## std. errors vs predicted
    ##
    par(pty = "m")
    plot(obj$pred, obj$std.error, type = "n", xlim = prelim, ylim = stdlim,
         xlab = "predicted", ylab = "std error")
    points(obj$pred, obj$std.error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
    abline(h=0)
    ##
    ## std. errors vs data
    ##
    par(pty = "m")
    plot(obj$data, obj$std.error, type = "n", xlim = datlim, ylim = stdlim,
         xlab = "data", ylab = "std error")
    points(obj$data, obj$std.error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
    abline(h=0)
    ##
  }
  ##
  return(invisible())
}
