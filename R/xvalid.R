"xvalid" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            model, reestimate = FALSE, variog.obj = NULL,
            output.reestimate = FALSE, locations.xvalid = "all",
            messages.screen = TRUE, ...) 
{
  n <- nrow(coords)
  data <- as.vector(data)
  if (length(data) != n) 
    stop("coords and data have incompatible dimentions")
  xmat <- trend.spatial(trend = model$trend, coords = coords)
  ##
  ## Locations to be used in the cross-validation
  ##
  if(locations.xvalid == "all") locations.xvalid <- 1:n
  else{
    if(is.matrix(locations.xvalid) | is.data.frame(locations.xvalid))
      locations.xvalid <- is.vector(locations.xvalid)
    if(!is.vector(locations.xvalid) | !is.numeric(locations.xvalid))
      stop("argument locations.xvalid must be a numeric vector with numbers indicating the locations to be cross-validated")
    if(any(locations.xvalid) > n | length(locations.xvalid) > n)
      stop("incorrect value to the argument locations.xvalid.\nThis must be a numeric vector with numbers indicating the locations to be cross-validated")
  }
  n.pt.xv <- length(locations.xvalid)
  if(messages.screen){
    cat(paste("xvalid: number of data locations                 =", n))
    cat("\n")
    cat(paste("xvalid: number of locations for cross-validation =", n.pt.xv))
    cat("\n")
    cat("xvalid: performing cross-validation at location ... ")
    }
  ##
  ## Defining a function to predict at one point
  ## 
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
          stop("an object with the estimated variogram must be provided in the argument variog.obj when argument reestimate = TRUE")
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
        if(is.null(model$call$lower)) model$call$lower <- 0
        if(model$method == "OLS")
          CVmod <- olsfit(vario = CVvar, ini=model$cov.pars, cov.model = model$cov.model,
                          fix.nugget = model$fix.nugget, nugget = model$nugget,
                          kappa = model$kappa, max.dist = model$max.dist,
                          minimisation.function = model$minimisation.function,
                          lower = model$call$lower, messages.screen = FALSE, ...)
        if(model$method == "WLS"){
          if(is.null(model$call$weight)) model$call$weight <- "npairs"
          CVmod <- wlsfit(vario = CVvar, ini=model$cov.pars, cov.model = model$cov.model,
                          fix.nugget = model$fix.nugget, nugget = model$nugget,
                          kappa = model$kappa, max.dist = model$max.dist,
                          minimisation.function = model$minimisation.function,
                          lower = model$call$lower, weight = model$call$weight,
                          messages.screen = FALSE, ...)
        }
        if(output.reestimate){
          CVpars <- CVmod$cov.pars
          if(CVmod$fix.nugget == FALSE) CVpars <- c(CVpars, CVmod$nugget)
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
      fix.pars <- c(CVmod$fix.nugget, T,T,T,T)
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
  if(messages.screen) cat("\nxvalid: end of cross-validation\n")
  if(output.reestimate){
    pars.names <- c("sigmasq", "phi")
    if(model$method == "ML" | model$method == "REML" | model$method == "RML"){
      fix.pars <- (model$parameters.summary[c("tausq", "kappa", 
                                              "psiA", "psiR", "lambda"), 1] == "fixed")
      pars.names <- c(pars.names,(c("tausq", "kappa", "psiA", "psiR", "lambda"))[fix.pars == FALSE])
    }
    if(model$method == "OLS" | model$method == "WLS")
      if(model$fix.nugget == FALSE) pars.names <- c(pars.names, "tausq")
    names(res) <- c(c("data", "pred", "krige.var"), pars.names)
  }
  else names(res) <- c("data", "predicted", "krige.var")
  res$error <- res$data - res$pred
  res$std.error <- res$err/sqrt(res$krige.var)
  res$prob <- pnorm(res$data, mean = res$pred, sd = sqrt(res$krige.var))
  if(output.reestimate){
    np <- length(pars.names)
    res <- res[,c((1:3), ((3+np+1):(6+np)),(4:(3+np)))] 
  }
  class(res) <- "xvalid"
  attr(res,"row.names") <- NULL
  return(res)
}

"plot.xvalid" <-
  function (obj, geodata, coords=geodata$coords, error.plots = TRUE, std.error.plots = TRUE, ask = TRUE)
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
