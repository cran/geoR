"xvalid" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            model, reestimate = FALSE, variog.obj = NULL,
            output.reestimate = FALSE, locations.xvalid = "all",
            data.xvalid = NULL, messages.screen = TRUE, ...) 
{
  call.fc <- match.call()
  if(missing(geodata))
    geodata <- list(coords = coords, data = data)
  n <- nrow(coords)
  data <- as.vector(data)
  if (length(data) != n) 
    stop("coords and data have incompatible sizes")
  xmat <- trend.spatial(trend = model$trend, geodata = geodata)
  if (nrow(xmat) != n) 
    stop("coords and trend have incompatible sizes")
  ##
  ## Locations to be used in the cross-validation
  ##
  if(all(locations.xvalid == "all") | is.vector(locations.xvalid)){
    autocross <- TRUE
    if(locations.xvalid == "all")
      locations.xvalid <- 1:n
    else
      if(any(locations.xvalid > n) | !is.numeric(locations.xvalid))
        stop("\nxvalid: vector indicating locations to be validated is not a numeric vector and/or has element(s) with value greater than the number of data loccations")
    crossvalid <- TRUE
  }
  else{
    autocross <- FALSE
    if(is.matrix(locations.xvalid) | is.data.frame(locations.xvalid))
      if(dim(locations.xvalid)[2] <= 2){
        if(dim(locations.xvalid)[2] == 1){
          locations.xvalid <- is.vector(locations.xvalid)
          crossvalid <- TRUE
          if(any(locations.xvalid) > n | length(locations.xvalid) > n)
            stop("incorrect value to the argument locations.xvalid.\nThis must be a numeric vector with numbers indicating the locations to be cross-validated")
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
        else
          if(any(locations.xvalid) > n | length(locations.xvalid) > n)
            stop("incorrect value to the argument locations.xvalid.\nThis must be a numeric vector with numbers indicating the locations to be cross-validated")
  }
  if(crossvalid == FALSE) n.pt.xv <- dim(locations.xvalid)[[1]]
  else n.pt.xv <- length(locations.xvalid)
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
                          messages.screen = FALSE, ...)
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
        fix.pars <- c(CVmod$fix.nugget, CVmod$fix.kappa,TRUE,TRUE,TRUE)
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
                         aniso.pars = val.pars[c("psiA", "psiR")]),
                       output = output.control(mess = FALSE))
      res <- c(data.out, kr$pred, kr$krige.var)
      if(output.reestimate) res <- c(res, CVpars)
      ##, err = (data.out - kr$pred), e.rel = (data.out - kr$pred)/sqrt(kr$krige.var), 
      ##pval = pnorm(data.out, mean = kr$pred, sd = sqrt(kr$krige.var)))
      return(res)
    }
    res <- as.data.frame(t(apply(matrix(locations.xvalid), 1, cv.f)))
  }
  else{
    xmat.val.loc <- trend.spatial(trend = model$trend, geodata = list(coords=locations.xvalid))
    if(model$method == "ML" | model$method == "REML" | model$method == "RML"){
      fix.pars <- (model$parameters.summary[c("tausq", "kappa", 
                                              "psiA", "psiR", "lambda"), 1] == "fixed")
      val.pars <- model$parameters.summary[c("tausq", "kappa", 
                                             "psiA", "psiR", "lambda"), 2]
    }
    if(model$method == "OLS" | model$method == "WLS"){
      fix.pars <- c(model$fix.nugget, model$fix.kappa,TRUE,TRUE,TRUE)
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
                       aniso.pars = val.pars[c("psiA", "psiR")]),
                      output = output.control(mess = FALSE))[1:2]
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
  if(autocross){
    if(!is.null(call.fc$geodata))
      attr(res,"geodata.xvalid") <- call.fc$geodata
    else
      attr(res,"locations.xvalid") <- call.fc$locations.xvalid
  }
  else
    if(!is.null(locations.xvalid))
      attr(res,"locations.xvalid") <- call.fc$locations.xvalid
  attr(res, "class") <- "xvalid"
  return(res)
}

"plot.xvalid" <-
  function (x, coords, borders = NULL, ask = TRUE,
            error = TRUE, std.error = TRUE,
            data.predicted = TRUE,
            pp = TRUE, map = TRUE, histogram = TRUE,
            error.predicted = TRUE, error.data = TRUE, ...)
{
  ##
  ## Saving original par() parameters
  ##
  if (is.R()) 
    par.ori <- par(no.readonly = TRUE)
  else par.ori <- par()
  on.exit(par(par.ori))
  ##
  ## checking input
  ##
  if(!is.null(borders)){
    if(!is.matrix(borders) & !is.data.frame(borders))
      stop("argument borders must be a two column matrix or a data frame with the coordinates of the borders")
    else
      if(ncol(borders) > 2)
        stop("argument borders must be a two column matrix or a data frame with the coordinates of the borders")
      else borders <- as.matrix(borders)
  }
  if(error | std.error){
    if(missing(coords)){
      if(!is.null(attr(x,"geodata.xvalid")))
        coords <- eval(attr(x,"geodata.xvalid"))$coords
      if(!is.null(attr(x,"locations.xvalid")))
        coords <- eval(attr(x,"locations.xvalid"))       
    }
    else{
      if(class(coords) == "geodata")
        coords <- coords$coords
    }
    if(!is.matrix(coords) & !is.data.frame(coords))
      stop("argument coords must be a two column matrix or a data frame with the data coordinates")
    else
      if(ncol(coords) > 2)
        stop("argument coords must be a two column matrix or a data frame with the data coordinates")
      else coords <- as.matrix(coords)
  }
  ##
  ## auxiliary computations for plots
  ##
  n <- length(x$data)
  xylim <- range(c(x$data, x$pred))
  prelim <- range(x$pred)
  datlim <- range(x$data)
  errlim <- max(abs(range(x$error)))
  errlim <- c(-errlim, errlim)
  err.std <- sqrt(var(x$error))
  if(n > 90){
    seqerr <- seq(-3.5*err.std, 3.5*err.std, l=15)
    seqstd <- seq(-3.5, 3.5, l=15)
  }
  else{
    seqerr <- seq(-4*err.std, 4*err.std, l=9)
    seqstd <- seq(-4, 4, l=9)
  }
  stdlim <- max(c(3, abs(range(x$std.error))))
  stdlim <- c(-stdlim, stdlim)
  ## indicator for negative and positive errors
  error.cut <- cut(x$error, breaks=c(errlim[1], 0, errlim[2]), include.l=TRUE, labels=FALSE)
  ##
  ## Data vs predicted
  ##
  if(data.predicted){
    par(pty = "s")
    plot(x$pred, x$data, type = "n", xlim = xylim, ylim = xylim,
         ylab = "data", xlab = "predicted")
    points(x$pred, x$data, ...)
    ##    points(x$pred, x$data, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
    abline(0,1)
  }
  ##
  par(ask = ask)
  ##
  if(!error | !std.error){
    ##
    ## P-P plot
    ##
    if(pp){
      par(pty = "s")  
      plot(ppoints(n), x$prob[order(x$prob)], xlim=c(0,1), ylim=c(0,1), xlab="theoretical prob", ylab="observed prob")
      abline(0,1)
    }
  }
  if(error){
    ##
    ## Plotting errors
    ##
    ## sizes proportional to errors values
    err.abs <- abs(x$error)
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
    ##
    if(map){
      par(pty = "s")
      ##
      plot(coords, xlab = "Coord X", ylab = "Coord Y",
           type = "n", 
           xlim = coords.lims[, 1], ylim = coords.lims[, 2])
      if (is.R()) {
        points(coords.order, pch = (c("x", "+"))[cut.order], col=(c("red", "blue"))[cut.order], cex = err.size)
      }
      else
        points(coords.order, pch = (c("x", "+"))[cut.order], col=(c(3, 4))[cut.order], cex = err.size)
      if(!is.null(borders))
        lines(borders)
    }
    ##
    ## errors histogram
    ##
    if(histogram){
      par(pty = "m")
      if(min(x$error) < min(seqerr)) seqerr <- c(min(x$error), seqerr)
      if(max(x$error) > max(seqerr)) seqerr <- c(seqerr, max(x$error))
      hist(x$error, prob=TRUE, main="", breaks=seqerr, xlab="data - predicted")
    }
    ##
    ## errors vs predicted
    ##
    if(error.predicted){
      par(pty = "m")
      plot(x$pred, x$error, type = "n", xlim = prelim, ylim = errlim,
           xlab = "predicted", ylab = "data - predicted")
    ##  points(x$pred, x$error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
      points(x$pred, x$error, ...)
      abline(h=0)
    }
    ##
    ## errors vs data
    ##
    if(error.data){
      par(pty = "m")
      plot(x$data, x$error, type = "n", xlim = datlim, ylim = errlim,
           xlab = "data", ylab = "data - predicted")
      ##      points(x$data, x$error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
      points(x$data, x$error, ...)
      abline(h=0)
      ##
    }
  }
  if(error & std.error){
    ##
    ## P-P plot
    ##
    if(pp){
      par(pty = "s")  
      plot(ppoints(n), x$prob[order(x$prob)], xlim=c(0,1), ylim=c(0,1), xlab="theoretical prob", ylab="observed prob")
      abline(0,1)
    }
  }
  if(std.error){
    ##
    ## Plotting std residuals
    ##
    ## sizes proportional to errors values
    err.abs <- abs(x$std.error)
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
    ##
    if(map){
      par(pty = "s")
      ##
      plot(coords, xlab = "Coord X", ylab = "Coord Y", type = "n", 
           xlim = coords.lims[, 1], ylim = coords.lims[, 2])
      if (is.R()) {
        points(coords.order, pch = (c("x", "+"))[cut.order], col=(c("red", "blue"))[cut.order], cex = err.size)
      }
      else
        points(coords.order, pch = (c("x", "+"))[cut.order], col=(c(3, 4))[cut.order], cex = err.size)
      if(!is.null(borders))
        lines(borders)
    }
    ##
    ## std. errors histogram
    ##
    if(histogram){
      par(pty = "m")
      if(min(x$std.error) < min(seqstd)) seqstd <- c(min(x$std.error), seqstd)
      if(max(x$std.error) > max(seqstd)) seqstd <- c(seqstd, max(x$std.error))
      hist(x$std.error, prob=TRUE, main="", breaks = seqstd, xlab="std residuals")
    }
    ##
    ## std. errors vs predicted
    ##
    if(error.predicted){
      par(pty = "m")
      plot(x$pred, x$std.error, type = "n", xlim = prelim, ylim = stdlim,
           xlab = "predicted", ylab = "std residuals")
      ##      points(x$pred, x$std.error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
      points(x$pred, x$std.error, ...)
      abline(h=0)
    }
    ##
    ## std. errors vs data
    ##
    if(error.data){
      par(pty = "m")
      plot(x$data, x$std.error, type = "n", xlim = datlim, ylim = stdlim,
           xlab = "data", ylab = "std residuals")
      ##      points(x$data, x$std.error, pch = (c("x", "+"))[error.cut], col=(c("red", "blue"))[error.cut])
      points(x$data, x$std.error,  ...)
      abline(h=0)
      ##
    }
  }
  ##
  return(invisible())
}
