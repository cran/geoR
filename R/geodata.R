##
## Basic data manipulation for the geoR package
## --------------------------------------------
##
## Functions for reading data and basic exploratory analysis. 
## These functions include
##    - creates objets of the class geodata
##    - methods for this class
##

"read.geodata" <-
  function(file, header = FALSE, coords.col= 1:2, data.col = 3,
           data.names = NULL, covar.col = NULL, covar.names = "header",
           realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar"),
           rep.data.action, rep.covar.action, ...)
{
  call.fc <- match.call()
  ##
  obj <- read.table(file = file, header = header, ...)
  if(covar.names == "header"){
    if(!is.null(covar.col)){
      col.names <- names(obj)
      covar.names <- col.names[covar.col]
    }
    else covar.names <- NULL
  }
  ##
  if(missing(rep.data.action)) rep.data.action <- "none"
  if(!is.function(rep.data.action))
    rep.data.action <- match.arg(rep.data.action, choices = c("none", "first")) 
  if(missing(rep.covar.action)) rep.covar.action <- rep.data.action
  if(!is.function(rep.covar.action))
    rep.covar.action <- match.arg(rep.covar.action, choices = c("none", "first")) 
  ##
  res <- as.geodata(obj = obj, coords.col = coords.col, data.col = data.col,
                    covar.col = covar.col, covar.names = covar.names,
                    realisations = realisations, rep.data.action = rep.data.action,
                    rep.covar.action = rep.covar.action)
  res$call <- call.fc
  return(res)
}

"as.geodata" <-
  function(obj, coords.col = 1:2, data.col = 3, data.names = NULL, 
           covar.col = NULL, covar.names = "obj.names", realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar", "none"),
           rep.data.action, rep.covar.action)
{
  if(!is.matrix(obj) & !is.data.frame(obj))
    stop("object must be a matrix or data.frame")
  if(!is.null(data.names) & length(data.col) < 2)
    stop("data.names allowed only if there is more than 1 column of data")
  res <- list()
  ##
  ## testing for NA's setting the coordinates of the data locations
  ##
  if(any(is.na(obj[,coords.col]))){
    warning("NA's not allowed in the coordinates")
    obj <- obj[complete.cases(obj),,drop = FALSE]
    warning("eliminating rows with NA's")
  }
  res$coords <- as.matrix(obj[,coords.col])
  ##
  ## setting the data
  ##
  res$data <- as.matrix(obj[,data.col])
  if(length(data.col) == 1) res$data <- as.vector(res$data)
  else if(!is.null(data.names)) colnames(res$data) <- data.names
  ##
  ## setting the covariates, if the case 
  ##
  if(!is.null(covar.col)){
    res[[3]] <- as.data.frame(obj[,covar.col])
    if(covar.names == "obj.names"){
      if(is.matrix(obj))      
        col.names <- dimnames(obj)[2]
      if(is.data.frame(obj))      
        col.names <- names(obj)
    }
    names(res)[3] <- "covariate"
    if(covar.names == "obj.names")
      if(is.null(col.names)) names(res[[3]]) <- paste("covar", 1:length(covar.col), sep="")
      else  names(res[[3]]) <- col.names[covar.col]
    else
      names(res[[3]]) <- covar.names
    covar.names <- names(res[[3]])
  }
  ##
  ## Dealing with NA's
  ##
  na.action <- match.arg(na.action)
  if(na.action != "none"){
    if(na.action == "ifany")
      na.data <- na.covar <- TRUE
    if(na.action == "ifdata")
      {na.data <- TRUE; na.covar <- FALSE}
    if(na.action == "ifcovar")
      {na.data <- FALSE; na.covar <- TRUE}
    not.na <- function(x) !any(is.na(x))
    if(na.data){
      ind <- apply(as.matrix(res$data), 1, not.na)
      if(!all(ind)){
        res$coords <- res$coords[ind,]
        res$data <- drop(as.matrix(res$data)[ind,])
        if(!is.null(covar.col))
          res[[3]] <- drop(as.matrix(res[[3]][ind,]))
        cat(paste("as.geodata:", sum(!ind), "points removed due to NA in the data\n")) 
      }
    }
    if(!is.null(covar.col) && na.covar){
      ind <- apply(as.matrix(res[[3]]), 1, not.na)
      if(!all(ind)){
        res$coords <- res$coords[ind,]
        res$data <- drop(as.matrix(res$data)[ind,])
        if(!is.null(covar.col))
          res[[3]] <- drop(res[[3]][ind,])
        cat(paste("as.geodata:", sum(!ind), "points removed due to NA in the covariate(s)\n")) 
      }
    }
  }
  ##
  ## Checking whether there are data from different realisations
  ##
  if(is.null(realisations)) realisations <- as.factor(rep(1, nrow(res$coords)))
  else{
    if(is.numeric(realisations) && length(realisations) == 1)
      realisations <- as.factor(obj[,realisations])
    res$realisations <- realisations
  }
  if(length(realisations) != nrow(res$coords))
    stop("realisations and coords have incompatible dimensions")
  ##
  ## Checking whether there are data at coincident locations
  ## and dealing with this acoording to the value of the argument
  ## rep.data.action 
  ##
  if(missing(rep.data.action)) rep.data.action <- "none"
  if(!is.function(rep.data.action))
    rep.data.action <- match.arg(rep.data.action, choices = c("none", "first")) 
  if(missing(rep.covar.action)) rep.covar.action <- rep.data.action
  if(!is.function(rep.covar.action))
    rep.covar.action <- match.arg(rep.covar.action, choices = c("none", "first")) 
  require(mva)
  rep.lev <- as.character(paste("x",res$coords[,1],"y",res$coords[,2], sep=""))
  rep.dup <- duplicated(rep.lev)
  if(sum(rep.dup) > 0)
    cat(paste("as.geodata:", sum(rep.dup), "redundant locations found\n"))
  if(is.function(rep.data.action) || rep.data.action == "first"){
    res$coords <- res$coords[!rep.dup,]
    measure.var.f <- function(x) return(summary(lm(x ~ as.factor(rep.lev)))$sigma^2)
    res$m.var <- drop(apply(as.matrix(res$data),2,measure.var.f))
    rep.action.f <- function(x, rep.action){ 
      if(!is.function(rep.action) && rep.action == "first")
        return(x[!rep.dup])
      else
        return((as.vector(by(x, rep.lev, rep.action))[codes(factor(rep.lev))])[!rep.dup])
    }
    res$data <- drop(apply(as.matrix(res$data), 2, rep.action.f, rep.action=rep.data.action))
    if(!is.null(covar.col))
      res[[3]] <- drop(apply(res[[3]], 2, rep.action.f, rep.action=rep.covar.action))
    if(!is.null(res$realisations))
      res$realisations <- res$realisations[!rep.dup]
  }
  else{
    check.coincide <- function(x){sum(dist(x) < 1e-12) > 0}
    any.coincide <- lapply(split(as.data.frame(res$coords), realisations), check.coincide)
    any.coincide <- as.vector(unlist(any.coincide))
    if(sum(any.coincide) > 0)
      cat("WARNING: there are data at coincident locations, some of the geoR's functions will not work.\n")
  }
  ##
  if(!is.null(covar.col)){
    res[[3]] <- as.data.frame(res[[3]])
    names(res[[3]]) <- covar.names
  }
  class(res) <- "geodata"
  return(res)
}

"summary.geodata" <- function(object, ...)
{
  x <- object
  res <- list()
  res$coords.summary <- apply(x$coords, 2, range)
  rownames(res$coords.summary) <- c("min", "max")
  require(mva)
  res$distances.summary <- range(dist(x$coords))
  names(res$distances.summary) <- c("min", "max")  
  if(!is.null(x$borders)){
    res$borders.summary <- apply(x$borders, 2, range)
    rownames(res$borders.summary) <- c("min", "max")
  }
  res$data.summary <- drop(apply(as.matrix(x$data), 2, summary))
  if(!is.null(x$units.m))
    res$units.m.summary <- drop(apply(as.matrix(x$units.m), 2, summary))
  if(!is.null(x$covariate))
    res$covariate.summary <- summary(x$covariate)
  return(res)
}

"print.summary.geodata" <- function(x, ...)
{
  cat("Coordinates summary\n")
  print(x$coords.summary)
  if(!is.null(x$borders.summary)){
    cat("\nBorders summary\n")
    print(x$borders.summary)
  }
  cat("\nData summary\n")
  print(x$data.summary)
  if(!is.null(x$units.m.summary)){
    cat("\nOffset variable summary\n")
    print(x$units.m.summary)
  }
  if(!is.null(x$covariate.summary)){
    cat("\nCovariates summary\n")
    print(x$covariate.summary)
  }
  return(invisible())
}

"points.geodata" <-
  function (x, coords = x$coords, data = x$data, 
            data.col = 1, borders = NULL,
            pt.divide = c("data.proportional",
              "rank.proportional", "quintiles",
              "quartiles", "deciles", "equal"),
            lambda=1, trend="cte", weights.divide=NULL,
            cex.min, cex.max, pch.seq, col.seq, add.to.plot = FALSE,
            round.quantiles = FALSE, graph.pars = FALSE, ...) 
{
  if(missing(x)) x <- list(coords = coords, data = data)
  # This is for compatibility with previously used argument pt.sizes
  if(!is.null(list(...)$pt.s)) pt.divide <- list(...)$pt.s
  #
  if(!is.numeric(pt.divide))
    pt.divide <- match.arg(pt.divide)
  if(!is.vector(data))
       data <- (as.data.frame(data))[,data.col]
  if(nrow(coords) != length(data))
    stop("coords and data have incompatible sizes")
    if (!is.null(weights.divide)) {
    if (length(weights.divide) != length(data)) 
      stop("length of weights.divide must be equals to the length of data")
    data <- data/weights.divide
  }
  ##
  ## data transformation (Box-Cox)
  ##
  if (lambda != 1) data <- BCtransform(data, lambda)$data
  ##
  ## trend removal
  ##
  xmat <- unclass(trend.spatial(trend = trend, geodata = x))
  if (nrow(xmat) != nrow(coords)) 
    stop("coords and trend have incompatible sizes")
  if (trend != "cte") {
    data <- lm(data ~ xmat + 0)$residuals
    names(data) <- NULL
  }
  ##
  attach(x)
  eval(borders)
  detach(x)
  if (!add.to.plot) {
    if(is.null(borders))
      coords.lims <- set.coords.lims(coords=coords)
    else{
      if(ncol(borders) != 2)
        stop("argument borders must be an object with 2 columns with the XY coordinates of the borders of the area")
      coords.lims <- set.coords.lims(coords=rbind(as.matrix(coords), as.matrix(borders)))
    }
    par(pty = "s")
    toplot <- apply(coords, 2, range)
    colnames(toplot) <- c("X Coord", "Y Coord")
    plot(toplot, type = "n",
         xlim = coords.lims[,1], ylim = coords.lims[, 2], ...)
    if(!is.null(borders))
      polygon(borders)
  }
  if (missing(cex.min)) cex.min <- 0.5
  if (missing(cex.max)) cex.max <- 1.5
  graph.list <- list()
  if(is.numeric(pt.divide)|pt.divide == "quintiles"|pt.divide == "quartiles"|pt.divide == "deciles") {
    if (pt.divide == "quintiles") {
      n.quant <- 5
      if (missing(col.seq)) 
        col.seq <- c("blue", "green", "yellow", "orange3", "red2")
    }
    if (pt.divide == "quartiles") {
      n.quant <- 4
      if (missing(col.seq)) 
        col.seq <- c("blue", "green", "yellow", "red") 
    }
    if (pt.divide == "deciles") {
      n.quant <- 10
      if (missing(col.seq)) 
        col.seq <- rainbow(13)[10:1]
    }
    if(is.numeric(pt.divide)){
      data.quantile <- pt.divide
      n.quant <- length(pt.divide) - 1
    }
    else
      data.quantile <- quantile(data, probs = seq(0, 1, by = (1/n.quant)))
    if(!missing(col.seq) && all(col.seq == "gray")) col.seq <- gray(seq(1,0, l=n.quant))
    if (missing(pch.seq)) pch.seq <- rep(21, n.quant)
    cex.pt <- seq(cex.min, cex.max, l = n.quant)
    if (round.quantiles == TRUE) {
      data.quantile[1] <- floor(data.quantile[1])
      data.quantile[n.quant + 1] <- ceiling(data.quantile[n.quant + 1])
      data.quantile <- round(data.quantile)
    }
    graph.list$quantiles <- data.quantile
    graph.list$cex <- cex.pt
    graph.list$col <- col.seq
    graph.list$pch <- pch.seq
    graph.list$data.group <- cut(data, breaks=data.quantile, include.l=TRUE)
    if (add.to.plot) 
      points(coords, pch = pch.seq, cex = cex.pt[as.numeric(graph.list$data.group)], bg = col.seq[as.numeric(graph.list$data.group)], ...)
    else
      points(coords, pch = pch.seq, cex = cex.pt[as.numeric(graph.list$data.group)], bg = col.seq[as.numeric(graph.list$data.group)])
  }
  else {
    n <- length(data)
    if (missing(pch.seq)) pch.seq <- 21
    if (missing(col.seq)) col.seq <- 0
    else if(all(col.seq == "gray")) col.seq <- gray(seq(1,0, l=n))
    coords.order <- coords[order(data), ]
    data.order <- data[order(data)]
    if (pt.divide == "rank.proportional") {
      data.quantile <- range(data.order)
      size <- seq(cex.min, cex.max, l = n)
      graph.list$cex <- range(size)
      graph.list$pch <- unique(range(pch.seq))
      graph.list$col <- col.seq
      if (length(col.seq) == 1) col.seq <- rep(col.seq, n)
      for (i in 1:n) {
        if (add.to.plot) 
          points(coords.order[i, , drop = FALSE], cex = size[i], 
                 pch = pch.seq, bg = col.seq[i], ...)
        else points(coords.order[i, , drop = FALSE], 
                    cex = size[i], pch = pch.seq, bg = col.seq[i])
      }
    }
    if (pt.divide == "data.proportional") {
      r.y <- range(data.order)
      size <- cex.min + ((data.order - r.y[1]) * (cex.max - 
                                                  cex.min))/(r.y[2] - r.y[1])
      graph.list$cex <- c(cex.min, cex.max)
      graph.list$pch <- unique(range(pch.seq))
      graph.list$col <- col.seq
      if (length(col.seq) == 1) col.seq <- rep(col.seq, n)
      for (i in 1:n) {
        if (add.to.plot) 
          points(coords.order[i, , drop = FALSE], cex = size[i], 
                 pch = pch.seq, bg = col.seq[i], ...)
        else points(coords.order[i, , drop = FALSE], 
                    cex = size[i], pch = pch.seq, bg = col.seq[i])
      }
    }
    if (pt.divide == "equal") {
      if (add.to.plot) 
        points(coords, pch = pch.seq, bg = col.seq, cex = cex.max, 
               ...)
      else points(coords, pch = pch.seq, bg = col.seq, 
                  cex = cex.max)
    }
  }
  if (graph.pars == TRUE) 
    return(graph.list)
  else return(invisible())
}

plot.geodata <- function (x, coords = x$coords, data = x$data, borders = NULL, 
    trend = "cte", lambda = 1, col.data = 1, weights.divide = NULL, 
    scatter3d = FALSE, ...) 
{
  if(missing(x)) x <- list(coords=coords, data = data)
  if (is.R()) par.ori <- par(no.readonly = TRUE)
  else par.ori <- par()
  on.exit(par(par.ori))
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  data <- data[, col.data]
  attach(x)
  eval(borders)
  detach(x)
  if (!is.null(weights.divide)) {
    if (length(weights.divide) != length(data)) 
      stop("length of weights.divide must be equals to the length of data")
    data <- data/weights.divide
  }
  if (lambda != 1) {
    if (lambda == 0) 
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  xmat <- unclass(trend.spatial(trend = trend, geodata = x))
  if (nrow(xmat) != nrow(coords)) 
    stop("coords and trend have incompatible sizes")
  if (trend != "cte") {
    data <- lm(data ~ xmat + 0)$residuals
    names(data) <- NULL
  }
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 0, 0.5))
  data.quantile <- quantile(data)
  if (is.null(borders)) 
    coords.lims <- set.coords.lims(coords = coords)
  else {
    if (ncol(borders) != 2) 
      stop("argument \"borders\" must be a 2 columns object with coordinates of the borders of the study area")
    coords.lims <- set.coords.lims(coords = rbind(as.matrix(coords), as.matrix(borders)))
  }
  par(pty = "s")
  plot(coords, xlab = "X Coord", ylab = "Y Coord ", type = "n", 
       xlim = coords.lims[, 1], ylim = coords.lims[, 2])
  if (!is.null(borders)) polygon(borders)
  if (is.R()) {
    data.breaks <- unique(quantile(data))
    n.breaks <- length(data.breaks)
    data.cut <- cut(data, breaks = data.breaks, include.l = TRUE, 
                    labels = FALSE)
    points(coords, pch = (1:4)[data.cut], col = c("blue", 
                                            "green", "yellow2", "red")[data.cut])
  }
  else {
    points(coords[(data <= data.quantile[2]), ], pch = 1, 
           cex = 0.6, col = 2)
    points(coords[((data > data.quantile[2]) & (data <= data.quantile[3])), 
                  ], pch = 2, cex = 1.4, col = 4)
    points(coords[((data > data.quantile[3]) & (data <= data.quantile[4])), 
                  ], pch = 3, cex = 1.7, col = 7)
    points(coords[(data > data.quantile[4]), ], pch = 4, 
           cex = 2, col = 8)
  }
  plot(data, coords[, 2], ylab = "Coord Y", cex = 1, ylim = coords.lims[, 2])
  if (!is.R()) 
    par(mar = c(5, 5, 1, 0.5))
  plot(coords[, 1], data, xlab = "Coord X", cex = 1, xlim = coords.lims[, 1], )
  par(pty = "m")
  if (is.R()) par(mar = c(4, 4, 1, 1))
  else par(mar = c(0, 1, 0, 0.5))
  if (scatter3d) {
    if (!require(scatterplot3d)) {
      cat("plot.geodata: the argument scatter3d=TRUE requires the package \"scatterplot3d\" \n              will plot an histogram instead")
      hist(data)
    }
    else scatterplot3d(x = coords[, 1], y = coords[, 2], 
                       z = data, box = FALSE, type = "h", pch = 16, xlab = "Coord X", 
                       ylab = "Coord Y", ...)
  }
  ##  else xyzplot(coords = coords, data = data, ...)
  else hist(data, main="", ...)
  return(invisible())
}
