"trend.spatial" <-
  function (trend, geodata) 
{
  if(!missing(geodata)){
    attach(geodata, pos=2)
    if(!is.null(geodata$covariate)){
      attach(geodata$covariate, pos=3)
      on.exit(detach("geodata$covariate"), add=TRUE)
    }
    on.exit(detach("geodata"), add=TRUE)
  }
  if (inherits(trend, "formula")) {
    trend.mat <- try(model.matrix(trend))
    if (inherits(trend.mat, "try-error")) 
      stop("\ntrend elements not found")
  }
  else {
    if(is.numeric(trend))
      trend.mat <- unclass(trend)
    else if (trend == "cte"){
      if(missing(geodata))
        stop("argument geodata must be provided with trend=\"cte\"")
      trend.mat <- as.matrix(rep(1, nrow(geodata$coords)))
    }
    else if (trend == "1st"){
      if(missing(geodata))
        stop("argument geodata must be provided with trend=\"1st\"")
      trend.mat <- cbind(1, geodata$coords)
    }
    else if (trend == "2nd"){ 
      if(missing(geodata))
        stop("argument geodata must be provided with trend=\"2nd\"")
      trend.mat <- cbind(1, geodata$coords, geodata$coords[,1]^2,
                         geodata$coords[,2]^2,
                         geodata$coords[,1] * geodata$coords[,2])
    }
    else stop("external trend must be provided for data locations to be estimated using the arguments trend.d and trend.l. Allowed values are \"cte\", \"1st\", \"2nd\" or  a model formula")
  }
  trend.mat <- as.matrix(trend.mat)
  dimnames(trend.mat) <- list(NULL, NULL)
  class(trend.mat) <- "trend.spatial"
  return(trend.mat)
}

#"trend.spatial" <-
#  function (trend, geodata) 
#{
#  print(search())
#  if(missing(geodata)) geodata <- list()
#  if (inherits(trend, "formula")) {
#    trend.frame <- geodata[attr(terms(trend), "term.labels")]
#    if(any(is.null(names(trend.frame)) || is.na(names(trend.frame)))){
#      trend.mat <- data.frame()
#      class(trend.mat) <- "try-error"
#    }
#    else
#      trend.mat <- try(model.matrix(trend, data = trend.frame))
#    if (inherits(trend.mat, "try-error")) {
#      if (!is.null(geodata$covariate))
#        trend.mat <-
#          try(model.matrix(trend, 
#                           data = as.data.frame(geodata$covariate)))
#      if (inherits(trend.mat, "try-error"))
#      	if(is.data.frame(geodata$data) | is.list(geodata$data))
#          trend.mat <- try(model.matrix(trend, data = as.data.frame(geodata$data)))
#      if (inherits(trend.mat, "try-error")) 
#        trend.mat <- try(model.matrix(trend))
#      if (inherits(trend.mat, "try-error")) 
#        stop("\ntrend elements not found")
#    }
#  }
#  else {
#    if(is.numeric(trend))
#      trend.mat <- unclass(trend)
#    else if (trend == "cte") 
#      trend.mat <- as.matrix(rep(1, nrow(geodata$coords)))
#    else if (trend == "1st") 
#      trend.mat <- cbind(1, geodata$coords)
#    else if (trend == "2nd") 
#      trend.mat <- cbind(1, geodata$coords, geodata$coords[,1]^2,
#                         geodata$coords[,2]^2,
#                         geodata$coords[,1] * geodata$coords[,2])
#    ## This would use orthogonal polynomials:
#    ##      trend.mat <- poly(as.matrix(geodata$coords), degree = 2)
#    else stop("external trend must be provided for data locations to be estimated using the arguments trend.d and trend.l. Allowed values are \"cte\", \"1st\", \"2nd\" or  a model formula")
#  }
#  trend.mat <- as.matrix(trend.mat)
#  dimnames(trend.mat) <- list(NULL, NULL)
#  class(trend.mat) <- "trend.spatial"
#  return(trend.mat)
#}

"locations.inside" <-
  function(locations, borders)
{
  if(is.list(borders))
    borders <- matrix(unlist(borders[1:2], ncol=2))
  borders <- as.matrix(borders)
  if(ncol(borders) != 2)
    stop("borders must be a matrix or data-frame with two columns")
  if (require(splancs) == FALSE)
    cat("package splancs in required to select points inside the borders\n")
  locations <- locations[as.vector(inout(pts = locations,
                                         poly = borders)),]
  return(locations)
}


"likfit.nospatial" <-
  function(temp.list, ...)
{
  results <- list()
  z <- temp.list$z
  n <- temp.list$n
  beta.size <- temp.list$beta.size
  xmat <- temp.list$xmat
  txmat <- temp.list$txmat
  ixx <- solve(crossprod(xmat))
  if(temp.list$fix.lambda == FALSE){
    if (temp.list$minimisation.function == "nlm"){
      assign(".temp.lower.lambda",-2, pos=1)
      assign(".temp.upper.lambda", 2, pos=1)
      results <- nlm(proflik.lambda, 1, ...)
      if(exists(".temp.lambda")){
        results$lambda <- .temp.lambda
        remove(".temp.lambda", pos=1, inherits = TRUE)
      }
      else{
        results$lambda <- results$estimate
      }
      rm(.temp.lower.lambda, .temp.upper.lambda, inherits = TRUE, pos=1)
    }
    if (temp.list$minimisation.function == "nlmP"){
      results <- nlmP(proflik.lambda, 1, lower=-2, upper=2,...)  
      results$lambda <- results$estimate
    }
    if (temp.list$minimisation.function == "optim"){
      results <- optim(1, proflik.lambda, method="L-BFGS-B", lower=-2, upper=2,...)
      results$minimum <- results$value
      results$lambda <- results$par
    }
    if(results$lambda == 1) {
      temp.list$log.jacobian <- 0
    }
    else {
      if(any(z <= 0))
        stop("Transformation option not allowed when there are zeros or negative data")
      if(any(z^(results$lambda - 1) <= 0))
        temp.list$log.jacobian <- log(prod(z^(results$lambda - 1)))
      else temp.list$log.jacobian <- sum(log(z^(results$lambda - 1)))
      if(results$lambda == 0)
        z <- log(z)
      else z <- ((z^results$lambda) - 1)/results$lambda
    }
  }
  else{
    results$lambda <- temp.list$lambda
    results$code <- 1
    if (temp.list$minimisation.function == "optim") results$convergence <- 0
  }
  ssres <- (z %*% (diag(n) - xmat %*%
                   solve(crossprod(xmat)) %*% txmat) %*% z)
  if(temp.list$method == "ML"){
    results$tausqhat <- ssres/n
    if(temp.list$fix.lambda)
      results$minimum <- as.vector(((n/2) * log(2 * pi) +
                          (n/2) * log(results$tausqhat) +
                          (n/2)  -
                          temp.list$log.jacobian))
  }
  if(temp.list$method == "RML") {
    results$tausqhat  <- (ssres/(n-beta.size))
    if(temp.list$fix.lambda){
      results$minimum <- as.vector((((n - beta.size)/2) * log(2 * pi) +
                          ((n - beta.size)/2) * log(results$tausqhat) +
                          (n/2) -
                          temp.list$log.jacobian
                          ))
    }
  }
  if (temp.list$minimisation.function == "optim") results$value <- results$minimum    
  return(results)
}

"loglik.spatial" <- function(pars)
{
  tausq <- pars[1]
  sigmasq <- pars[2]
  sill.total <- tausq + sigmasq
  phi <- pars[3]
  lambda <- pars[4]
  z <- .temp.list$z
  n <- .temp.list$n
  if(.temp.list$fix.lambda == FALSE) {
    if(lambda == 1) {
      .temp.list$log.jacobian <- 0
    }
    else {
      if(any(z < 0))
        stop("Transformation option not allowed when there are zeros or negative data"
             )
      if(any(z^(lambda - 1) <= 0))
        .temp.list$log.jacobian <- log(prod(z^(lambda -
                                               1)))
      else .temp.list$log.jacobian <- sum(log(z^(lambda - 1)))
      if(lambda == 0)
        z <- log(z)
      else z <- ((z^lambda) - 1)/lambda
    }
  }
  beta.size <- .temp.list$beta.size
  kappa <- .temp.list$kappa
  covinf <- varcov.spatial(dists.lowertri = .temp.list$
                           dists.lowertri, cov.model = .temp.list$cov.model,
                           kappa = kappa, nugget = tausq,
                           cov.pars = c(sigmasq, phi), scaled = FALSE,
                           inv = TRUE, det = TRUE,
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
            res = xix, PACKAGE = "geoR")$res
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
            res = xiy, PACKAGE = "geoR")$res
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
            res = yiy, PACKAGE = "geoR")$res
  ssresmat <- as.vector(yiy - 2*crossprod(beta.hat,xiy) +  beta.hat %*% xix %*% beta.hat)
  if(.temp.list$method == "ML") {
    loglik <- ( - (n/2) * log(2 * pi) -
               covinf$log.det.to.half -
               0.5 * ssresmat + 
               .temp.list$log.jacobian)
  }
  if(.temp.list$method == "RML") {
    xx.eigen <- eigen(crossprod(.temp.list$xmat), symmetric = TRUE, only.values = TRUE)
    loglik <- ( - ((n - beta.size)/2) * log(2 * pi) +
               0.5 * sum(log(xx.eigen$values)) -
               covinf$log.det.to.half -
               (0.5) * ssresmat -
               choldet +
               .temp.list$log.jacobian)
  }
  return(as.vector(loglik))
}

"matern" <-
  function (u, phi, kappa) 
{
  if(is.vector(u)) names(u) <- NULL
  if(is.matrix(u)) dimnames(u) <- list(NULL, NULL)
  uphi <- u/phi
  uphi <- ifelse(u > 0,
                 (((2^(-(kappa-1)))/gamma(kappa)) *
                  (uphi^kappa) *
                  besselK(x=uphi, nu=kappa)), 1)    
  uphi[u > 600*phi] <- 0 
  return(uphi)
}


"nlmP" <-
  function(objfunc, params, lower = rep( -Inf, length(params)),
           upper = rep(+Inf, length(params)), ... )
{
  ## minimizer, using nlm with transformation of variables
  ##
  ## objfunc is a function to be optimised
  ## params is a starting value for the parameters
  Nparams <- length(params)
  if(length(lower) != Nparams)
    stop(" lower boundry different length than params")
  if(length(upper) != Nparams)
    stop(" upper boundry different length than params")
  checklimits <- upper - lower
  if(any(checklimits <= 0))
    stop(" bad boundries")
  if(any(params < lower))
    stop(" starting params too low")
  if(any(params > upper))
    stop(" starting params too high")
  
  bothlimQQ <- (lower != (-Inf)) & (upper != +Inf)
  loweronlyQQ <- (lower != (-Inf)) & (upper == +Inf)
  upperonlyQQ <- (lower == (-Inf)) & (upper != +Inf)
  ubothQQ <- upper[bothlimQQ]
  lbothQQ <- lower[bothlimQQ]
  dbothQQ <- ubothQQ - lbothQQ
  loneQQ <- lower[loweronlyQQ]
  uoneQQ <- upper[upperonlyQQ]
  
  .bounds.list <- list(bothlimQQ = bothlimQQ, 
                       loweronlyQQ = loweronlyQQ,
                       upperonlyQQ = upperonlyQQ,
                       ubothQQ = ubothQQ,
                       lbothQQ = lbothQQ,
                       dbothQQ = dbothQQ,
                       loneQQ = loneQQ,
                       uoneQQ = uoneQQ)
  
  assign(".objfuncQQ", objfunc, pos=1)
  assign(".bounds.list", .bounds.list, pos=1)
  
  ## reduce the parameter space by a scale to keep parameters
  ## away from the boundries
  
  normaltomad <- function(normalparamsX)
    {
      madparamsX <- normalparamsX
      if(any(.bounds.list$bothlimQQ)) {
        noughtone <- (normalparamsX[.bounds.list$bothlimQQ] -
                      .bounds.list$lbothQQ)/.bounds.list$dbothQQ
        madparamsX[.bounds.list$bothlimQQ] <- log(noughtone/(1 - noughtone))
      }
      
      if(any(.bounds.list$loweronlyQQ))
        madparamsX[.bounds.list$loweronlyQQ] <-
          log(normalparamsX[.bounds.list$loweronlyQQ] - .bounds.list$loneQQ)
      
      if(any(.bounds.list$upperonlyQQ))
        madparamsX[.bounds.list$upperonlyQQ] <-
          log(.bounds.list$uoneQQ - normalparamsX[.bounds.list$upperonlyQQ])
      
      return(madparamsX)
    }
  
  madtonormalQQ <<- function(madparamsX)
    {
      normalparamsX <- madparamsX
      
      if(any(.bounds.list$bothlimQQ)) {
###        madparamsX[((.bounds.list$bothlimQQ) & (madparamsX > 300))] <- 300
        emad <- exp(madparamsX[.bounds.list$bothlimQQ])
        normalparamsX[.bounds.list$bothlimQQ] <-
          .bounds.list$dbothQQ * (emad/(1 + emad)) + .bounds.list$lbothQQ
      }
      
      if(any(.bounds.list$loweronlyQQ)){
        normalparamsX[.bounds.list$loweronlyQQ] <-
          exp(madparamsX[.bounds.list$loweronlyQQ]) + .bounds.list$loneQQ
      }
      
      if(any(.bounds.list$upperonlyQQ))
        normalparamsX[.bounds.list$upperonlyQQ] <-
          - exp(madparamsX[.bounds.list$upperonlyQQ]) + .bounds.list$uoneQQ
      
      if(exists(".ind.prof.phi"))
        if(is.nan(normalparamsX[.ind.prof.phi]))
          normalparamsX[.ind.prof.phi] <- 0
      
      return(normalparamsX)
    }
  
  newobjfunc <- function(madparams) {
    normalparams <-  madtonormalQQ(madparams)
    
    .objfuncQQ(normalparams)
    
  }
  
  startmadparams <- normaltomad(params)
  result <- nlm(newobjfunc, startmadparams, ...)
  result$madestimate <- result$estimate
  result$estimate <- madtonormalQQ(result$madestimate)
  remove(".bounds.list", pos=1, inherits=TRUE)
  remove(".objfuncQQ", pos=1, inherits=TRUE)
  remove("madtonormalQQ", pos=1, inherits=TRUE)
  
###  return(result, madtonormalQQ(normaltomad(params)),params)
  return(result)
}



"dinvchisq" <-
  function(x, df, scale=1/df, log = FALSE)
{
  if(df <= 0)
    stop("df must be greater than zero")
  if(scale <= 0)
    stop("scale must be greater than zero")
  nu <- df/2
  if(log)
    return(ifelse(x > 0, nu*log(nu) - log(gamma(nu)) + nu*log(scale) -
                  (nu+1)*log(x) - (nu*scale/x), NA))
  else
    return(ifelse(x > 0,
                  (((nu)^(nu))/gamma(nu)) * (scale^nu) *
                  (x^(-(nu+1))) * exp(-nu*scale/x), NA))
}

"rinvchisq" <- 
  function (n, df, scale = 1/df)
{
  if((length(scale)!= 1) & (length(scale) != n))
    stop("scale should be a scalar or a vector of the same length as x")
  if(df <= 0)
    stop("df must be greater than zero")
  if(any(scale <= 0))
    stop("scale must be greater than zero")
  return((df*scale)/rchisq(n, df=df)) 
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
  if(missing(x))
    x <- list(coords = coords, data = data)
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
  if (lambda != 1)
    data <- BCtransform(data, lambda)$data
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
  if (add.to.plot == FALSE) {
    if(is.null(borders))
      coords.lims <- set.coords.lims(coords=coords)
    else{
      if(ncol(borders) != 2)
        stop("argument borders must have 2 columns with the XY coordinates of the borders of the area")
      coords.lims <- set.coords.lims(coords=rbind(coords, borders))
    }
    par(pty = "s")
    plot(apply(coords, 2, range), type = "n",
         xlim = coords.lims[,1], ylim = coords.lims[, 2], ...)
    if(!is.null(borders))
      polygon(borders)
  }
  if (missing(cex.min)) 
    cex.min <- 0.5
  if (missing(cex.max)) 
    cex.max <- 1.5
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
    if (missing(pch.seq)) 
      pch.seq <- rep(21, n.quant)
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
    if (missing(pch.seq)) 
      pch.seq <- 21
    if (missing(col.seq)) 
      col.seq <- 0
    n <- length(data)
    coords.order <- coords[order(data), ]
    data.order <- data[order(data)]
    if (pt.divide == "rank.proportional") {
      data.quantile <- range(data.order)
      size <- seq(cex.min, cex.max, l = n)
      graph.list$cex <- range(size)
      graph.list$pch <- unique(range(pch.seq))
      graph.list$col <- col.seq
      if (length(col.seq) == 1) 
        col.seq <- rep(col.seq, n)
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
      if (length(col.seq) == 1) 
        col.seq <- rep(col.seq, n)
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

