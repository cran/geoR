"variog" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            uvec = "default", trend = "cte", lambda = 1,
            option = c("bin", "cloud", "smooth"),
            estimator.type = c("classical", "modulus"), 
            nugget.tolerance = 0, max.dist = NULL, pairs.min = 2,
            bin.cloud = FALSE, direction = "omnidirectional", tolerance = pi/8,
            unit.angle = c("radians","degrees"), messages.screen = TRUE, ...) 
{
  if (is.R()){
    require(mva)
    require(modreg)
  }
  call.fc <- match.call()
  ##
  ## Directional variogram
  ##
  unit.angle <- match.arg(unit.angle)
  if(is.numeric(direction)){
    if(length(direction) > 1)
      stop("only one direction is allowed")
    if(length(tolerance) > 1)
      stop("only one tolerance value is allowed")
    if(unit.angle == "degrees"){
      ang.deg <- direction
      ang.rad <- (ang.deg * pi)/180
      tol.deg <- tolerance
      tol.rad <- (tol.deg * pi)/180
    }
    else{
      ang.rad <- direction
      ang.deg <- (ang.rad * 180)/pi
      tol.rad <- tolerance
      tol.deg <- (tol.rad * 180)/pi
    }
    if(ang.rad > pi | ang.rad < 0)
      stop("direction must be an angle in the interval [0,pi] radians")
    if(tol.rad > pi/2 | tol.rad < 0)
      stop("tolerance must be an angle in the interval [0,pi/2] radians")
    if(tol.deg >= 90){
      direction <- "omnidirectional"
      cat("variog: computing omnidirectional variogram\n")
    }(0, ymax, legend=c(substitute(a*degree, list(a=n.o[1])),
                        substitute(a*degree, list(a=n.o[2])),
                        substitute(a*degree, list(a=n.o[3])),
                        substitute(a*degree, list(a=n.o[4])),
                        "omnidirectional", expression()),...)
    else{
      if(messages.screen){
        cat(paste("variog: computing variogram for direction = ", round(ang.deg, dig=3), " degrees (", round(ang.rad, dig=3), " radians)\n", sep=""))
        cat(paste("        tolerance angle = ", round(tol.deg, dig=3), " degrees (", round(tol.rad, dig=3), " radians)\n", sep=""))
      }
    }
  }
  else
    if(messages.screen)
      cat("variog: computing omnidirectional variogram\n")
  ##   
  ##
  coords <- as.matrix(coords)
  data <- as.matrix(data)
  data.var <- apply(data, 2, var)
  n.data <- nrow(coords)
  n.datasets <- ncol(data)
  if (ncol(data) == 1) 
    data <- as.vector(data)
  ##
  ## variogram estimator
  ##
  option <- match.arg(option)
  if (estimator.type == "robust") 
    estimator.type <- "modulus"
  estimator.type <- match.arg(estimator.type)
  if (estimator.type == "modulus") 
    estimator.type <- "robust"
  if (lambda != 1) {
    if (lambda == 0) 
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  ##
  ## trend removal
  ##
  xmat <- trend.spatial(trend = trend, coords = coords)
  if (trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ x + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  ##
  ## 
  ##
  u <- as.vector(dist(as.matrix(coords)))
  if(direction != "omnidirectional"){
    u.ang <- .C("tgangle",
                as.double(as.vector(coords[,1])),
                as.double(as.vector(coords[,2])),
                as.integer(dim(coords)[1]),
                res = as.double(rep(0, length(u))))$res
    u.ang <- atan(u.ang)
    u.ang[u.ang < 0] <- u.ang[u.ang < 0] + pi
  }
  if (option == "bin" & bin.cloud == FALSE & direction == "omnidirectional") {
    if (!is.null(max.dist)) 
      umax <- max(u[u < max.dist])
    else umax <- max(u)
    if (all(uvec == "default")) 
      uvec <- seq(0, umax, l = 15)
    ubin <- c(0, uvec)
    nvec <- length(ubin)
    d <- 0.5 * diff(ubin[2:nvec])
    bins.lim <- c(0, (ubin[2:(nvec - 1)] + d), (d[nvec - 
                                                  2] + ubin[nvec]))
    if (uvec[1] == 0 & nugget.tolerance == 0) 
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (nugget.tolerance > 0) {
      bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim > 
                                                  nugget.tolerance])
      uvec <- c(0, (bins.lim[-(1:2)] - 0.5 * diff(bins.lim)[-1]))
    }
    nbins <- length(bins.lim) - 1
    if (is.null(max.dist)) 
      max.dist <- max(bins.lim)
    bin.f <- function(data) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      result <- .C("binit", as.integer(n.data),
                   as.double(as.vector(coords[, 1])),
                   as.double(as.vector(coords[, 2])), as.double(as.vector(data)), 
                   as.integer(nbins), as.double(as.vector(bins.lim)), 
                   as.integer(estimator.type == "robust"), as.double(max.dist), 
                   cbin = as.integer(cbin), vbin = as.double(vbin), 
                   as.integer(TRUE), sdbin = as.double(sdbin))[c("vbin", 
                                       "cbin", "sdbin")]
    }
    result <- array(unlist(lapply(as.data.frame(data), bin.f)), 
                    dim = c(nbins, 3, n.datasets))
    indp <- (result[, 2, 1] >= pairs.min)
    result <- list(u = uvec[indp], v = result[indp, 1, ], 
                   n = result[indp, 2, 1], sd = result[indp, 3, ], bins.lim = bins.lim, 
                   ind.bin = indp)
  }
  else {
    if (is.matrix(data)) {
      v <- matrix(0, nrow = length(u), ncol = ncol(data))
      for (i in 1:ncol(data)) {
        v[, i] <- as.vector(dist(data[, i]))
        if (estimator.type == "robust") 
          v[, i] <- v[, i]^(0.5)
        else v[, i] <- (v[, i]^2)/2
      }
      if (!is.null(max.dist)) {
        v <- v[u <= max.dist, ]
        if(direction != "omnidirectional")
          u.ang <- u.ang[u <= max.dist]
        u <- u[u <= max.dist]
      }
      if(direction != "omnidirectional"){
        ang.ind <- ((u.ang >= ang.rad - tol.rad) & (u.ang <= ang.rad + tol.rad))
        v <- v[ang.ind,]
        u <- u[ang.ind]
      }
    }
    else {
      v <- as.vector(dist(data))
      if (estimator.type == "robust") 
        v <- v^(0.5)
      else v <- (v^2)/2
      if (is.numeric(max.dist)) {
        v <- v[u <= max.dist]
        if(direction != "omnidirectional")
          u.ang <- u.ang[u <= max.dist]
        u <- u[u <= max.dist]
      }
      if(direction != "omnidirectional"){
        ang.ind <- ((u.ang >= ang.rad - tol.rad) & (u.ang <= ang.rad + tol.rad))
        v <- v[ang.ind]
        u <- u[ang.ind]
      }
    }
    if (option == "cloud") {
      result <- list(u = u, v = v)
    }
    if (option == "bin") {
      if (!is.null(max.dist)) 
        umax <- max(u[u < max.dist])
      else umax <- max(u)
      result <- rfm.bin(cloud = list(u = u, v = v),
                        estimator.type = estimator.type, 
                        uvec = uvec, nugget.tolerance = nugget.tolerance, 
                        bin.cloud = bin.cloud, max.dist = umax, ...)
      indp <- rep(TRUE, length(u))
      if (pairs.min > 0) {
        indp <- (result$n >= pairs.min)
        if (is.matrix(result$v)) {
          result$v <- result$v[indp, ]
          result$sd <- result$sd[indp, ]
        }
        else {
          result$v <- result$v[indp]
          result$sd <- result$sd[indp]
        }
        result$u <- result$u[indp]
        result$n <- result$n[indp]
      }
      result$ind.bin <- indp
    }
    if (option == "smooth") {
      if (is.R()) 
        require(modreg)
      if (is.matrix(v)) 
        stop("smooth not yet available for several variables")
      temp <- ksmooth(u, v, ...)
      result <- list(u = temp[[1]], v = temp[[2]])
    }
  }
  result <- c(result, list(var.mark = data.var, output.type = option, 
                           estimator.type = estimator.type, n.data = n.data,
                           lambda = lambda, trend = trend))
  result$nugget.tolerance <- nugget.tolerance
  if(direction != "omnidirectional") result$direction <- ang.rad
  else result$direction <- "omnidirectional"
  if(direction != "omnidirectional") result$tolerance <- tol.rad
  else result$tolerance <- "none" 
  result$uvec <- uvec
  result$call <- call.fc
  class(result) <- "variogram"
  return(result)
}

"variog4" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            uvec = "default", trend = "cte", lambda = 1,
            option = c("bin", "cloud", "smooth"),
            estimator.type = c("classical", "modulus"), 
            nugget.tolerance = 0, max.dist = NULL, pairs.min = 2,
            bin.cloud = FALSE, direction = c(0, pi/4, pi/2, 3*pi/4), tolerance = pi/8,
            unit.angle = c("radians", "degrees"), messages.screen = TRUE, ...) 
{
  if(length(direction) != 4)
    stop("argument direction must be a vector with 4 values. For different specifications use the functio variog()")
  if(length(tolerance) != 1)
    stop("only 1 values can be provided to the argument tolerance . For different specifications use the functio variog()")
  res <- list()
  if(unit.angle == "radians")
    dg <- direction * 180/pi
  else dg <- direction
  for(angle in direction){
    res[[as.character(round(dg[which(direction == angle)], dig=1))]] <-
      variog(coords=coords, data=data,
             uvec=uvec, trend = trend,
             lambda = lambda, option = option,
             estimator.type = estimator.type,
             nugget.tolerance = nugget.tolerance,
             max.dist = max.dist,
             pairs.min = pairs.min,
             bin.cloud = bin.cloud,
             direction = angle,
             tolerance = tolerance,
             unit.angle = unit.angle,
             messages.screen = TRUE, ...)
  }
  res$omnidirectional <- variog(coords=coords, data=data,
                                uvec=uvec, trend = trend,
                                lambda = lambda, option = option,
                                estimator.type = estimator.type,
                                nugget.tolerance = nugget.tolerance,
                                max.dist = max.dist,
                                pairs.min = pairs.min,
                                bin.cloud = bin.cloud,
                                direction = "omnidirectional",
                                tolerance = tolerance,
                                unit.angle = unit.angle,
                                messages.screen = TRUE, ...)
  class(res) <- "variog4"
  return(res)
  
}

"plot.variog4" <-
  function (obj, omnidirectional = FALSE, same.plot = TRUE, legend = TRUE,...)
{
     ymax <- max(obj[[1]]$v, obj[[2]]$v, obj[[3]]$v, obj[[4]]$v)
     n.o <- names(obj)[1:4]
     if (same.plot) {
         xx <- obj[[1]]$u
         yy <- cbind(obj[[1]]$v, obj[[2]]$v, obj[[3]]$v, obj[[4]]$v)
         if (omnidirectional)
             yy <- cbind(obj[[5]]$v, yy)
         GP <- list(...)
         if (is.null(GP$lty))
             GP$lty <- 1:5
         if (is.null(GP$lwd))
             GP$lwd <- 1
         if (is.null(GP$col))
             GP$col <- 1:5
         if (is.null(GP$pch))
             GP$pch <- NULL
         if (is.null(GP$type))
             GP$type <- "l"
         matplot(x = xx, y = yy, type = GP$type, ...)
         if (legend) {
             if (omnidirectional) {
                 legend(0, ymax, legend = c("omnid.", substitute(a *
                   degree, list(a = n.o[1])), substitute(a * degree,
                   list(a = n.o[2])), substitute(a * degree, list(a = n.o[3])),
                   substitute(a * degree, list(a = n.o[4])), expression()),
                   lty = GP$lty, lwd = GP$lwd, col = GP$col)
             }
             else {
                 legend(0, ymax, legend = c(substitute(a * degree,
                   list(a = n.o[1])), substitute(a * degree, list(a = n.o[2])),
                   substitute(a * degree, list(a = n.o[3])), substitute(a *
                     degree, list(a = n.o[4])), expression()),
                   lty = GP$lty, lwd = GP$lwd, col = GP$col)
             }
         }
     }
     else {
         temp.mf <- par()$mfrow
         par(mfrow = c(2, 2))
         GP <- list(...)
         if (is.null(GP$lty)) {
             GP$lty <- rep(1, 4)
             if (omnidirectional)
                 GP$lty <- c(GP$lty, 2)
         }
         else {
             if (length(GP$lty) == 1)
                 if (omnidirectional)
                   GP$lty <- rep(GP$lty, 5)
                 else GP$lty <- rep(GP$lty, 4)
             if (length(GP$lty) == 2)
                 if (omnidirectional)
                   GP$lty <- c(rep(GP$lty[1], 4), GP$lty[2])
                 else GP$lty <- c(rep(GP$lty, 4))
             if (length(GP$lty) == 4 & omnidirectional)
                 GP$lty <- c(rep(GP$lty, 2))
         }
         if (is.null(GP$lwd)) {
             GP$lwd <- rep(1, 4)
             if (omnidirectional)
                 GP$lwd <- c(GP$lwd, 1)
         }
         else {
             if (length(GP$lwd) == 1)
                 if (omnidirectional)
                   GP$lwd <- rep(GP$lwd, 5)
                 else GP$lwd <- rep(GP$lwd, 4)
             if (length(GP$lwd) == 2)
                 if (omnidirectional)
                   GP$lwd <- c(rep(GP$lwd[1], 4), GP$lwd[2])
                 else GP$lwd <- rep(GP$lwd, 4)
             if (length(GP$lwd) == 4 & omnidirectional)
                 GP$lwd <- c(rep(GP$lwd, 1))
         }
         if (is.null(GP$col)) {
             GP$col <- rep(1, 4)
             if (omnidirectional)
                 GP$col <- c(GP$col, 1)
         }
         else {
             if (length(GP$col) == 1)
                 if (omnidirectional)
                   GP$col <- rep(GP$col, 5)
                 else GP$col <- rep(GP$col, 4)
             if (length(GP$col) == 2)
                 if (omnidirectional)
                   GP$col <- c(rep(GP$col[1], 4), GP$col[2])
                 else GP$col <- rep(GP$col, 2)
             if (length(GP$col) == 4 & omnidirectional)
                 GP$col <- c(rep(GP$col, 1))
         }
         if (is.null(GP$pch)) {
             GP$pch <- rep(1, 4)
             if (omnidirectional)
                 GP$pch <- c(GP$pch, 1)
         }
         else {
             if (length(GP$pch) == 1)
                 if (omnidirectional)
                   GP$pch <- rep(GP$pch, 5)
                 else GP$pch <- rep(GP$pch, 4)
             if (length(GP$pch) == 2)
                 if (omnidirectional)
                   GP$pch <- c(rep(GP$pch[1], 4), GP$pch[2])
                 else GP$pch <- rep(GP$pch, 2)
             if (length(GP$pch) == 4 & omnidirectional)
                 GP$pch <- c(rep(GP$pch, 2))
         }
         if (is.null(GP$type)) {
             GP$type <- rep("l", 4)
             if (omnidirectional)
                 GP$type <- c(GP$type, "l")
         }
         else {
             if (length(GP$type) == 1)
                 if (omnidirectional)
                   GP$type <- rep(GP$type, 5)
                 else GP$type <- rep(GP$type, 4)
             if (length(GP$type) == 2 & omnidirectional)
                 GP$type <- c(rep(GP$type[1], 4), GP$type[2])
             if (length(GP$type) == 4 & omnidirectional)
                 GP$type <- c(rep(GP$type, 2))
         }
         for (i in 1:4) {
             plot(obj[[i]], ylim = c(0, ymax), type = GP$type[i],
                 col = GP$col[i], lwd = GP$lwd[i], lty = GP$lty[i],
                 pch = GP$pch[i])
             if (omnidirectional) {
                 lines(obj$omnidirectional, type = GP$type[5],
                   col = GP$col[5], lwd = GP$lwd[5], lty = GP$lty[5])
                 legend(0, ymax, legend = c(substitute(a * degree,
                   list(a = n.o[i])), "omn.", expression()),
		  			lty = c(GP$lty[i], GP$lty[5]),
					col=c(GP$col[i],  GP$col[5]),
					lwd=c(GP$lwd[i],  GP$lwd[5]),
					bty = "n")
             }
             else title(main = substitute(a * degree, list(a = n.o[i])),
                 cex = 1.3)
         }
         par(mfrow = temp.mf)
     }
     return(invisible())
}
