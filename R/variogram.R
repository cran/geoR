"variog" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            uvec = "default", trend = "cte", lambda = 1,
            option = c("bin", "cloud", "smooth"),
            estimator.type = c("classical", "modulus"), 
            nugget.tolerance = 0, max.dist = NULL, pairs.min = 2,
            bin.cloud = FALSE, direction = "omnidirectional", tolerance = pi/8,
            unit.angle = c("radians","degrees"),
            messages.screen = TRUE, ...) 
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
    }
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
      uvec <- seq(0, umax, l = 12)
    if (is.numeric(uvec) & length(uvec) == 1) 
      uvec <- seq(0, umax, l = uvec)
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
      keep <- list(...)
      if(is.null(keep$keep.NA)) keep.NA <- FALSE
      else keep.NA <- keep$keep.NA
      result <- rfm.bin(cloud = list(u = u, v = v),
                        estimator.type = estimator.type, 
                        uvec = uvec, nugget.tolerance = nugget.tolerance, 
                        bin.cloud = bin.cloud, max.dist = umax, keep.NA = keep.NA)
      if(keep.NA){
        if (pairs.min > 0) {
          indp <- (result$n < pairs.min)
          if (is.matrix(result$v)) {
            result$v[indp, ] <- NA
            result$sd[indp, ] <- NA
          }
          else {
            result$v[indp] <- NA
            result$sd[indp] <- NA
          }
        }
        result$ind.bin <- indp
      }
      else{
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
  require(mva)
  u <- as.vector(dist(as.matrix(coords)))
  if(length(direction) != 4)
    stop("argument direction must be a vector with 4 values. For different specifications use the functio variog()")
  if(length(tolerance) != 1)
    stop("only 1 values can be provided to the argument tolerance . For different specifications use the functio variog()")
  res <- list()
  if(unit.angle == "radians")
    dg <- direction * 180/pi
  else dg <- direction
  if (!is.null(max.dist)) 
    umax <- max(u[u < max.dist])
  else umax <- max(u)
  if (all(uvec == "default")) 
    uvec <- seq(0, umax, l = 12)
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
  u <- NULL
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
             messages.screen = TRUE, keep.NA = TRUE)
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
                                messages.screen = TRUE,
                                keep.NA = TRUE 
                                )
  class(res) <- "variog4"
  return(res)
  
}

"plot.variog4" <-
  function (obj, omnidirectional = FALSE, same.plot = TRUE, legend = TRUE,...)
{
  ymax <- max(c(obj[[1]]$v, obj[[2]]$v, obj[[3]]$v, obj[[4]]$v), na.rm=T)
  n.o <- names(obj)[1:4]
  GP <- list(...)
  if(is.null(GP$xlab)) GP$xlab <- "distance"
  if(is.null(GP$ylab)) GP$ylab<- "semi-variance"
  if (same.plot) {
    xx <- obj[[5]]$u
    yy <- cbind(obj[[1]]$v, obj[[2]]$v, obj[[3]]$v, obj[[4]]$v)
    if (omnidirectional)
      yy <- cbind(obj[[5]]$v, yy)
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
    matplot(x = xx, y = yy, type = GP$type, lty=GP$lty, lwd=GP$lwd, col=GP$col, pch=GP$pch, xlab=GP$xlab, ylab=GP$ylab, ylim=c(0,max(yy)))
    if (legend) {
      if (omnidirectional) {
        legend(0, ymax,
               legend = c("omnid.",
                 substitute(a * degree, list(a = n.o[1])),
                 substitute(a * degree, list(a = n.o[2])),
                 substitute(a * degree, list(a = n.o[3])),
                 substitute(a * degree, list(a = n.o[4])),
                 expression()),
               lty = GP$lty, lwd = GP$lwd, col = GP$col)
      }
      else {
        legend(0, ymax,
               legend = c(substitute(a * degree,
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
      plot.default(obj[[i]]$u, obj[[i]]$v,
                   ylim = c(0, ymax), type = GP$type[i],
           col = GP$col[i], lwd = GP$lwd[i], lty = GP$lty[i],
           pch = GP$pch[i], xlab=GP$xlab, ylab=GP$ylab)
      if (omnidirectional) {
        lines(obj$omnidirectional, type = GP$type[5],
              col = GP$col[5], lwd = GP$lwd[5], lty = GP$lty[5])
        legend(0, ymax, legend = c(substitute(a * degree,
                          list(a = n.o[i])), "omn.", expression()),
               lty = c(GP$lty[i], GP$lty[5]),
               col=c(GP$col[i],  GP$col[5]),
               lwd=c(GP$lwd[i],  GP$lwd[5]))
      }
      else title(main = substitute(a * degree, list(a = n.o[i])),
                 cex = 1.3)
    }
    par(mfrow = temp.mf)
  }
  return(invisible())
}


"rfm.bin" <-
  function (cloud, l = 15, uvec = "default", nugget.tolerance = 0, 
            estimator.type = c("classical", "robust"), bin.cloud = FALSE,
            max.dist, keep.NA = FALSE)
{
  if (all(uvec == "default")) 
    uvec <- seq(0, max(cloud$u), l = l)
  estimator.type <- match.arg(estimator.type)
  ##  if(nugget.tolerance > 0) {
  ##    dnug <- mean(cloud$u[cloud$u <= nugget.tolerance])
  ##    cloud$u[cloud$u <= nugget.tolerance] <- 0
  ##    uvec <- uvec[uvec > nugget.tolerance]
  ##  }
  ##  u <- c(0, uvec)
  ##  n <- length(u)
  if(all(uvec == "default"))
    uvec <- seq(0, max.dist, l = 15)
  ubin <- c(0, uvec)
  nvec <- length(ubin)
  d <- 0.5 * diff(ubin[2:nvec])
  bins.lim <- c(0, (ubin[2:(nvec - 1)] + d), (d[nvec - 2] + ubin[
                                                                 nvec]))
  if(uvec[1] == 0 & nugget.tolerance == 0)
    uvec[1] <- (bins.lim[1] + bins.lim[2])/2
  if(nugget.tolerance > 0) {
    bins.lim <- c(0, nugget.tolerance, bins.lim[bins.lim >
                                                nugget.tolerance])
    uvec <- c(0, (bins.lim[ - (1:2)] - 0.5 * diff(bins.lim)[
                                                            -1]))
  }
  nbins <- nc <- length(bins.lim) - 1
  if(is.null(max.dist))
    max.dist <- max(bins.lim)
  min.dist <- min(cloud$u)
  ##	d <- 0.5 * (u[3:n] - u[2:(n - 1)])
  ##	low <- c(0, (u[2:(n - 1)] + d))
  ##	high <- c((u[3:n] - d), (d[n - 2] + u[n]))
  ##	nc <- n - 1
  if (!is.matrix(cloud$v)) {
    vbin <- rep(0, nc)
    nbin <- rep(0, nc)
    sdbin <- rep(0, nc)
    if (bin.cloud == TRUE) 
      bins.clouds <- list()
    for (i in 1:nc) {
      ind <- (cloud$u > bins.lim[i]) & (cloud$u <= bins.lim[i+1])
      vbin[i] <- mean(cloud$v[ind])
      if (bin.cloud == TRUE) 
        bins.clouds[[i]] <- cloud$v[ind]
      nbin[i] <- sum(ind)
      if (estimator.type == "robust") 
        vbin[i] <- ((vbin[i])^4)/(0.914 + (0.988/nbin[i]))
      if (nbin[i] > 0) 
        sdbin[i] <- sqrt(var(cloud$v[ind]))
      else sdbin[i] <- NA
      NULL
    }
    if (uvec[1] == 0) 
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (min.dist == 0) {
      ind <- (cloud$u == 0)
      n.zero <- sum(ind)
      v.zero <- mean(cloud$v[ind])
      if (bin.cloud == TRUE) {
        bins.clouds[2:(length(bins.clouds) + 1)] <- bins.clouds[1:nc]
        bins.clouds[[1]] <- cloud$v[ind]
      }
      if (estimator.type == "robust") 
        v.zero <- ((v.zero)^4)/(0.914 + (0.988/n.zero))
      if (n.zero > 0) 
        sd.zero <- sqrt(var(cloud$v[ind]))
      else sd.zero <- NA
      uvec <- c(0, uvec)
      vbin <- c(v.zero, vbin)
      nbin <- c(n.zero, nbin)
      sdbin <- c(sd.zero, sdbin)
    }
    if(keep.NA == FALSE){
      u <- uvec[!is.na(vbin)]
      v <- vbin[!is.na(vbin)]
      n <- nbin[!is.na(vbin)]
      sd <- sdbin[!is.na(vbin)]
    }
    else{
      u <- uvec
      v <- vbin
      n <- nbin
      sd <- sdbin
    }
    if (bin.cloud == TRUE) 
      bins.clouds <- bins.clouds[!is.na(vbin)]
  }
  else {
    if (bin.cloud == TRUE) 
      stop("option bins.cloud=T allowed only for 1 variable")
    nvcols <- ncol(cloud$v)
    vbin <- matrix(0, nrow = nc, ncol = nvcols)
    nbin <- rep(0, nc)
    sdbin <- matrix(0, nrow = nc, ncol = nvcols)
    for (i in 1:nc) {
      ind <- (cloud$u >= bins.lim[i]) & (cloud$u < bins.lim[i+1])
      nbin[i] <- sum(ind)
      for (j in 1:nvcols) {
        vbin[i, j] <- mean(cloud$v[ind, j])
        if (estimator.type == "robust") 
          vbin[i, j] <- ((vbin[i, j])^4)/(0.914 + (0.988/nbin[i]))
        if (nbin[i] > 0) 
          sdbin[i, j] <- sqrt(var(cloud$v[ind, j]))
        else sdbin[i, j] <- NA
      }
      NULL
    }
    if (uvec[1] == 0) 
      uvec[1] <- (bins.lim[1] + bins.lim[2])/2
    if (min.dist == 0) {
      v.zero <- rep(0, nvcols)
      n.zero <- rep(0, nvcols)
      sd.zero <- rep(0, nvcols)
      for (j in 1:nvcols) {
        ind <- (cloud$u == 0)
        n.zero[j] <- sum(ind)
        v.zero[j] <- mean(cloud$v[ind, j])
        if (estimator.type == "robust") 
          v.zero[j] <- ((v.zero[j])^4)/(0.914 + (0.988/n.zero[j]))
        if (n.zero[j] > 0) 
          sd.zero[j] <- sqrt(var(cloud$v[ind, j]))
        else sd.zero[j] <- NA
        uvec <- c(0, uvec)
        vbin <- rbind(v.zero, vbin)
        nbin <- c(n.zero, nbin)
        sdbin <- rbind(sd.zero, sdbin)
      }
    }
    if(keep.NA == FALSE){
      u <- uvec[!is.na(vbin[, 1])]
      n <- nbin[!is.na(vbin[, 1])]
      v <- matrix(0, nrow = length(u), ncol = nvcols)
      sd <- matrix(0, nrow = length(u), ncol = nvcols)
    }
    else{
      u <- uvec
      n <- nbin
      v <- matrix(0, nrow = length(u), ncol = nvcols)
      sd <- matrix(0, nrow = length(u), ncol = nvcols)
    }
    for (j in 1:nvcols) {
      if(keep.NA == FALSE){
        v[, j] <- vbin[!is.na(vbin[, j]), j]
        sd[, j] <- sdbin[!is.na(vbin[, j]), j]
      }
      else{
        v[, j] <- vbin[, j]
        sd[, j] <- sdbin[, j]
      }
    }
  }
  if (nugget.tolerance > 0) {
    u[1] <- nugget.tolerance
  }
  result <- list(u = u, v = v, n = n, sd = sd, output.type = "bin", bins.lim = bins.lim)
  if (!is.matrix(cloud$v) && bin.cloud == TRUE) 
    result$bin.cloud <- bins.clouds
  if (!is.null(class(cloud))) 
    class(result) <- class(cloud)
  return(result)
}

"plot.variogram" <-
  function (obj, max.dist, vario.col = "all", scaled = FALSE,  
            var.lines = FALSE,  envelope.obj = NULL,
            bin.cloud = FALSE,  ...) 
{
  if(missing(max.dist)) max.dist <- max(obj$u)
  Ldots <- list(...)
  if(is.null(Ldots$xlab)) Ldots$xlab <- "distance"
  if(is.null(Ldots$ylab)) Ldots$ylab <- "semi-variance"
  if(is.null(Ldots$ty)){
    if (obj$output.type == "bin") Ldots$type <- "b"
    if (obj$output.type == "smooth") Ldots$type <- "l"
    if (obj$output.type == "cloud") Ldots$type <- "p"
  }
 if (bin.cloud == TRUE &&  Ldots$type != "b") 
    stop("plot.variogram: object must be a binned variogram with option bin.cloud=T")
  if (bin.cloud == TRUE && all(is.na(obj$bin.cloud))) 
    stop("plot.variogram: object must be a binned variogram with option bin.cloud=T")
  if (bin.cloud == TRUE && any(!is.na(obj$bin.cloud))) 
    boxplot(obj$bin.cloud, varwidth = TRUE, 
            xlab = "midpoints of distance class",
            ylab = paste("variogram values / ", 
              obj$estimator.type, "estimator"))
  else {
    u <- obj$u[obj$u <= max.dist]
    v <- obj$v
    if(is.vector(v) | length(v) == length(obj$u))
      v <- matrix(v, ncol=1)
    v <- v[obj$u <= max.dist,, drop=FALSE]
    if(vario.col == "all")
      vario.col <- 1:dim(v)[2]
    else
      if(!is.numeric(vario.col) | any(vario.col > ncol(v)))
        stop("argument vario.col must be equals to \"all\" or a vector indicating the column numbers to be plotted")
    v <- v[, vario.col, drop=F]
    if (scaled)
      v <- t(t(v)/obj$var.mark[vario.col])
    if (is.null(Ldots$ylim)){
      ymax <- max(v)
      if (!is.null(envelope.obj)) 
        ymax <- max(c(envelope.obj$v.upper, ymax))
      Ldots$ylim <- c(0, ymax)
    }
    if(ncol(v) == 1){
      v <- as.vector(v)
      plot(x= u, y= v, xlim = c(0, max.dist), ylim = Ldots$ylim, 
           xlab = Ldots$xlab, ylab = Ldots$ylab, type = Ldots$type)
    }
    else
      matplot(x=u, y= v, xlim = c(0, max.dist), ylim = Ldots$ylim, 
           xlab = Ldots$xlab, ylab = Ldots$ylab, type = Ldots$type)
    if (var.lines) {
      if (scaled) abline(h = 1, lty = 3)
      else abline(h = obj$var.mark, lty = 3)
    }
    if (!is.null(envelope.obj)) {
      lines(u, envelope.obj$v.lower, lty = 4)
      lines(u, envelope.obj$v.upper, lty = 4)
    }
  }
  return(invisible())
}

"lines.variogram" <-
function (obj, max.dist, type = "o", scaled = FALSE, ...) 
{
  if(missing(max.dist)) max.dist <- max(obj$u)
  if (scaled) 
    obj$v <- obj$v/obj$var.mark
  if (!is.matrix(obj$v)) 
    lines(obj$u[obj$u <= max.dist], obj$v[obj$u <= max.dist], 
          type = type, ...)
  else {
    for (j in 1:ncol(obj$v)) lines(obj$u[obj$u <= max.dist], 
                                   obj$v[obj$u <= max.dist, j], type = type, ...)
  }
}
