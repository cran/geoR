"grf" <-
  function(n, grid = "irreg", 
           nx, ny, xlims = c(0, 1), ylims = c(0, 1), nsim = 1, 
           cov.model = "matern",
           cov.pars = stop("covariance parameters (sigmasq and phi) needed"),
           kappa = 0.5,  nugget=0, lambda=1, aniso.pars = NULL,
           method = c("cholesky", "svd", "eigen", "circular.embedding"),
           messages)
{
  ##
  ## reading and checking input
  ##
  call.fc <- match.call()
  if(missing(messages))
    messages.screen <- ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages"))
  else messages.screen <- messages
  method <- match.arg(method)
  if((method == "circular.embedding") & messages.screen)
    cat("grf: for simulation of fields with large number of points the consider the package RandomFields.\n")
  ##
  ## defining the model to simulate from
  ##
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (cov.model == "matern" && kappa == 0.5) cov.model <- "exponential"
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
    nst <- 1
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
    nst <- nrow(cov.pars)
  }
  sill.total <- tausq + sum(sigmasq)
  messa <- grf.aux1(nst, nugget, sigmasq, phi, kappa, cov.model)
  if (messages.screen) {
    cat(messa$nst)
    cat(messa$nugget)
    cat(messa$cov.structures)
    cat(paste("grf: decomposition algorithm used is: ", method, "\n"))
  }
  ##
  ##
  ##
  rseed <- get(".Random.seed", envir=.GlobalEnv, inherits = FALSE)
  results <- list()
  ##
  ## defining the locations for the simulated data
  ##
  if(is.character(grid)) grid <- match.arg(grid, choices=c("irreg", "reg"))
  if (is.matrix(grid) | is.data.frame(grid)) {
    results$coords <- as.matrix(grid)
    if (messages.screen) 
      cat("grf: simulation on locations provided by the user\n")
  }
  else {
    ##
    ## checking whether it is a 1D simulation
    ##
    if((!missing(nx) && nx == 1) | (!missing(ny) && ny == 1) |
       diff(xlims) == 0 | diff(ylims) == 0){
      sim1d <- TRUE
      if (messages.screen) 
        cat("simulations in 1D\n")
    }
    else sim1d <- FALSE
    ##
    ## defining number of points in each direction
    ##
    if(missing(nx)){
      if(sim1d)
        if(diff(xlims) == 0) nx <- 1
        else nx <- n
      else
        if(is.character(grid) && grid == "reg") nx <- round(sqrt(n))
        else nx <- n
    }
    if(missing(ny)){
      if(sim1d)
        if(diff(ylims) == 0) ny <- 1
        else ny <- n
      else
        if(is.character(grid) && grid == "reg") ny <- round(sqrt(n))
        else ny <- n
    }
    ##
    ## defining the grid
    ##
    if (is.character(grid) && grid == "irreg") {
      results$coords <- cbind(x = runif(nx, xlims[1], xlims[2]),
                              y = runif(ny, ylims[1], ylims[2]))
      if (messages.screen) 
        cat(paste("grf: simulation(s) on randomly chosen locations with ", n, " points\n"))
    }
    else {
      xpts <- seq(xlims[1], xlims[2], l = nx)
      ypts <- seq(ylims[1], ylims[2], l = ny)
      results$coords <- as.matrix(expand.grid(x = xpts, y = ypts))
      if(length(xpts) == 1) xspacing <- 0
      else xspacing <- xpts[2] - xpts[1] 
      if(length(ypts) == 1) yspacing <- 0
      else yspacing <- ypts[2] - ypts[1] 
      if(abs(xspacing - yspacing) < 1e-12) equal.spacing <- TRUE
      else equal.spacing <- FALSE
      if (messages.screen) 
        cat(paste("grf: generating grid ", nx, " * ", ny, 
                  " with ", (nx*ny), " points\n"))
    }
  }
  n <- nrow(results$coords)
  if(length(unique(results$coords[,1])) == 1 |
     length(unique(results$coords[,2])) == 1)
    sim1d <- TRUE
  else sim1d <- FALSE
  ##
  ## transforming to the isotropic space 
  ##
  if(!is.null(aniso.pars)) {
    if(method == "circular.embedding")
      stop("anisotropic models not implemented for the circular embedding method. \nConsider using the package \"RandomFields")
    if(length(aniso.pars) != 2 | !is.numeric(aniso.pars))
      stop("anisotropy parameters must be provided as a numeric vector with two elements: the rotation angle (in radians) and the anisotropy ratio (a number greater than 1)")
    if(messages.screen)
      cat("grf: transforming to the isotropic space \n")
    results$coords <- coords.aniso(coords = results$coords,
                                   aniso.pars = aniso.pars)
  }
  ##
  ## simulating data at locations defined by the matrix results$coords
  ##
  if (all(phi) == 0) {
    results$data <- matrix(rnorm((n * nsim), mean = 0, sd = sqrt(sill.total)), 
                           nrow = n, ncol = nsim)
  }
  else {
    if (method == "circular.embedding") {
      if (is.character(grid) && grid == "irreg") 
        stop("Option for \"circular.embedding\" algorithm only allowed for regular grids. You might have to include the argument grid=\"reg\"")
      if(cov.model == "power")
        stop("power covariance model not implemented for the circular embedding method") 
      stepx <- (xlims[2] - xlims[1])/(nx - 1)
      stepy <- (ylims[2] - ylims[1])/(ny - 1)
      if (round(1e+08 * stepx) != round(1e+08 *stepy)) 
        stop("grf: distance between grid points must be the same in X and Y directions")
      temp <- list(n = n, nst = nst, sigmasq = sigmasq, 
                   xlims = xlims, ylims = ylims, stepx = stepx, 
                   cov.model = cov.model, phi = phi, kappa = kappa)
      if(messages.screen)
        cat("\ngrf: WARNING:\nmessages of the type mtot=XXXXX will appear on your screen. \nIf there are many (3 or more, say) or they run indefinitely, you should stop the simulation and try again with a different grid (e.g. try to add 1 point in each direction)\n")
      grf.aux3 <- function (nsim, temp) {
        realiz <- rep(0, temp$n)
        for (i in 1:temp$nst) {
          realiz <- realiz + sqrt(temp$sigmasq[i]) *
            grf.aux2(xlim = temp$xlims, 
                     ylim = temp$ylims, step = temp$stepx,
                     cov.model = temp$cov.model,
                     phi = temp$phi[i], kappa = temp$kappa)
          NULL
        }
        return(realiz)
      }      
      results$data <- apply(as.matrix(1:nsim), 1, grf.aux3, temp = temp)
      if (nugget != 0) {
        results$data <- results$data +
          matrix(rnorm((n * nsim), sd = sqrt(nugget)), ncol = nsim)
      }
    }
    else{
      results$data <- matrix(rnorm((n * nsim)), nrow = n, ncol = nsim)
      cov.decomp <- t(varcov.spatial(coords = results$coords, 
                                   cov.model = cov.model, kappa = kappa,
                                   nugget = nugget, cov.pars = cov.pars, 
                                   only.decomposition = TRUE,
                                   func.inv = method)$sqrt.varcov)
      results$data <- cov.decomp %*% results$data
    }
    if (nsim == 1) 
      results$data <- as.vector(results$data)
  }
  ##
  ## transforming data (Box - Cox)
  ##
  if (lambda != 1){
    if (lambda != 0)
      results$data <- (results$data * lambda + 1)^(1/lambda)
    else
      results$data <- exp(results$data)
    messa$transformation <- paste("grf: Data transformed (Box-Cox), for lambda =", lambda)
    if (messages.screen) 
      cat(messa$transformation); cat("\n")
  }
  ##
  ## back-transformation to the anisotropic space 
  ##
  if(!is.null(aniso.pars)) {
    if(messages.screen)
      cat("grf: back-transforming to the anisotropic space \n")
    results$coords <- coords.aniso(coords = results$coords,
                                   aniso.pars = aniso.pars, reverse=TRUE)
  }
  else{aniso.pars <- "no anisotropy parameters provided/used"}
  ##
  ## preparing output
  ##
  if (messages.screen) 
    cat(paste("grf: End of simulation procedure. Number of realizations:",
              nsim, "\n"))
  results  <- c(results, list(cov.model = cov.model, 
                              nugget = nugget, cov.pars = cov.pars,
                              kappa = kappa, lambda = lambda,
                              aniso.pars = aniso.pars, method = method,
#                              sim.dim = ifelse(sim1d, "1d", "2d"),
                              .Random.seed = rseed, messages = messa,
                              call = call.fc))
  if(is.character(grid) && grid == "reg"){
    if(equal.spacing) attr(results, "spacing") <- xspacing
    else{
      attr(results, "xspacing") <- xspacing
      attr(results, "yspacing") <- yspacing
    }
  }
  attr(results, 'sp.dim') <- ifelse(sim1d, "1d", "2d")
  class(results) <- c("grf", "geodata", "variomodel")
  return(results)
}


"grf.aux1" <-
  function (nst, nugget, sigmasq, phi, kappa, cov.model) 
{
  cov.nst <- paste("grf: process with ", nst, " covariance structure(s)\n")
  cov.nugget <- paste("grf: nugget effect is: tausq=", nugget,"\n")
  cov.message <- NULL
  for (i in 1:nst) {
    if (phi[i] == 0) 
      cov.message[i] <- paste("grf: covariance model", i, "is a pure nugget effect\n")
    else {
      if (cov.model == "matern" | cov.model == "powered.exponential" | 
          cov.model == "cauchy" | cov.model == "gneiting-matern") 
        cov.message[i] <- paste("grf: covariance model ", 
                                i, " is: ", cov.model, "(sigmasq=", sigmasq[i], 
                                ", phi=", phi[i], ", kappa = ", kappa, ")\n", sep = "")
      else cov.message[i] <- paste("grf: covariance model ", 
                                   i, " is: ", cov.model, "(sigmasq=", sigmasq[i], 
                                   ", phi=", phi[i], ")\n", sep = "")
    }
  }
  return(list(nst = cov.nst, nugget = cov.nugget, cov.structures = cov.message))
}

"grf.aux2" <-
  function (xlim, ylim, step, cov.model, phi, kappa = 0.5) 
{
  if(is.null(kappa)) kappa <- 1
  nx <- c(diff(xlim)/step + 1, diff(ylim)/step + 1)
  res <- double(nx[1] * nx[2])
  ln <- as.integer(2)
  mm <- integer(ln)
  x <- .C("woodandchan",
          as.integer(cor.number(cov.model)),
          as.integer(nx),
          ln,
          as.double(step),
          as.double(phi),
          as.double(kappa), 
          res = res,
          m = mm, PACKAGE = "geoR")$res
  cat("\n")
  return(x)
}

"lines.variomodel.grf" <-
  function (x, max.dist = max(dist(x$coords)), length = 100, 
            lwd = 2, ...) 
{
  if(! "package:stats" %in% search()) require(mva)
  if (x$cov.model == "matern" | x$cov.model == "powered.exponential" | 
      x$cov.model == "cauchy" | x$cov.model == "gneiting-matern") 
    kappa <- x$kappa
  else kappa <- NULL
  distance <- seq(0, max.dist, length = length)
  if (is.vector(x$cov.pars)) 
    sill.total <- x$nugget + x$cov.pars[1]
  else sill.total <- x$nugget + sum(x$cov.pars[, 1])
  gamma <- sill.total - cov.spatial(distance, cov.model = x$cov.model, 
                                  kappa = kappa, cov.pars = x$cov.pars)
  lines(distance, gamma, lwd = lwd, ...)
  return(invisible())
}

"plot.1d" <-
  function(x, xlim, ylim, x1vals, ...)
{
  #cat("data in 1-D\n")
  if(length(x1vals) == 1) col.ind <- 2
  else col.ind <- 1
  order.it <- order(x$coords[,col.ind])
  if(is.null(list(...)$xla)) xlabel <- "locations"
  else xlabel <- list(...)$xla
  if(is.null(list(...)$yla)) ylabel <- "data"
  else ylabel <- list(...)$yla
  pty.prev <- par()$pty
  par(pty="m")
  plot(x$coords[order.it,col.ind], x$data[order.it],
       xlab = xlabel, ylab = ylabel, xlim = xlim,
       ylim = ylim, ...)
  par(pty=pty.prev)
  return(invisible())
}

"image.grf" <-
  function (x, sim.number = 1, xlim, ylim,
            x.leg, y.leg, ...) 
{
  pty.prev <- par()$pty
  x1vals <- unique(round(x$coords[,1], dig=12))
  x2vals <- unique(round(x$coords[,2], dig=12))
  if(missing(xlim)) xlim <- NULL
  if(missing(ylim)) ylim <- NULL
  ##
  ## Plotting simulations in 1-D
  ##
  if(attr(x, 'sp.dim') == "1d" | length(x1vals) == 1 | length(x2vals) == 1)
    plot.1d(x, xlim=xlim, ylim = ylim, x1vals = x1vals, ...)
  else{
    ##
    ## Plotting simulations in 2-D
    ##
    ldots <- match.call(expand.dots = FALSE)$...
    ldots[[match(names(ldots), "offset.leg")]] <- NULL
    if(length(ldots[!is.na(match(names(ldots), "xlab"))])==0)
      ldots$xlab <- "X Coord"
    if(length(ldots[!is.na(match(names(ldots), "ylab"))])==0)
      ldots$ylab <- "Y Coord"
    ##
    ## Checking for retangular grid
    ##
    nx <- length(as.numeric(levels(as.factor(round(x$coords[, 1], dig=12)))))
    ny <- length(as.numeric(levels(as.factor(round(x$coords[, 2], dig=12)))))
    x$data <- as.matrix(x$data)
    n <- nrow(x$data)
    if (nx * ny != n) 
      stop("cannot produce image plot probably due to irregular grid of locations")
    ##
    ## Preparing image plot elements
    ##
    locations <- prepare.graph.kriging(locations=x$coords,
                                       values=x$data[, sim.number],
                                       borders =  NULL,
                                       borders.obj = eval(attr(x, "borders")),
                                       xlim = xlim, ylim = ylim) 
    ##
    par(pty = "s")
    do.call("image", c(list(x=locations$x, y=locations$y,
                            z=locations$values,
                            xlim = locations$coords.lims[,1],
                            ylim = locations$coords.lims[,2]),
                       ldots))
    ##
    ## Adding the legend (if the case)
    ##
    if(!missing(x.leg) && !missing(y.leg)){
      if(is.null(ldots$col)) ldots$col <- heat.colors(12)
      legend.krige(x.leg=x.leg, y.leg=y.leg,
                   values=locations$values[!is.na(locations$values)],
                   vertical = vertical, cex=cex.leg,
                   col=ldots$col, ...)
    }
  }
  par(pty = pty.prev)
  return(invisible())
}

#"image.grf" <-
#  function (x, sim.number = 1, ...) 
#{
#  x1vals <- unique(x$coords[,1])
#  x2vals <- unique(x$coords[,2])
#  if(x$sim.dim == "1d" | length(x1vals) == 1 | length(x2vals) == 1)
#    plot.1d(x, ...)
#  else{
#    xl <- as.numeric(levels(as.factor(round(x$coords[, 1], dig=12))))
#    nx <- length(xl)
 #   yl <- as.numeric(levels(as.factor(round(x$coords[, 2], dig=12))))
#    ny <- length(yl)
 #   x$data <- as.matrix(x$data)
 #   n <- nrow(x$data)
 #   if (nx * ny != n) 
 #     stop("cannot produce image plot probably due to data on irregular grid")
 ##   m <- matrix(x$data[, sim.number], ncol = ny)
 #   coords.lims <- set.coords.lims(coords=x$coords)
 #   x.ex <- diff(range(coords.lims[,1]))/(2*(nx-1))
 #   y.ex <- diff(range(coords.lims[,2]))/(2*(ny-1))
 #   xlim.ex <- coords.lims[,1] + c(-x.ex, x.ex)
 #   ylim.ex <- coords.lims[,2] + c(-y.ex, y.ex)
 #   pty.prev <- par()$pty
 #   par(pty = "s")
 #   image(xl, yl, m, xlim= xlim.ex, ylim=ylim.ex,...)
 #   par(pty=pty.prev)
 # }
 # return(invisible())
#}

"persp.grf" <- 
  function(x, sim.number = 1, ...)
{
  x1vals <- unique(round(x$coords[,1], dig=12))
  x2vals <- unique(round(x$coords[,2], dig=12))
  ldots <- list(...)
  if(is.null(ldots$xlim)) xlim <- NULL
  if(is.null(ldots$ylim)) ylim <- NULL
  if(attr(x, 'sp.dim') == "1d" | length(x1vals) == 1 | length(x2vals) == 1)
    plot.1d(x, xlim=xlim, ylim = ylim, x1vals = x1vals, ...)
  else{
    xl <- as.numeric(levels(as.factor(round(x$coords[, 1], dig=12))))
    nx <- length(xl)
    yl <- as.numeric(levels(as.factor(round(x$coords[, 2], dig=12))))
    ny <- length(yl)
    x$data <- as.matrix(x$data)
    n <- nrow(x$data)
    if(nx * ny != n)
      stop("cannot produce perspective plot, probably irregular grid")
    m <- matrix(x$data[, sim.number], ncol = ny)
    persp(xl, yl, m, ...)
  }
  return(invisible())
}

"plot.grf" <-
  function (x, model.line = TRUE, plot.locations = FALSE, ...) 
{
  nsim <- ncol(x$data)
  if (plot.locations){
    points.geodata(x, pt.divide="equal", xlab = "Coord X", ylab = "Coord Y")
    if(is.null(list(...)$ask)){
      ask.now <- par()$ask
      par(ask = TRUE)
      on.exit(par(ask=ask.now)) 
    }
  }
  if (is.vector(x$cov.pars)) 
    sill.total <- x$nugget + x$cov.pars[1]
  else sill.total <- x$nugget + sum(x$cov.pars[, 1])
  if (x$lambda != 1){
    if (x$lambda == 0) data <- log(x$data)
    else data <- ((x$data^x$lambda)-1)/x$lambda
  }
  else
    data <- x$data          
  sim.bin <- variog(x, data=data)
  plot(sim.bin, ...)
  if (model.line){
    var.model <- list(nugget = x$nugget, cov.pars = x$cov.pars, 
                      kappa = x$kappa, max.dist = max(sim.bin$u),
                      cov.model = x$cov.model)
    lines.variomodel(var.model, lwd = 3)
  }
  return(invisible())
}

"print.grf" <-
  function(x, ...)
{
  print.default(x, ...)
}

"lines.grf" <- function(x, ...){
  if(attr(x, "sp.dim") != "1d")
    stop("can only be used for simulations in  1-D")
  if(is.matrix(x$data))
    matplot(x$coords[,1], x$data, add=T, ...)
  else
    lines(x$coords[,1], x$data, ...)
  return(invisible())
}
