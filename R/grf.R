"grf" <-
  function(n, grid = "irreg",
           nx = round(sqrt(n)), ny = round(sqrt(n)),
           xlims = c(0, 1), ylims = c(0, 1), nsim = 1, 
           cov.model = "matern",
           cov.pars = stop("covariance parameters (sigmasq and phi) needed"),
           kappa = 0.5,  nugget=0, lambda=1, aniso.pars = NULL,
           method = c("cholesky", "svd", "eigen", "circular.embedding"),
           messages.screen = TRUE)
{
  ##
  ## reading and checking input
  ##
  call.fc <- match.call()
  method <- match.arg(method)
  if((method == "circular.embedding") & messages.screen)
    cat("grf: for simulation of fields with large number of points the consider the packages RandomFields.\n") 
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (cov.model == "matern" && kappa == 0.5)
    cov.model <- "exponential"
  rseed <- .Random.seed
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
  results <- list()
  ##
  ## defining the locations for the simulated data
  ##
  if (is.matrix(grid) | is.data.frame(grid)) {
    results$coords <- as.matrix(grid)
    if (messages.screen) 
      cat("grf: simulation(s) on a grid provided by the user\n")
  }
  else {
    if (grid == "irreg") {
      results$coords <- cbind(x = runif(n, xlims[1], xlims[2]),
                              y = runif(n, ylims[1], ylims[2]))
      if (messages.screen) 
        cat(paste("grf: simulation(s) on random locations with ", n, " points\n"))
    }
    else {
      xpts <- seq(xlims[1], xlims[2], l = nx)
      ypts <- seq(ylims[1], ylims[2], l = ny)
      results$coords <- as.matrix(expand.grid(x = xpts, y = ypts))
      xspacing <- xpts[2] - xpts[1] 
      yspacing <- ypts[2] - ypts[1]
      if((xspacing - yspacing) < 1e-12) equal.spacing <- TRUE
      if (messages.screen) 
        cat(paste("grf: generating grid ", nx, " * ", ny, 
                  " with ", (nx*ny), " points\n"))
    }
  }
  n <- nrow(results$coords)
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
      if (grid == "irreg") 
        stop("Option for \"circular.embedding\" algorithm only allowed for regular grids. You might have to include the argument grid=\"reg\"")
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
                              .Random.seed = rseed, messages = messa,
                              call = call.fc))
  if(grid == "reg"){
    if(equal.spacing) attr(results, "spacing") <- xspacing
    else{
      attr(results, "xspacing") <- xspacing
      attr(results, "yspacing") <- yspacing
    }
  }
  class(results) <- c("grf", "geodata")
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
          m = mm)$res
  cat("\n")
  return(x)
}

"lines.grf" <-
  function (x, max.dist = max(dist(x$coords)), length = 100, 
            lwd = 2, ...) 
{
  if(is.R()) require(mva)
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

"image.grf" <-
  function (x, sim.number = 1, ...) 
{
  xl <- as.numeric(levels(as.factor(x$coords[, 1])))
  nx <- length(xl)
  yl <- as.numeric(levels(as.factor(x$coords[, 2])))
  ny <- length(yl)
  x$data <- as.matrix(x$data)
  n <- nrow(x$data)
  if (nx * ny != n) 
    stop("Probably irregular grid")
  m <- matrix(x$data[, sim.number], ncol = ny)
  coords.lims <- set.coords.lims(coords=x$coords)
  pty.prev <- par()$pty
  par(pty = "s")
  image(xl, yl, m, xlim= coords.lims[,1], ylim=coords.lims[,2],...)
  par(pty=pty.prev)
  return(invisible())
}

"persp.grf" <- 
function(x, sim.number = 1, ...)
{
	x <- as.numeric(levels(as.factor(x$coords[, 1])))
	nx <- length(x)
	y <- as.numeric(levels(as.factor(x$coords[, 2])))
	ny <- length(y)
	x$data <- as.matrix(x$data)
	n <- nrow(x$data)
	if(nx * ny != n)
		stop("Probably irregular grid")
	m <- matrix(x$data[, sim.number], ncol = ny)
	persp(x, y, m, ...)
	return(invisible())
}

"plot.grf" <-
  function (x, model.line = TRUE, plot.grid = FALSE, ylim="default", ...) 
{
  nsim <- ncol(x$data)
  if (plot.grid) 
    points.geodata(x, pt.siz="equal", xlab = "Coord X", ylab = "Coord Y")
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


