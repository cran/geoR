"lines.grf" <-
  function (obj, max.dist = max(dist(obj$coords)), length = 100, 
            lwd = 2, ...) 
{
  if(is.R()) require(mva)
  if (obj$cov.model == "matern" | obj$cov.model == "powered.exponential" | 
      obj$cov.model == "cauchy" | obj$cov.model == "gneiting-matern") 
    kappa <- obj$kappa
  else kappa <- NULL
  distance <- seq(0, max.dist, length = length)
  if (is.vector(obj$cov.pars)) 
    sill.total <- obj$nugget + obj$cov.pars[1]
  else sill.total <- obj$nugget + sum(obj$cov.pars[, 1])
  gamma <- sill.total - cov.spatial(distance, cov.model = obj$cov.model, 
                                  kappa = kappa, cov.pars = obj$cov.pars)
  lines(distance, gamma, lwd = lwd, ...)
  return(invisible())
}

"cov.spatial" <-
  function(obj, cov.model = c("matern", "exponential", "gaussian", "spherical",
                  "circular", "cubic", "wave", "power", "powered.exponential",
                  "cauchy", "gneiting", "gneiting.matern", "pure.nugget"),
           cov.pars = stop("no cov.pars argument provided"), kappa = 0.5)
{
  ##
  ## checking/reading input
  ##
  cov.model <- match.arg(cov.model)
  if(cov.model == "matern" | cov.model == "powered.exponential" | 
     cov.model == "cauchy" | cov.model == "gneiting.matern"){
    if(is.null(kappa))
      stop("for matern, powered.exponential, cauchy and gneiting.matern covariance functions the parameter kappa must be provided")
    if(cov.model == "gneiting.matern" & length(kappa) != 2)
      stop("gneiting.matern correlation function model requires a vector with 2 parameters in the argument kappa")
    if((cov.model == "matern" | cov.model == "powered.exponential" | 
        cov.model == "cauchy") & length(kappa) != 1)
      stop("for this choice of  correlation function model kappa should be a scalar parameter")
    if(cov.model == "matern" & kappa == 0.5)
      cov.model == "exponential"
  }
  if(is.vector(cov.pars))
    sigmasq <- cov.pars[1]
  else sigmasq <- cov.pars[, 1]
  if(is.vector(cov.pars))
    phi <- cov.pars[2]
  else phi <- cov.pars[, 2]
  if(is.vector(cov.pars))
    ns <- 1
  else ns <- nrow(cov.pars)
  covs <- array(0, dim = dim(obj))
  ##
  ## computing correlations/covariances
  ##
  for(i in 1:ns) {
    if(phi[i] < 1e-12)
      cov.model <- "pure.nugget"
    cov.values <- switch(cov.model,
                  pure.nugget = rep(0, length(obj)),
                  wave = (1/obj) * (phi[i] * sin(obj/phi[i])),
                  exponential = exp( - (obj/phi[i])),
                  matern = matern(u = obj, phi = phi[i], kappa = kappa),
                  gaussian = exp( - ((obj/phi[i])^2)),
                  spherical = ifelse(obj < phi[i], (1 - 1.5 * (obj/phi[i]) +
                    0.5 * (obj/phi[i])^3), 0),
                  circular = {
                    obj.sc <- obj/phi[i];
                    obj.sc[obj.sc > 1] <- 1;
                    ifelse(obj < phi[i], (1 - (2 * ((obj.sc) *
                                                    sqrt(1 - ((obj.sc)^2)) +
                                                    asin(obj.sc)))/pi), 0)
                  },
                  cubic = {
                    obj.sc <- obj/phi[i];
                    ifelse(obj < phi[i], (1 - (7 * (obj.sc^2) -
                                               8.75 * (obj.sc^3) +
                                               3.5 * (obj.sc^5) -
                                               0.75 * (obj.sc^7))), 0)
                  },
                  power = (obj)^phi,
                  powered.exponential = {
                    if(kappa > 2 | kappa <= 0)
                      stop("for power exponential correlation model the parameter kappa must be in the intervel (0,2]"
                           );
                    exp( - ((obj/phi[i])^kappa))
                  },
                  cauchy = (1 + (obj/phi[i])^2)^(-1 * kappa),
                  gneiting = {
                    obj.sc <- obj/phi[i];
                    t2 <- (1 - obj.sc);
                    t2 <- ifelse(t2 > 0, (t2^8), 0);
                    (1 + 8 * obj.sc + 25 * (obj.sc^2) + 32 * (obj.sc^
                                                              3)) * t2
                  },
                  gneiting.matern = { 
                    obj.sc <- obj/(phi[i] * kappa[2]);
                    t2 <- (1 - obj.sc);
                    t2 <- ifelse(t2 > 0, (t2^8), 0);
                    cov.values <- (1 + 8 * obj.sc + 25 * (obj.sc^2) + 32 * (obj.sc^3)) * t2;
                    cov.values * matern(u = obj, phi = phi[i], kappa = kappa[1])
                  },
                  stop("wrong or no specification of cov.model")
                  )
    cov.values <- sigmasq[i] * cov.values
    covs <- covs + cov.values
  }
  covs[obj < 1e-15] <- sum(sigmasq)
  return(covs)
}

"cor.number" <- 
  function(cov.model= c("exponential", "matern", "gaussian",
             "spherical", "circular", "cubic", "wave", "power",
             "powered.exponential", "cauchy", "gneiting",
             "gneiting.matern", "pure.nugget"))
{
###	WARNING: codes above must be the same as in the C code
###              "cor_diag"
  
  cov.model <- match.arg(cov.model)
  cornumber <- switch(cov.model,
                      pure.nugget = as.integer(1),
                      exponential = as.integer(2),
                      spherical = as.integer(3),
                      gaussian = as.integer(4),
                      wave = as.integer(5),
                      cubic = as.integer(6),
                      power = as.integer(7),
                      powered.exponential = as.integer(8),
                      cauchy = as.integer(9),
                      gneiting = as.integer(10),
                      circular = as.integer(11),
                      matern = as.integer(12),
                      gneiting.matern = as.integer(13),
                      stop("wrong or no specification of cov.model")
                      )
  return(cornumber)
}

"read.geodata" <-
  function(file, header = FALSE, coords.col= 1:2, data.col = 3,
           data.names = NULL, covar.col = NULL, covar.names = "header", ...)
{
  call.fc <- match.call()
  obj <- read.table(file = file, header = header, ...)
  if(covar.names == "header"){
    if(!is.null(covar.col)){
      col.names <- names(obj)
      covar.names <- col.names[covar.col]
    }
    else covar.names <- NULL
  }
  res <- as.geodata(obj = obj, coords.col = coords.col, data.col = data.col,
                    covar.col = covar.col, covar.names = covar.names)
  res$call <- call.fc
  return(res)
}

"as.geodata" <-
  function(obj, coords.col = 1:2, data.col = 3, data.names = NULL, 
           covar.col = NULL, covar.names = "obj.names")
{
  if(!is.matrix(obj) & !is.data.frame(obj))
    stop("object must be a matrix or data.frame")
  if(!is.null(data.names) & length(data.col) < 2)
    stop("data.names allowed only if there is more than 1 column of data")
  res <- list()
  res$coords <- as.matrix(obj[,coords.col])
  res$data <- as.matrix(obj[,data.col])
  if(length(data.col) == 1) res$data <- as.vector(res$data)
  else{
    res$data <- as.data.frame(res$data)
    if(!is.null(data.names)) names(res$data) <- data.names
  }
  if(!is.null(covar.col)){
    res[[3]] <- as.data.frame(as.matrix(obj[,covar.col]))
    if(covar.names == "obj.names"){
      if(is.matrix(obj))      
        col.names <- dimnames(obj)[2]
      if(is.data.frame(obj))      
        col.names <- names(obj)
    }
    if(length(covar.col) == 1){
      if(covar.names == "obj.names")
        names(res)[3] <- col.names[covar.col]
      else
        names(res)[3] <- covar.names
    }
    else{
      names(res)[3] <- "covariates"
      if(covar.names == "obj.names")
        names(res[[3]]) <- col.names[covar.col]
      else
        names(res[[3]]) <- covar.names
    }
  }
  require(mva)
  if(min(dist(res$coords)) < 1e-12)
    cat("WARNING: there are data at coincident locations, several geoR functions will not work\n") 
  class(res) <- "geodata"
  return(res)
}


"image.grf" <-
  function (obj, sim.number = 1, ...) 
{
  x <- as.numeric(levels(as.factor(obj$coords[, 1])))
  nx <- length(x)
  y <- as.numeric(levels(as.factor(obj$coords[, 2])))
  ny <- length(y)
  obj$data <- as.matrix(obj$data)
  n <- nrow(obj$data)
  if (nx * ny != n) 
    stop("Probably irregular grid")
  m <- matrix(obj$data[, sim.number], ncol = ny)
  coords.lims <- apply(obj$coords,2,range)
  coords.diff <- diff(coords.lims)
  if(coords.diff[1] != coords.diff[2]){
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[,ind.min] <- coords.lims[,ind.min] + c(-coords.diff.diff, coords.diff.diff)/2
  }
  par(pty = "s")
  image(x, y, m, xlim= coords.lims[,1], ylim=coords.lims[,2],...)
  return(invisible())
}

"image.kriging" <-
  function (kriging.obj, locations = kriging.obj$locations, values = kriging.obj$predict, coords.data, ...) 
{
##  par.ori <- par(no.readonly = TRUE)
##  on.exit(par(par.ori))
  if(is.null(locations))
    stop("prediction locations must be provided")
  if(is.null(values))
    stop("data.to be plotted must be provided")
  locations <- locations[order(locations[, 2], locations[,1]), ]
  x <- as.numeric(levels(as.factor(locations[, 1])))
  nx <- length(x)
  y <- as.numeric(levels(as.factor(locations[, 2])))
  ny <- length(y)
  if  (nx == ny) 
    par(pty = "s")
  values <- matrix(values, ncol = ny)
  coords.lims <- apply(locations,2,range)
  coords.diff <- diff(coords.lims)
  if(coords.diff[1] != coords.diff[2]){
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[,ind.min] <- coords.lims[,ind.min] + c(-coords.diff.diff, coords.diff.diff)/2
  }
  par(pty = "s")
  image(x, y, values, xlim= coords.lims[,1], ylim=coords.lims[,2],...)
  if(!missing(coords.data))
    points(coords.data)
  return(invisible())
}

"lines.variogram" <-
function (obj, max.dist = max(obj$u), type = "o", scaled = FALSE, ...) 
{
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
"lines.variomodel" <-
function (obj, max.dist = obj$max.dist, n.points = 100, scaled = FALSE,...)
{
  if (is.null(max.dist)) 
    stop("argument max.dist needed for this object")
  if (obj$cov.model == "matern" | obj$cov.model == "powered.exponential" | 
      obj$cov.model == "cauchy" | obj$cov.model == "gneiting-matern") 
    kappa <- obj$kappa
  else kappa <- NULL
  distance <- seq(0, max.dist, length = n.points)
  if (is.vector(obj$cov.pars)) 
    sill.total <- obj$nugget + obj$cov.pars[1]
  else sill.total <- obj$nugget + sum(obj$cov.pars[, 1])
  if (scaled){
    if(is.vector(obj$cov.model)) obj$cov.model[1] <-  obj$cov.model[1]/sill.total
    else obj$cov.model[,1] <-  obj$cov.model[,1]/sill.total
    sill.total <- 1
  }
  if (obj$cov.pars[2] < 1e-12)
    gamma <- rep(sill.total, n.points)
  else
    gamma <- sill.total - cov.spatial(distance, cov.model = obj$cov.model, 
                                    kappa = kappa, cov.pars = obj$cov.pars)
  lines(distance, gamma, ...)
}

"persp.kriging" <- 
  function(kriging.obj, locations = kriging.obj$locations, values = kriging.obj$
           predict, ...)
{
  locations <- locations[order(locations[, 2], locations[, 1]),  ]
  x <- as.numeric(levels(as.factor(locations[, 1])))
  nx <- length(x)
  y <- as.numeric(levels(as.factor(locations[, 2])))
  ny <- length(y)
  if(nx == ny)
    par(pty = "s")
  values <- matrix(values, ncol = ny)
  persp(x, y, values, ...)
  return(invisible())
}

"persp.grf" <- 
function(obj, sim.number = 1, ...)
{
	x <- as.numeric(levels(as.factor(obj$coords[, 1])))
	nx <- length(x)
	y <- as.numeric(levels(as.factor(obj$coords[, 2])))
	ny <- length(y)
	obj$data <- as.matrix(obj$data)
	n <- nrow(obj$data)
	if(nx * ny != n)
		stop("Probably irregular grid")
	m <- matrix(obj$data[, sim.number], ncol = ny)
	persp(x, y, m, ...)
	return(invisible())
}

"coords.aniso" <- 
  function(coords, aniso.pars, reverse=FALSE)
{
  coords <- as.matrix(coords)
  n <- nrow(coords)
  if(length(aniso.pars) != 2)
    stop("argument aniso.pars must be a vector with 2 elementsm the anisotropy angle and anisotropy ratio, respectively")
  psiA <- aniso.pars[1]
  psiR <- aniso.pars[2]
  if(psiR < 1){
    psiR <- round(psiR, dig=8)
    if(psiR < 1)
      stop("anisotropy ratio must be greater than 1")
  }
  rm <- matrix(c(cos(psiA), -sin(psiA),
                 sin(psiA), cos(psiA)),
               ncol = 2)
  tm <- diag(c(1, 1/psiR))
  if(reverse)
    coords.mod <- coords %*% solve(rm %*% tm)
  else
    coords.mod <- coords %*% rm %*% tm
  return(coords.mod)
}


"plot.grf" <-
  function (obj, model.line = TRUE, plot.grid = FALSE, ylim="default", ...) 
{
  nsim <- ncol(obj$data)
  if (plot.grid) 
    plot(obj$coords[, 1], obj$coords[, 2], xlab = "Coord X", ylab = "Coord Y")
  if (is.vector(obj$cov.pars)) 
    sill.total <- obj$nugget + obj$cov.pars[1]
  else sill.total <- obj$nugget + sum(obj$cov.pars[, 1])
  if (obj$lambda != 1){
    if (obj$lambda == 0) data <- log(obj$data)
    else data <- ((obj$data^obj$lambda)-1)/obj$lambda
  }
  else
    data <- obj$data          
  sim.bin <- variog(obj)
  plot(sim.bin, var.lines = FALSE, ylim = c(0,max(sim.bin$v, sill.total)), ...)
  if (model.line){
    var.model <- list(nugget = obj$nugget, cov.pars = obj$cov.pars, 
                      kappa = obj$kappa, max.dist = max(sim.bin$u),
                      cov.model = obj$cov.model)
    lines.variomodel(var.model, lwd = 3)
  }
  return(invisible())
}


"trend.spatial" <-
  function(trend, coords=NULL)
{
  if(inherits(trend, "formula")){
    trend.mat <- model.matrix(trend)
  }
  else{
    if(trend == "cte")
      trend.mat <- as.matrix(rep(1, nrow(coords)))
    else
      if(trend == "1st")
        trend.mat <- cbind(1, coords)
      else
        if(trend == "2nd")
          trend.mat <- cbind(1, coords, coords[, 1]^2, coords[
                                                              , 2]^2, coords[, 1] * coords[, 2])
	else
          stop("external trend must be provided for data locations to be estimated using the argments `trend.d` and `trend.l`. Both (trend.d and trend.l) must be the same `cte`, `1st`, `2nd` or given by a formula of the type ~X where X is a matrix or vector of covariates."
               )
  }
  trend.mat <- as.matrix(trend.mat)
  dimnames(trend.mat) <- list(NULL, NULL)
  return(trend.mat)
}

#"dist0.krige" <-
#function (x0, coords) 
#{
#  if (length(x0) != 2) 
#    stop(paste("length of x0 is", length(x0), "(it must be 2)"))
#  coords[, 1] <- coords[, 1] - x0[1]
#  coords[, 2] <- coords[, 2] - x0[2]
#  return(sqrt(coords[, 1]^2 + coords[, 2]^2))
#}


"polygrid" <- 
  function(xgrid, ygrid, poly, vec.inout = F)
{
  if(is.R())
  	require(splancs)
  if(exists("inout")){
    xygrid <- expand.grid(x = xgrid, y = ygrid)
    ind <- as.vector(inout(xygrid, poly))
    xypoly <- xygrid[ind == T,  ]
    if(vec.inout == F)
      return(xypoly)
    else return(list(xypoly = xypoly, vec.inout = ind))
  }
  else{
    cat("ERROR: cannot run the function\n")
    cat("package \"splancs\" should be installed/loaded")
    return(invisible())
  }
}

#"variog.env" <-
#  function (obj.variog, coords, model.pars, nsim = 99, messages.screen = TRUE)  
#{
#  cat("This function has been made obsolete\n")
#  cat("There are now two functions for variogram envelops:\n")
#  cat(" - variog.env.model:\n")
#  cat("       the same as the previous variog.env, based on the model")
#  cat(" - variog.env.mc:\n")
#  cat("       the new one based on permutations of the data")
#  return(invisible())
#}

"variog.model.env" <-
  function(geodata, coords = geodata$coords, obj.variog,
           model.pars, nsim = 99, save.sim = FALSE,
           messages.screen = TRUE) 
{
  call.fc <- match.call()
  obj.variog$v <- NULL
  ##
  ## reading input
  ##
  if(!is.null(model.pars$beta)) beta <- model.pars$beta
  else beta <- 0
  if(!is.null(model.pars$cov.model))
    cov.model <- model.pars$cov.model
  else cov.model <- "exponential"
  if(!is.null(model.pars$kappa)) kappa <- model.pars$kappa
  else kappa <- 0.5
  if(!is.null(model.pars$nugget)) nugget <- model.pars$nugget
  else nugget <- 0
  cov.pars <- model.pars$cov.pars
  if(!is.null(obj.variog$estimator.type))
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (obj.variog$output.type != "bin") 
    stop("envelops can be computed only for binned variogram")
  ##
  ## generating simulations from the model with parameters provided
  ##
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations (with ",
              obj.variog$n.data, 
              "points each) using the function grf\n"))
  simula <- grf(obj.variog$n.data, grid = as.matrix(coords),
                cov.model = cov.model, cov.pars = cov.pars,
                nugget = nugget, kappa = kappa, nsim = nsim,
                messages.screen = FALSE, lambda = obj.variog$lambda)
  if(messages.screen)
    cat("variog.env: adding the mean or trend\n")
  x.mat <- trend.spatial(trend=obj.variog$trend, coords=coords)
  simula$data <- as.vector(x.mat %*% beta) + simula$data
  ##
  ## computing empirical variograms for the simulations
  ##
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if(is.R()){
    bin.f <- function(sim){
      cbin <- vbin <- sdbin <- rep(0, nbins)  
      temp <- .C("binit",
                 as.integer(obj.variog$n.data),
                 as.double(as.vector(coords[,1])),
                 as.double(as.vector(coords[,2])),
                 as.double(as.vector(sim)),
                 as.integer(nbins),
                 as.double(as.vector(obj.variog$bins.lim)),
                 as.integer(estimator.type == "robust"),
                 as.double(max(obj.variog$u)),
                 as.double(cbin),
                 vbin = as.double(vbin),
                 as.integer(FALSE),
                 as.double(sdbin)
                 )$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
    simula.bins <- simula.bins[obj.variog$ind.bin,]
  }
  else{
    bin.f <- function(sim, nbins, n.data, coords, bins.lim, estimator.type, max.u){
      cbin <- vbin <- sdbin <- rep(0, nbins)  
      temp <- .C("binit",
                 as.integer(n.data),
                 as.double(as.vector(coords[,1])),
                 as.double(as.vector(coords[,2])),
                 as.double(as.vector(sim)),
                 as.integer(nbins),
                 as.double(as.vector(bins.lim)),
                 as.integer(estimator.type == "robust"),
                 as.double(max.u),
                 as.double(cbin),
                 vbin = as.double(vbin),
                 as.integer(FALSE),
                 as.double(sdbin)
                 )$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f, nbins=nbins, n.data=obj.variog$n.data, coords=coords, bins.lim=obj.variog$bins.lim, estimator.type=estimator.type, max.u=max(obj.variog$u))
    simula.bins <- simula.bins[obj.variog$ind.bin,]
  }
  if(save.sim == FALSE) simula$data <- NULL
  ##
  ## computing envelops
  ##
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, range)
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ],
                  v.upper = limits[2,])
  if(save.sim) res.env$simulated <- simula$data
  res.env$call <- call.fc
  class(res.env) <- "variogram.envelope"
  return(res.env)
}

"variog.mc.env" <-
  function (geodata, coords = geodata$coords, data = geodata$data,
            obj.variog, nsim = 99, save.sim = FALSE,
            messages.screen = TRUE) 
{
  call.fc <- match.call()
  ##
  ## Checking input
  ##
  obj.variog$v <- NULL
  if((is.matrix(data) | is.data.frame(data)))
    if(ncol(data) > 1)
      stop("envelops can be computed for only one data set at once")
  if(!is.null(obj.variog$estimator.type))
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  ##
  ## generating several "data-sets" by permutating data values
  ##
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim,
              "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data){return(data[sample(1:n.data)])}
  simula$data <- apply(as.matrix(1:nsim),1, perm.f, data=data, n.data=n.data) 
  ##
  ## computing empirical variograms for the simulations
  ##
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if(is.R()){
    bin.f <- function(sim){
      cbin <- vbin <- sdbin <- rep(0, nbins)  
      temp <- .C("binit",
                 as.integer(obj.variog$n.data),
                 as.double(as.vector(coords[,1])),
                 as.double(as.vector(coords[,2])),
                 as.double(as.vector(sim)),
                 as.integer(nbins),
                 as.double(as.vector(obj.variog$bins.lim)),
                 as.integer(estimator.type == "robust"),
                 as.double(max(obj.variog$u)),
                 as.double(cbin),
                 vbin = as.double(vbin),
                 as.integer(FALSE),
                 as.double(sdbin)
                 )$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
    simula.bins <- simula.bins[obj.variog$ind.bin,]
  }
  else{
    bin.f <- function(sim, nbins, n.data, coords, bins.lim, estimator.type, max.u){
      cbin <- vbin <- sdbin <- rep(0, nbins)  
      temp <- .C("binit",
                 as.integer(n.data),
                 as.double(as.vector(coords[,1])),
                 as.double(as.vector(coords[,2])),
                 as.double(as.vector(sim)),
                 as.integer(nbins),
                 as.double(as.vector(bins.lim)),
                 as.integer(estimator.type == "robust"),
                 as.double(max.u),
                 as.double(cbin),
                 vbin = as.double(vbin),
                 as.integer(FALSE),
                 as.double(sdbin)
                 )$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f, nbins=nbins, n.data=obj.variog$n.data, coords=coords, bins.lim=obj.variog$bins.lim, estimator.type=estimator.type, max.u=max(obj.variog$u))
    simula.bins <- simula.bins[obj.variog$ind.bin,]
  }
  if(save.sim == FALSE) simula$data <- NULL
  ##
  ## computing envelops
  ##
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, range)
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ],
                  v.upper = limits[2,])
  if(save.sim) res.env$simulations <- simula$data
  res.env$call <- call.fc
  class(res.env) <- "variogram.envelope"
  return(res.env)
}

"lines.variogram.envelope" <-
  function(obj, lty=3, ...)
{
  lines(obj$u, obj$v.lower, ...)
  lines(obj$u, obj$v.upper, ...)
  return(invisible())
}

"varcov.spatial" <-
  function (coords = NULL, dists.lowertri = NULL, cov.model = "matern",
            kappa = 0.5, nugget = 0, cov.pars = stop("no cov.pars argument"), 
            inv = FALSE, det = FALSE,
            func.inv = c("cholesky", "eigen", "svd", "solve"),
            scaled = FALSE, only.decomposition = FALSE, 
            sqrt.inv = FALSE, try.another.decomposition = TRUE,
            only.inv.lower.diag = FALSE) 
{
  if (is.R()) 
    require(mva)
  func.inv <- match.arg(func.inv)
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if (only.inv.lower.diag) 
    inv <- TRUE
  if (is.null(coords) & is.null(dists.lowertri)) 
    stop("one of the arguments, coords or dists.lowertri must be provided")
  if (!is.null(coords) & !is.null(dists.lowertri)) 
    stop("only ONE argument, either coords or dists.lowertri must be provided")
  if (!is.null(coords)) 
    n <- nrow(coords)
  if (!is.null(dists.lowertri)) {
    n <- as.integer(round(0.5 * (1 + sqrt(1 + 8 * length(dists.lowertri)))))
  }
  tausq <- nugget
  if (is.vector(cov.pars)) {
    sigmasq <- cov.pars[1]
    phi <- cov.pars[2]
  }
  else {
    sigmasq <- cov.pars[, 1]
    phi <- cov.pars[, 2]
  }
##  print(c(tausq=tausq, sigmasq=sigmasq, phi=phi, kappa=kappa))
  if (!is.null(coords)) {
    dists.lowertri <- as.vector(dist(coords))
  }
  if (round(1e+12 * min(dists.lowertri)) == 0) 
    warning("Two or more pairs of data at coincident (or very close) locations. \nThis can cause matrices operations to crash!\n")
  varcov <- matrix(0, n, n)
  if (scaled) {
    if (all(phi < 1e-12)) 
      varcov <- diag(x = (1 + (tausq/sum(sigmasq))), n)
    else {
      if (is.vector(cov.pars)) 
        cov.pars.sc <- c(1, phi)
      else cov.pars.sc <- cbind(1, phi)
      covvec <- cov.spatial(obj = dists.lowertri, cov.model = cov.model, 
                            kappa = kappa, cov.pars = cov.pars.sc)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      if (is.R()) remove("covvec")
      else remove("covvec", frame = sys.nframe())
      diag(varcov) <- 1 + (tausq/sum(sigmasq))
    }
  }
  else {
    if (all(sigmasq < 1e-10) | all(phi < 1e-10)) {
      varcov <- diag(x = (tausq + sum(sigmasq)), n)
    }
    else {
      covvec <- cov.spatial(obj = dists.lowertri, cov.model = cov.model, 
                            kappa = kappa, cov.pars = cov.pars)
      varcov[lower.tri(varcov)] <- covvec
      varcov <- t(varcov)
      varcov[lower.tri(varcov)] <- covvec
      if (is.R()) remove("covvec")
      else remove("covvec", frame = sys.nframe())
      diag(varcov) <- tausq + sum(sigmasq)
    }
  }
  if (inv | det | only.decomposition | sqrt.inv | only.inv.lower.diag) {
    if (func.inv == "cholesky") {
      options(show.error.messages = FALSE)
      varcov.sqrt <- try(chol(varcov))
      options(show.error.messages = TRUE)
      if (!is.numeric(varcov.sqrt)) {
        if (try.another.decomposition) 
          func.inv <- "eigen"
        else {
          print(varcov.sqrt[1])
          stop()
        }
      }
      else {
        if (only.decomposition | inv) 
          if (is.R()) remove("varcov")
          else remove("varcov", frame = sys.nframe())
        if (only.decomposition == FALSE) {
          if (det) 
            cov.logdeth <- sum(log(diag(varcov.sqrt)))
          if (sqrt.inv) 
            inverse.sqrt <- solve(varcov.sqrt)
          if (inv) {
            if (is.R()) {
              invcov <- chol2inv(varcov.sqrt)
              if (!sqrt.inv)
                remove("varcov.sqrt")
            }
            else {
              invcov.sqrt <- solve.upper(varcov.sqrt)
              invcov <- invcov.sqrt %*% t(invcov.sqrt)
              if (!sqrt.inv) 
                remove("varcov.sqrt", frame = sys.nframe())
            }
          }
        }
      }
    }
    if (func.inv == "svd") {
      varcov.svd <- svd(varcov, nv = 0)
      options(show.error.messages = FALSE)
      cov.logdeth <- try(sum(log(sqrt(varcov.svd$d))))
      options(show.error.messages = TRUE)
      if (is.numeric(cov.logdeth)) {
        if (try.another.decomposition) 
          func.inv <- "eigen"
        else {
          print(cov.logdeth[1])
          stop()
        }
      }
      else {
        if (only.decomposition | inv) 
          if (is.R())  remove("varcov")
          else remove("varcov", frame = sys.nframe())
        if (only.decomposition) 
          varcov.sqrt <- t(varcov.svd$u %*% (t(varcov.svd$u) * 
                                             sqrt(varcov.svd$d)))
        if (inv) {
          invcov <- t(varcov.svd$u %*% (t(varcov.svd$u) * 
                                        (1/varcov.svd$d)))
        }
        if (sqrt.inv) 
          inverse.sqrt <- t(varcov.svd$u %*% (t(varcov.svd$u) * 
                                              (1/sqrt(varcov.svd$d))))
      }
    }
    if (func.inv == "solve") {
      if (det) 
        stop("the option func.inv == \"solve\" does not allow computation of determinants. \nUse func.inv = \"chol\",\"svd\" or \"eigen\"\n")
      options(show.error.messages = FALSE)
      invcov <- try(solve(varcov))
      options(show.error.messages = TRUE)
      if (is.numeric(cov.logdeth)) {
        if (try.another.decomposition) 
          func.inv <- "eigen"
        else {
          print(invcov[1])
          stop()
        }
      }
      if (is.R()) remove("varcov")
      else remove("varcov", frame = sys.nframe())
    }
    if (func.inv == "eigen") {
      options(show.error.messages = FALSE)
      varcov.eig <- try(eigen(varcov, symmetric = TRUE))
      cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))))
      options(show.error.messages = TRUE)
      if (!is.numeric(cov.logdeth) | is.numeric(varcov.eig)) {
        diag(varcov) <- 1.0001 * diag(varcov)
        options(show.error.messages = FALSE)
        varcov.eig <- try(eigen(varcov, symmetric = TRUE))
        cov.logdeth <- try(sum(log(sqrt(varcov.eig$val))))
        options(show.error.messages = TRUE)
        if (!is.numeric(cov.logdeth) | !is.numeric(varcov.eig)) {
          return(list(crash.parms = c(tausq=tausq, sigmasq=sigmasq, phi=phi, kappa=kappa)))
        }
      }
      else {
        if (only.decomposition | inv) 
          if (is.R()) remove("varcov")
          else remove("varcov", frame = sys.nframe())
        if (only.decomposition) 
          varcov.sqrt <- (varcov.eig$vec %*% diag(sqrt(varcov.eig$val)) %*% 
                          t(varcov.eig$vec))
        if (inv) 
          invcov <- (varcov.eig$vec %*% diag(1/varcov.eig$val) %*% 
                     t(varcov.eig$vec))
        if (sqrt.inv) 
          inverse.sqrt <- (varcov.eig$vec %*% diag(1/sqrt(varcov.eig$val)) %*% 
                           t(varcov.eig$vec))
      }
    }
  }
  if (only.decomposition == FALSE) {
    if (det) {
      if (inv) {
        if (only.inv.lower.diag) 
          result <- list(lower.inverse = invcov[lower.tri(invcov)], 
                         diag.inverse = diag(invcov), log.det.to.half = cov.logdeth)
        else result <- list(inverse = invcov, log.det.to.half = cov.logdeth)
      }
      else {
        result <- list(varcov = varcov, log.det.to.half = cov.logdeth)
      }
      if (sqrt.inv) 
        result$sqrt.inverse <- inverse.sqrt
    }
    else {
      if (inv) {
        if (only.inv.lower.diag) 
          result <- list(lower.inverse = invcov[lower.tri(invcov)], 
                         diag.inverse = diag(invcov))
        else {
          if (sqrt.inv) 
            result <- list(inverse = invcov, sqrt.inverse = inverse.sqrt)
          else result <- list(inverse = invcov)
        }
      }
      else result <- list(varcov = varcov)
    }
  }
  else result <- list(sqrt.varcov = varcov.sqrt)
  result$crash.parms <- NULL
  return(result)
}

"plot.geodata" <-
  function (geodata, coords = geodata$coords, data = geodata$data, 
            trend = "cte", lambda = 1, col.data = 1,
            weights.divide = NULL, window.new = FALSE, ...) 
{
  if (is.R()) 
    par.ori <- par(no.readonly = TRUE)
  else par.ori <- par()
  on.exit(par(par.ori))
  coords <- as.matrix(geodata$coords)
  data <- as.matrix(data)
  data <- data[, col.data]
  if (window.new) {
    if (is.R()) 
      X11()
    else trellis.device()
  }
  if (lambda != 1) {
    if (lambda == 0) 
      data <- log(data)
    else data <- ((data^lambda) - 1)/lambda
  }
  if (!is.null(weights.divide)) {
    if (length(weights.divide) != length(data)) 
      stop("length of weights.divide must be equals to the length of data")
    data <- data/weights.divide
  }
  par(mfrow = c(2, 2))
  par(mar = c(4, 4, 0, 0.5))
  data.quantile <- quantile(data)
  coords.lims <- apply(coords, 2, range)
  coords.diff <- diff(coords.lims)
  if (coords.diff[1] != coords.diff[2]) {
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[, ind.min] <- coords.lims[, ind.min] + c(-coords.diff.diff, 
                                                         coords.diff.diff)/2
  }
  par(pty = "s")
  plot(coords, xlab = "Coord X", ylab = "Coord Y", type = "n", 
       xlim = coords.lims[, 1], ylim = coords.lims[, 2])
  if (is.R()) {
    data.cut <- cut(data, breaks=quantile(data), include.l=TRUE, labels=FALSE)
#    points(coords, cex = c(0.4,0.6, 0.9, 1.2)[data.cut], pch=21, bg=c("blue", "green", "yellow2", "red")[data.cut])
    points(coords, pch = (1:4)[data.cut], col=c("blue", "green", "yellow2", "red")[data.cut])
  }
  else {
    points(coords[(data <= data.quantile[2]), ], pch = 1, 
           cex = 0.6, col = 2)
    points(coords[((data > data.quantile[2]) & (data <= 
                                                data.quantile[3])), ], pch = 2, cex = 1.4, col = 4)
    points(coords[((data > data.quantile[3]) & (data <= 
                                                data.quantile[4])), ], pch = 3, cex = 1.7, col = 7)
    points(coords[(data > data.quantile[4]), ], pch = 4, 
           cex = 2, col = 8)
  }
  par(pty = "m")
  if (is.R()) 
    par(mar = c(4, 4, 1, 1))
  else par(mar = c(0, 1, 0, 0.5))
  if (is.R()) {
    if (require(scatterplot3d) == FALSE) {
      hist(data)
      warning("plot.geodata: a 3d plot will be draw instead of the histogram if the package \"scatterplot3d\" is available")
    }
    else scatterplot3d(x = coords[, 1], y = coords[, 2], 
                       z = data, box = F, type = "h", pch = 16, xlab = "Coord X", 
                       ylab = "Coord Y", ...)
  }
  else xyzplot(coords = coords, data = data, ...)
  if (!is.R()) 
    par(mar = c(5, 5, 1, 0.5))
  plot(coords[, 1], data, xlab = "Coord X", cex = 1)
  plot(coords[, 2], data, xlab = "Coord Y", cex = 1)
  return(invisible())
}


