##
## Auxiliary functions for the geoR library
## ----------------------------------------
##
## These functions are typically called by other main functions
## to perform internal calculations
##

"solve.geoR" <-
  function (a, b = NULL, ...) 
{
  require(methods)
  a <- eval(a)
  b <- eval(b)
  if(exists("trySilent")){
    if (is.null(b)) res <- trySilent(solve(a, ...))
    else res <- trySilent(solve(a, b, ...))
  }
  else{
    error.now <- options()$show.error.messages
    if (is.null(error.now) | error.now) 
      on.exit(options(show.error.messages = TRUE))
    options(show.error.messages = FALSE)
    if (is.null(b)) res <- try(solve(a, ...))
    else res <- try(solve(a, b, ...))
  }
  if (inherits(res, "try-error")) {
    test <- all.equal.numeric(a, t(a), 100 * .Machine$double.eps)
    if(!(is.logical(test) && test)){
      ##      options(show.error.messages = TRUE)
      stop("matrix `a' is not symmetric")
    }
    t.ei <- eigen(a, symmetric = TRUE)
    if(exists("trySilent")){
      if (is.null(b)) res <- trySilent(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec))
      else res <- trySilent(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec) %*% b)
    }
    else{
      if (is.null(b)) res <- try(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec))
      else res <- try(t.ei$vec %*% diag(t.ei$val^(-1)) %*% t(t.ei$vec) %*% b)
    }
    if (any(is.na(res)) | any(is.nan(res)) | any(is.infinite(res))) 
      class(res) <- "try-error"
  }
  if (inherits(res, "try-error")) 
    stop("Singular matrix. Covariates may have different orders of magnitude.")
  return(res)
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

#"dist0.krige" <-
#function (x0, coords) 
#{
#  if (length(x0) != 2) 
#    stop(paste("length of x0 is", length(x0), "(it must be 2)"))
#  coords[, 1] <- coords[, 1] - x0[1]
#  coords[, 2] <- coords[, 2] - x0[2]
#  return(sqrt(coords[, 1]^2 + coords[, 2]^2))
#}


"set.coords.lims" <-
  function(coords, xlim, ylim)
{
  coords.lims <- apply(coords, 2, range)
  if(!missing(xlim) && is.numeric(xlim)) coords.lims[,1] <- xlim[order(xlim)]
  if(!missing(ylim) && is.numeric(ylim)) coords.lims[,2] <- ylim[order(ylim)]
  coords.diff <- diff(coords.lims)
  if (coords.diff[1] != coords.diff[2]) {
    coords.diff.diff <- abs(diff(as.vector(coords.diff)))
    ind.min <- which(coords.diff == min(coords.diff))
    coords.lims[, ind.min] <-
      coords.lims[, ind.min] +
        c(-coords.diff.diff, coords.diff.diff)/2
  }
  return(coords.lims)
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


"polygrid" <- 
  function(xgrid, ygrid, borders, vec.inout = FALSE)
{
  ## checking for splancs
  if(is.R())
    if(!require(splancs))
      cat("ERROR: cannot run the function\npackage \"splancs\" should be installed/loaded")
    else library(splancs)
  ## checking input
  if(!is.list(xgrid) && is.vector(drop(xgrid))){
    if(missing(ygrid)) stop("xgrid must have x and y coordinates or a vector must be provided for ygrid")
    if(!is.vector(ygrid)) stop("ygrid must be a vector")
    xygrid <- expand.grid(x = xgrid, y = ygrid)
  }
  if(is.matrix(xgrid) || is.data.frame(xgrid)){
    if(ncol(xgrid) != 2) stop("xgrid must be a vector or a 2 column matrix or data-frame")
    xygrid <- xgrid
    if(!missing(xgrid)) stop("xgrid has 2 column, ygrid was ignored")
  }
  else
    if(is.list(xgrid)){
      if(length(xgrid) != 2) stop("if xgrid is a list it must have 2 elements")
      xygrid <- expand.grid(x = xgrid[[1]], y = xgrid[[2]])
      if(!missing(xgrid)) stop("xgrid is a list, ygrid was ignored")
    }
  if(nrow(borders) < 3) stop("borders must have at least 3 points")
  if(exists("inout")){
    ind <- as.vector(inout(pts=xygrid, poly=borders))
    xypoly <- xygrid[ind == TRUE,  ]
    if(vec.inout == FALSE)
      return(xypoly)
    else return(list(xypoly = xypoly, vec.inout = ind))
  }
  else{
    cat("ERROR: cannot run the function\n")
    cat("package \"splancs\" should be installed/loaded")
    return(invisible())
  }
}

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
    require(methods)
    if(exists("trySilent")){
      trend.mat <- trySilent(model.matrix(trend))
    }
    else{
      error.now <- options()$show.error.messages
      if (is.null(error.now) | error.now) 
        on.exit(options(show.error.messages = TRUE))
      options(show.error.messages = FALSE)
      trend.mat <- try(model.matrix(trend))
    }    
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

"nlmP" <- function(objfunc, params, lower = rep( -Inf, length(params)),
                   upper = rep(+Inf, length(params)), ... )
{
  ## minimizer, using nlm with transformation of variables
  ## to allow for limits for the parameters   
  ##
  ## objfunc is a function to be optimised
  ## params is a starting value for the parameters
  ##
  ## NOTE: this function was used before optim() becomes available for R
  ##       It has limited usage now.
  ##
  ## Adapted from a function from Patrick E. Brown, Lancaster University 
  ##
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



