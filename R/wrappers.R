##
## "wrappers" for pieces of C code used in geoR/geoS
##
bilinearformXAY <-
  function(X, lowerA, diagA, Y)
  {
    nA <- length(diagA)
    nX <- length(X)/nA
    nY <- length(Y)/nA
    if(length(lowerA) != (nA * (nA -1)/2))
      stop("lowerA and diagA have incompatible dimentions")
    out <- .C("bilinearform_XAY",
              as.double(as.vector(lowerA)),
              as.double(as.vector(diagA)),
              as.double(as.vector(X)),
              as.double(as.vector(Y)),
              as.integer(nX),
              as.integer(nY),
              as.integer(nA),
              res=as.double(rep(0,(nX*nY))))$res
    attr(out, "dim") <- c(nX, nY)
    return(out)
  }

diagquadraticformXAX <-
  function(X, lowerA, diagA)
  {
    nA <- length(diagA)
    nX <- length(X)/nA
    if(length(lowerA) != (nA * (nA -1)/2))
      stop("lowerA and diagA have incompatible dimentions")
    out <- .C("diag_quadraticform_XAX",
              as.double(as.vector(lowerA)),
              as.double(as.vector(diagA)),
              as.double(as.vector(X)),
              as.integer(nX),
              as.integer(nA),
              res = as.double(rep(0,nX)))$res
    return(out)
  }

loccoords <-
  function(coords, locations)
  {
    ## Computes a matrix for which each row has the distances between
    ## each point in 'locations' to all the points in 'coords'
    coords <- as.matrix(coords)
    locations <- as.matrix(locations)
    dimc <- dim(coords)
    diml <- dim(locations)
    if((dimc[2] != 2) | (diml[2] != 2))
      stop("coords and locations must have two columns")
    nc <- dimc[1]
    nl <- diml[1]
    out <- as.double(rep(0, nc*nl))
    .C("loccoords",
       as.double(as.vector(locations[,1])),
       as.double(as.vector(locations[,2])),
       as.double(as.vector(coords[,1])),
       as.double(as.vector(coords[,2])),
       as.integer(nl),
       as.integer(nc),
       out, DUP=FALSE)
    attr(out, "dim") <- c(nc, nl)
    return(out)
  }

distdiag <-
  function(coords)
  {
    ## returns the lower triangle of the matrix with euclidean distances
    ## between pairs of points, including the diagonal. 
    ##
    coords <- as.matrix(coords)
    dimc <- dim(coords)
    if(dimc[2] == 1 & dimc[1] == 2)
      return(0)
    else{
      if(dimc[2] != 2)
        stop("coords must have two columns")
      nc <- dimc[1]
      out <- as.double(rep(0, (nc * (nc+1)/2)))
      .C("distdiag",
         as.double(coords[,1]),
         as.double(coords[,2]),
         as.integer(nc),
         out, DUP = FALSE)
      return(out)
    }
  }

"corr.diaglowertri" <-
  function(coords, cov.model, phi, kappa)
{
  cov.model <- match.arg(cov.model,
                         choices = c("matern", "exponential", "gaussian",
                           "spherical", "circular", "cubic", "wave", "power",
                           "powered.exponential", "cauchy", "gneiting",
                           "gneiting.matern", "pure.nugget"))
  if(any(cov.model == c("cauchy", "matern", "powered.exponential", "power", "gneiting.matern"))){
    if(missing(kappa))
      stop("argument kappa is needed for this choice of correlation function")
  }
  else kappa <- 1
  coords <- as.matrix(coords)
  dimc <- dim(coords)
  if(dimc[2] == 1 & dimc[1] == 2)
    return(0)
  else{
    if(dimc[2] != 2)
      stop("coords must have two columns")
    nc <- dimc[1]
    out <- as.double(rep(0, (nc * (nc+1)/2)))
    .C("cor_diag",
       as.double(coords[,1]),
       as.double(coords[,2]),
       as.integer(nc),
       as.integer(cor.number(cov.model)),
       as.double(phi),
       as.double(kappa),
       out, DUP = FALSE)
    return(out)
  }
}

"cond.sim" <-
  function(env.loc, env.iter, loc.coincide, tmean, Rinv, mod, vbetai,
           fixed.sigmasq)
  {
    NTOT <- mod$nloc * mod$Nsims
    if(fixed.sigmasq)
      invchisc <- rep(1, NTOT)
    else
      invchisc <- sqrt(mod$df.model/rchisq(mod$Nsims, df=mod$df.model))
    ##
    if(mod$beta.size == 1){
      Blower <- 0
      Bdiag <- vbetai
    }
    else{
      Blower <- vbetai[lower.tri(vbetai)]
      Bdiag <- diag(vbetai)
    }
    ##
    if((length(tmean) %% mod$nloc) > 0)
      stop("cond.sim: wrong size of tmean")
    tmean <- matrix(tmean, nrow = mod$nloc)
    ncol.tmean <- ncol(tmean)
    if(ncol(tmean) > 1){
      if(ncol.tmean != mod$Nsims)
        stop("cond.sim: size of tmean does not matches with Nsims")
      require(geoRglmm)
      diff.mean <- as.integer(1)
    }      
    else
      diff.mean <- as.integer(0)
    ##
    normalsc <- rnorm(NTOT)
    if(is.null(loc.coincide))
      locations <- get("locations", envir=env.loc)
    else
      locations <- get("locations", envir=env.loc)[-loc.coincide,,drop=FALSE]
    ##
    ##
##    R0 <- varcov.spatial(coords = locations,
##                         cov.pars=c(mod$Dval, mod$phi))[[1]]
##    iR <- matrix(0, mod$n,mod$n)
##    iR[lower.tri(iR)] <- Rinv$lower
##    iR <- iR  + t(iR)
##    diag(iR) <- Rinv$diag
##    v0iRv0 <- crossprod(get("v0", envir=env.iter), iR) %*%
##      get("v0", envir=env.iter)
##    V <- vbetai
##    bVb <- t(get("b", envir=env.iter)) %*% V %*% get("b", envir=env.iter)
##    Vmat <- R0 - v0iRv0  + bVb
##    Vchol <- try(chol(Vmat))
##    if(!is.numeric(Vchol)) print(try(chol(Vmat)))
    ##
    ##
    if(is.null(loc.coincide))
      loccoin <- TRUE
    else loccoin <- -loc.coincide
    normalsc <- .C("kb_sim_new",
                   as.double(as.vector(tmean)),
                   out = as.double(normalsc),
                   as.double(as.vector(Rinv$lower)),
                   as.double(as.vector(Rinv$diag)),
                   as.double(as.vector(get("v0", envir=env.iter))),
                   as.integer(mod$nloc),
                   as.integer(mod$n),
                   as.double(mod$Dval),
                   as.integer(mod$Nsims),
                   as.double(invchisc),
                   as.double(mod$s2),
                   as.double(Blower),
                   as.double(Bdiag),
                   as.double(as.vector(get("b", envir=env.iter))),
                   as.integer(mod$beta.size),
                   as.double(get("locations", envir=env.loc)[loccoin,1]),
                   as.double(get("locations", envir=env.loc)[loccoin,2]),
                   as.integer(mod$cov.model.number),
                   as.double(mod$phi),
                   as.double(mod$kappa),
                   as.integer(diff.mean))$out      
    attr(normalsc, "dim") <- c(mod$nloc, mod$Nsims)
    return(normalsc)
  }

