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













