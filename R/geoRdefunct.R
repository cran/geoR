"geoRdefunct" <-
  function()
  {
    cat("\n")
    cat("The following functions are no longer used in geoR:")
    cat("---------------------------------------------------")
    cat("\nolsfit: use variofit() instead")
    cat("\nwlsfit: use variofit() instead")
    cat("\nlikfit.old: use likfit() instead")
    cat("\n")
  }
    
"olsfit" <- function(...)
  stop("this function is now obsolete.\nuse variofit() instead.")

"wlsfit" <- function(...)
  stop("this function is now obsolete.\nuse variofit() instead.")


"distdiag" <-
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
         out, DUP = FALSE,
         PACKAGE = "geoR")
      return(out)
    }
  }

