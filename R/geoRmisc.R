##
## Miscelaneous geoR functions
##

"cite.geoR" <- function()
{
    cat("\n")
    cat("To cite geoR in publications, use\n\n")
    msg <- "RIBEIRO Jr., P.J. & DIGGLE, P.J. (2001) geoR: A package for geostatistical analysis. R-NEWS, Vol 1, No 2, 15-18. ISSN 1609-3631."
    writeLines(strwrap(msg, prefix = "  "))
    cat("\n")
    msg <- paste("Please cite geoR when using it for data analysis!")
    writeLines(strwrap(msg))
    cat("\nA BibTeX entry for LaTeX users is\n\n")
    cat("  @Article{,\n")
    cat("     title	   = {{geoR}: a package for geostatistical analysis},\n")
    cat("     author        = {Ribeiro Jr., P.J. and Diggle, P.J.},\n")
    cat("     journal       = {R-NEWS},\n")
    cat("     year	   = {2001},\n")
    cat("     volume	   = {1},\n")
    cat("     number	   = {2},\n")
    cat("     pages	   = {15--18},\n")
    cat("     issn          = {1609-3631},\n")
    cat("     url           = {http://cran.R-project.org/doc/Rnews}\n")
    cat("   }\n\n")
}

#geoR.options <- function(messages = TRUE, ...)
#{
#  res <- list(...)
#  res$messages <- messages
#  .geoR.options <<- res
#  return(invisible())
#}

"coords2coords" <-
  function(coords, xlim, ylim, xlim.ori, ylim.ori)
{
  if(missing(ylim.ori)) xlim.ori <- range(coords[,1], na.rm=TRUE)
  if(missing(ylim.ori)) ylim.ori <- range(coords[,2], na.rm=TRUE)
  coords[,1] <- xlim[1] + (coords[,1] - xlim.ori[1]) * diff(xlim)/diff(xlim.ori)
  coords[,2] <- ylim[1] + (coords[,2] - ylim.ori[1]) * diff(ylim)/diff(ylim.ori)
  return(coords)
}

"zoom.coords" <-
  function(x, ...)
{
  UseMethod("zoom.coords")
}

"zoom.coords.default" <-
    function(x, xzoom, yzoom=xzoom, xlim.ori, ylim.ori, xoff=0, yoff=0, ...)
{
  if(missing(ylim.ori)) xlim.ori <- range(x[,1], na.rm=TRUE)
  if(missing(ylim.ori)) ylim.ori <- range(x[,2], na.rm=TRUE)
  xlim <- xlim.ori + c(-1,1) * (diff(xlim.ori)/2) * (xzoom - 1)
  ylim <- ylim.ori + c(-1,1) * (diff(ylim.ori)/2) * (yzoom - 1)
  res <- coords2coords(x, xlim=xlim, ylim=ylim, xlim.ori = xlim.ori, ylim.ori=ylim.ori)
  res[,1] <- res[,1] + xoff
  res[,2] <- res[,2] + yoff
  return(res)
}

"zoom.coords.geodata" <-
  function(x, ...)
{
  x$coords <- zoom.coords.default(x$coords, ...)
  if(!is.null(x$borders)) x$borders <- zoom.coords.default(x$borders, ...)
  if(!is.null(x$subarea.lims)) x$subarea.lims <- zoom.coords.default(x$subarea.lims, ...)
  return(x)
}

"rect.coords" <-
  function(coords, xzoom = 1, yzoom=xzoom, add.to.plot=TRUE, quiet=FALSE, ...)
{
  rx <- range(coords[,1], na.rm=TRUE)
  ry <- range(coords[,2], na.rm=TRUE)
  res <- cbind(c(rx,rev(rx)), rep(ry,c(2,2)))
  res <- zoom.coords(res, xzoom=xzoom, yzoom=yzoom)
  if(add.to.plot) rect(res[1,1], res[1,2], res[3,1], res[3,2], ...)
  if(quiet)
    return(invisible())
  else
    return(res)  
}

# make it generic with a method for geoR
"dup.coords" <-
  function(x)
{
  UseMethod("dup.coords")
}

"dup.coords.default" <-
  function(x)
{
  ap1 <- unclass(factor(paste("x",x[,1],"y",x[,2], sep="")))
  ap2 <- table(ap1)
  ap2 <- ap2[ap2 > 1]
  takecoords <- function(n){
    if(!is.null(rownames(x))) rownames(x[ap1 == n,])
    else (1:nrow(x))[ap1 == n]
    }
  res <- sapply(as.numeric(names(ap2)), takecoords)
  if(length(res) == 0) res <- NULL
  return(res)
}

"dup.coords.geodata" <-
  function(x)
{
  return(dup.coords.default(x$coords))
}

"subarea" <-
  function(geodata, xlim, ylim, ...)
{
  if(class(geodata) != "geodata")
    stop("an object of the class geodata must be provided")
  if(missing(xlim) & !missing(ylim)) xlim <- c(-Inf, +Inf)
  if(!missing(xlim) & missing(ylim)) ylim <- c(-Inf, +Inf)
  if(missing(xlim) & missing(ylim)){
    cat("Enter 2 points defining the corners of the subarea\n")
    pt <- locator(2)
    xlim <- sort(pt[[1]])
    ylim <- sort(pt[[2]])
  }
  if(!is.vector(xlim) || length(xlim) != 2)
    stop("xlim must be a vector with 2 elements")
  if(!is.vector(ylim) || length(ylim) != 2)
    stop("ylim must be a vector with 2 elements")
  geo <- geodata
  ind <- (geodata$coords[,1] > xlim[1] & geodata$coords[,1] < xlim[2] &  
          geodata$coords[,2] > ylim[1] & geodata$coords[,2] < ylim[2])  
  geo$coords <- geodata$coords[ind,]
  xlim.all <- c(xlim, range(geo$coords[,1]))
  ylim.all <- c(ylim, range(geo$coords[,2]))
  if(is.vector(geodata$data)) geo$data <- geodata$data[ind]
  else geo$data <- geodata$data[ind,]
  geo$units.m <- geodata$units.m[ind]
  if(!is.null(geodata$covariate)){
    if(is.vector(geodata$covariate))
      geo$covariate <- geodata$covariate[ind]
    if(is.matrix(geodata$covariate) |  is.data.frame(geodata$covariate))   
      geo$covariate <- geodata$covariate[ind,]
  }
  if(!is.null(geodata$borders)){
    geo$borders <- geodata$borders[(geodata$borders[,1]>xlim[1] & geodata$borders[,1]<xlim[2] &  
                                    geodata$borders[,2]>ylim[1] & geodata$borders[,2] < ylim[2]),]
    xlim.all <- c(xlim.all, range(geo$borders[,1]))
    ylim.all <- c(ylim.all, range(geo$borders[,2]))
  }
  geo$subarea.lims <- cbind(range(xlim.all[is.finite(xlim.all)]),
                            range(ylim.all[is.finite(ylim.all)]))
  return(geo)
}

