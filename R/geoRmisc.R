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
    function(coords, xzoom, yzoom=xzoom, xlim.ori, ylim.ori, xoff=0, yoff=0)
{
  if(missing(ylim.ori)) xlim.ori <- range(coords[,1], na.rm=TRUE)
  if(missing(ylim.ori)) ylim.ori <- range(coords[,2], na.rm=TRUE)
  xlim <- xlim.ori + c(-1,1) * (diff(xlim.ori)/2) * (xzoom - 1)
  ylim <- ylim.ori + c(-1,1) * (diff(ylim.ori)/2) * (yzoom - 1)
  res <- coords2coords(coords, xlim=xlim, ylim=ylim, xlim.ori = xlim.ori, ylim.ori=ylim.ori)
  res[,1] <- res[,1] + xoff
  res[,2] <- res[,2] + yoff
  return(res)
}

"rect.coords" <-
  function(coords, xzoom = 1, yzoom=xzoom, add.to.plot=TRUE, ...)
{
  rx <- range(coords[,1], na.rm=TRUE)
  ry <- range(coords[,2], na.rm=TRUE)
  res <- cbind(c(rx,rev(rx)), rep(ry,c(2,2)))
  res <- zoom.coords(res, xzoom=xzoom, yzoom=yzoom)
  if(add.to.plot) rect(res[1,1], res[1,2], res[3,1], res[3,2], ...)
  return(res)  
}
