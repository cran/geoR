##
## Miscelaneous geoR functions
##

cite.geoR <- function()
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
