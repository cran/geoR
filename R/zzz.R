".First.lib" <-
  function(lib, pkg)
{
  library.dynam("geoR", package = pkg, lib.loc = lib)  
  cat("\n")
  cat("------------------------------------------------\n")
  cat("geoR: a package for geostatistical analysis in R\n")
  cat("geoR is now loaded\n")
  cat("------------------------------------------------\n")
  cat("\n")
  cat("\n")
}

