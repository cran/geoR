".First.lib" <-
  function(lib, pkg)
{
  library.dynam("geoR", package = pkg, lib.loc = lib)  
  cat("\n")
  cat("-------------------------------------------------\n")
  if(is.R()){
    cat(package.description("geoR", lib = lib, field="Title"))
    cat("\n")
    ver <- package.description("geoR", lib = lib, field="Version")
    dat <- package.description("geoR", lib = lib, field="Date")
    cat(paste("geoR version ", ver," (", dat, ") is now loaded\n", sep=""))
  }
  else{
    cat("geoS: a package for geostatistical analysis in S-PLUS\n")
    cat("geoS is now loaded\n")
  }
  cat("-------------------------------------------------\n")
  cat("\n")
  return(invisible(0))
}

