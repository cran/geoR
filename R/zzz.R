".First.lib" <-
  function(lib, pkg)
{
  library.dynam("geoR", package = pkg, lib.loc = lib)  
  cat("\n")
  cat("-------------------------------------------------\n")
  cat(package.description("geoR", lib = lib, field="Title"))
  cat("\n")
  ver <- package.description("geoR", lib = lib, field="Version")
  dat <- package.description("geoR", lib = lib, field="Date")
  cat(paste("geoR version ", ver," (", dat, ") is now loaded\n", sep=""))
  cat("-------------------------------------------------\n")
  cat("\n")
  return(invisible(0))
}

