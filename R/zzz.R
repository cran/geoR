".First.lib" <-
  function(lib, pkg)
{
  library.dynam("geoR", package = pkg, lib.loc = lib)
  cat("\n")
  cat("-------------------------------------------------\n")
  ## from 1.9-0, package.description is deprecated in favour of
  ## packageDescription (which doesn't exist in previous versions)
  if(!exists("packageDescription",mode="function")){
    pkg.info <- package.description("geoR", lib.loc = lib,
                                    fields=c("Title","Version","Date"))
    pkg.info <- list(Title=pkg.info[1], Version=pkg.info[2],
                     Date=pkg.info[3])
  }
  else pkg.info <- packageDescription("geoR", lib.loc = lib,
                                      fields=c("Title","Version","Date"))
  cat(pkg.info$Title)
  cat("\n")
  cat(paste("geoR version ", pkg.info$Version, " (", pkg.info$Date, ") is now loaded\n", sep=""))
  cat("-------------------------------------------------\n")
  cat("\n")
  return(invisible(0))
}

