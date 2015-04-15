\name{dup.coords}
\alias{dup.coords}
\alias{dup.coords.default}
\alias{dup.coords.geodata}
\alias{duplicated.geodata}

\title{Locates duplicated coordinates}
\description{
  This funtions takes an object with 2-D coordinates and returns the
  positions of the duplicated coordinates. Also sets a method
  for \code{duplicated}
}
\usage{
dup.coords(x, \dots)
\method{dup.coords}{default}(x, \dots)
\method{dup.coords}{geodata}(x, incomparables, \dots)
\method{duplicated}{geodata}(x, incomparables, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a two column numeric matrix or data frame.}
  \item{incomparables}{unused. Just for compatibility with
    the generic function \code{\link{duplicated}}.  }
  \item{\dots}{arguments passed to \code{\link{sapply}}. If \code{simplify = TRUE} (default) results
  are returned as an array if possible (when the number of replicates are the same at each replicated location) }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  Function and methods returns \code{NULL} if there are no duplicates
  locations.
  
  Otherwise, the default method returns a
  list where each component is a vector with the positions  or the
  rownames, if available, of the duplicates coordinates.

  The method for \code{geodata} returns a data-frame with
  rownames equals to the positions of the duplicated coordinates,
  the first column is a factor indicating duplicates and
  the remaning are output of \code{\link{as.data.frame.geodata}}.
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br} \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{as.geodata}} for the definition of geodata class,
  \code{\link{duplicated}} for the base function to identify duplicated
  values and \code{\link{jitterDupCoords}} for a function which jitters
  duplicated coordinates.  }

\examples{
## simulating data
dt <- grf(30, cov.p=c(1, .3)) 
## "forcing" some duplicated locations
dt$coords[4,] <- dt$coords[14,] <- dt$coords[24,] <- dt$coords[2,]
dt$coords[17,] <- dt$coords[23,] <- dt$coords[8,]
## output of the method for geodata
dup.coords(dt)
## which is the same as a method for duplicated()
duplicated(dt)
## the default method:
dup.coords(dt$coords)
}

\keyword{spatial}
\keyword{manip}
%\keyword{logic}
