\name{nearloc}
\alias{nearloc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Near location to a point}
\description{
  For a given set of points and locations identified by 2D coordinates this
  function finds the nearest location of each point
}
\usage{
nearloc(points, locations, positions = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{points}{a matrix, data-frame or list with the 2D coordinates of
    a set of points for which you want to find the nearest location.  }
  \item{locations}{a matrix, data-frame or list with the 2D coordinates of
    a set of locations.  }
  \item{positions}{logical defining what to be returned. If \code{TRUE}
    the function returns the positions of the locations, otherwise the
    coordinates of the locations. Defaults to \code{FALSE}.  }
}
%details{
% ~~ If necessary, more details than the __description__  above ~~
%
\value{
  If \code{positions = FALSE} the function returns 
  a matrix, data-frame or list of the same type and size
  as the object provided in the argument \code{points} with the
  coordinates of the nearest locations.

  If \code{positions = FALSE} the function returns a vector with the
  position of the nearest points in the \code{locations} object.
}
%references{ ~put references to the literature/web site here ~ }
\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ }

\seealso{\code{\link{loccoords}}}
\examples{
set.seed(276)
gr <- expand.grid(seq(0,1, l=11), seq(0,1, l=11))
plot(gr, asp=1)
pts <- matrix(runif(10), nc=2)
points(pts, pch=19)
near <- nearloc(points=pts, locations=gr)
points(near, pch=19, col=2)
rownames(near)
nearloc(points=pts, locations=gr, pos=TRUE)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
