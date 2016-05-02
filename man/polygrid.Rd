\name{polygrid}

\alias{polygrid}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Coordinates of Points Inside a Polygon }

\description{
  This function builds a rectangular grid and extracts points which are
  inside of an internal polygonal region.
}

\usage{
polygrid(xgrid, ygrid, borders, vec.inout = FALSE, ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{xgrid}{grid values in the \emph{x}-direction.  }
  \item{ygrid}{grid values in the \emph{y}-direction.  }
  \item{borders}{a matrix with polygon coordinates defining the
    borders of the region.  }
  \item{vec.inout}{logical. If \code{TRUE} a logical vector is included
    in the output indicating whether each point of the grid is inside
    the polygon.  Defaults to \code{FALSE}.  }
  \item{...}{currently not used (kept for back compatibility). }
}

\details{
  The function works as follows:
  First it creates a grid using the \R{} function
  \code{\link[base]{expand.grid}} and then it uses the geoR'
  internal function
  \code{.geoR_inout()} which wraps usage of \code{\link[sp]{SpatialPoints}}
  and \code{\link[sp]{over}}
  from the package \pkg{sp} to extract the points
  of the grid which are inside the polygon.
  
  Within the package \pkg{geoR}
  this function is typically used to select points in a non-rectangular
  region to perform spatial prediction
  using \code{\link{krige.bayes}}, \code{\link{krige.conv}} or
  \code{\link{ksline}}. It is also useful to produce
  image or perspective plots of the prediction results.
}

\value{
  A list with components:
  \item{xypoly}{an \eqn{n \times 2}{n x 2} matrix with the coordinates of the points inside the polygon.  }
  \item{vec.inout}{logical, a vector indicating whether each point of
    the rectangular grid is inside the polygon. Only returned if \code{vec.inout = TRUE}.  }
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{pred_grid}}, \code{\link{expand.grid}}, \code{\link[sp]{over}},
    \code{\link[sp]{SpatialPoints}}. 
}

\examples{
 poly <- matrix(c(.2, .8, .7, .1, .2, .1, .2, .7, .7, .1), ncol=2)
 plot(0:1, 0:1, type="n")
 lines(poly)
 poly.in <- polygrid(seq(0,1,l=11), seq(0,1,l=11), poly, vec=TRUE)
 points(poly.in$xy)
}
\keyword{spatial}%-- one or more ...
