\name{polygrid}

\alias{polygrid}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Coordinates of Points Inside a Polygon }

\description{
  This function builds a rectangular grid and extracts points which are
  inside of an internal polygonal region.
}

\usage{
polygrid(xgrid, ygrid, poly, vec.inout = FALSE)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{xgrid}{grid values in the \emph{x}-direction.  }
  \item{ygrid}{grid values in the \emph{y}-direction.  }
  \item{poly}{a matrix with the polygon coordinates.  }
  \item{vec.inout}{logical. If \code{TRUE} a logical vector is included
    in the output indicating whether each point of the grid is inside
    the polygon.  Defaults to \code{FALSE}.  }
}

\details{
  This function requires the package \code{splancs}.

  The function works as follows:
  First it creates a grid using the \R function
  \code{\link{expand.grid}} and then it uses the function
  \code{\link[splancs]{inout}} from the package \code{splancs} to extract the points
  of the grid which are inside the polygon.
  
  Within the package \bold{geoR}
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
    the rectangular grid is inside the polygon. Only returned if \code{vec.inout =
      TRUE}.  }
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{expand.grid}}, \code{\link[splancs]{inout}. }
}

\examples{
poly <- matrix(c(.2, .8, .7, .1, .2, .1, .2, .7, .7, .1), ncol=2)
plot(0:1, 0:1, type="n")
lines(poly)
poly.in <- polygrid(seq(0,1,l=11), seq(0,1,l=11), poly, vec=T)
points(poly.in$xy)
}
\keyword{spatial}%-- one or more ...
