\name{wrappers}

\alias{wrappers}
%\alias{distdiag}
\alias{diffpairs}
\alias{loccoords}
\alias{.diagquadraticformXAX}
\alias{.bilinearformXAY}
\alias{.corr.diaglowertri}
\alias{.Ccor.spatial}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Wrappers for the C functions used in geoR}
\description{
  These functions are \emph{wrappers} for some (but not all)
  the C functions
  included in the \pkg{geoR} package.\cr
  Typically the C code is directly  called from the \pkg{geoR}
  functions but these functions allows independent calls.
}

\usage{
%distdiag(coords)
diffpairs(coords, data)
loccoords(coords, locations)
.diagquadraticformXAX(X, lowerA, diagA) 
.bilinearformXAY(X, lowerA, diagA, Y) 
.corr.diaglowertri(coords, cov.model, phi, kappa)
.Ccor.spatial(x, phi, kappa, cov.model)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix with the data
    coordinates.  }
  \item{data}{an vector with the data values.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix with the
    coordinates of the prediction locations. }
  \item{lowerA}{a vector with the diagonal terms of the symmetric matrix 
    A. }
  \item{diagA}{a vector with the diagonal terms of the symmetric matrix A. }
  \item{X}{a matrix with conforming dimensions. }
  \item{Y}{a matrix with conforming dimensions. }
  \item{cov.model}{covariance model, see \code{\link{cov.spatial}} for options and
    more details.  }
  \item{phi}{numerical value of the correlation function parameter phi.}
  \item{kappa}{numerical value of the correlation function parameter
    kappa.}
  \item{x}{a vector of distances.  }
}
%\details{
%}
\value{
  The outputs for the different functions are: 
  
  \item{diffpairs}{returns a list with elements \code{dist} - the
    distance between pairs of points, and \code{diff} - the difference
    between the values of the attributes.  }
  \item{loccoords}{returns a \eqn{n \times N}{n x N}  matrix with distances
    between data points and prediction locations.  }
  %    \item{distdiag}{returns a vector with distances between data
  %      locations, including the diagonal zero values.  }
  \item{diagquadraticformXAX}{returns a vector with the diagonal term of the
    quadratic form \eqn{X' A X}{t(X) \%*\% A \%*\% X}.  }
  \item{bilinearformXAY}{returns a vector or a matrix with the terms of the
    quadratic form \eqn{X' A Y}{t(X) \%*\% A \%*\% Y}.  }
  \item{corr.diaglowertri}{returns the lower triangle of the
    correlation matrix, including the diagonal.  }
  \item{Ccor.spatial}{returns a vector of values of spatial
    correlations.  }
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


%\note{ ~~further notes~~ }


%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

%\examples{
%}

\keyword{spatial} 
\keyword{programming}
\keyword{interface}

