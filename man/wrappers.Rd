\name{wrappers}

\alias{wrappers}
\alias{distdiag}
\alias{loccoords}
\alias{diagquadraticformXAX}
\alias{bilinearformXAY}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Wrappers for the C functions used in geoR}
\description{
  These functions are \emph{wrappers} for some, but no tall,
  the C functions
  included in the \bold{geoR} package.\cr
  Typically the C code is directly  called from the \bold{geoR}
  functions but these functions allows independent calls.
}

\usage{
distdiag(coords)
loccoords(coords, locations)
diagquadraticformXAX(X, lowerA, diagA) 
bilinearform(X, lowerA, diagA, Y) 
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix with the data
    coordinates.  }
  \item{locations}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates of the prediction locations. }
  \item{lowerA}{a vector with the diagonal terms of the symmetric matrix 
    A. }
  \item{diagA}{a vector with the diagonal terms of the symmetric matrix A. }
  \item{X}{a matrix with conforming dimensions. }
  \item{Y}{a matrix with conforming dimensions. }
}
%\details{
%}
\value{
  The outputs are: 
  
  \itemize{
    \item{loccoords}{returns a vector with distances between data points and
      prediction locations.  }
    \item{distdiag}{returns a vector with distances between data
      locations, including 
      the diagonal zero values.  }
    \item{diagquadraticformXAX}{returns a vector with the diagonal term of the
      quadratic form \eqn{X' A X}{t(X) \%*\% A \%*\% X}.  }
    \item{bilinearform}{returns a vector or a matrix with the terms of the
      quadratic form \eqn{X' A Y}{t(X) \%*\% A \%*\% Y}.  }
  }
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


%\note{ ~~further notes~~ }


%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

%\examples{
%}

\keyword{spatial} 
\keyword{programming}
\keyword{interface}
