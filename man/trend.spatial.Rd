\name{trend.spatial}

\alias{trend.spatial}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Prepare Trend Matrix}

\description{
  Builds the trend matrix according to the specification
  of the mean part of the model. 
}

\usage{
trend.spatial(trend, coords = NULL)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{trend}{specifies the mean part of the model. The options are:
    \code{"cte"} (constant mean - default option),
    \code{"1st"} (a first degree polynomial
    on the coordinates), \code{"2nd"} (a second degree polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend). }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing in each row Euclidean
    coordinates. Needs to be provided only if \code{trend = "1st"} or
    \code{trend = "2nd"}.  }
}

%\details{
%}

\value{
  An \eqn{n \times p}{n x p} trend matrix where \eqn{n}
  is the number of spatial
  locations and \eqn{p} is the number of mean parameters in the model. 
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\note{This is an auxiliary function called by
  other \bold{geoR} functions.  }

%\seealso{}
%
%\examples{
%
%}

\keyword{spatial}%-- one or more ...

