\name{as.geodata}
\alias{as.geodata}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Converts an Object to the Class "geodata"}

\description{
  Converts a matrix or a data-frame to an object of the
  \code{\link{class}}  \code{"geodata"}.
  Objects of the class \code{"geodata"} are lists with two obligatory
  components, \code{coords} and \code{data}.  Optional components are
  allowed and a typical example is a vector or matrix with  values
  of the covariate(s).
}

\usage{
as.geodata(obj, coords.col = 1:2, data.col = 3, data.names = NULL, 
           covar.col = NULL, covar.names = "obj.names")
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{a matrix or data-frame. Each line corresponds to one
    spatial location. It should contain values of 
    coordinates, data and, optionally, covariates at the locations.  }
  \item{coords.col}{a vector with the column numbers corresponding to the
    spatial coordinates.  }
  \item{data.col}{a scalar or vector with column number(s) corresponding
    to the data.  }
  \item{data.names}{a string or vector of strings with names for
    the data columns. Only valid if there is more than one column of
    data.
    By default, the names in the original object are used.  }
  \item{covar.col}{optional. A scalar or vector
    of the column number(s) corresponding to the covariate(s).  }
  \item{covar.names}{a string or vector of strings with the name(s) of
    the covariates. By default, the names in the original object are used.  }
}

\details{
  Objects of the class \code{"geodata"} contain data for
  geostatistical analysis using the package \bold{geoR}.
  Storing data in this format facilitates the usage of the functions in \bold{geoR}.
  However, conversion of objects to this class is not obligatory
  to carry out the analysis.
  
}

\value{
  An object of the \code{\link{class}} \code{"geodata"} which is a list
  with two obligatory and one optional component:
  \item{coords}{an \eqn{n \times 2}  matrix where \eqn{n} is the number
    of spatial locations.}
  \item{data}{a vector of length \eqn{n}, for the univariate case or,
  an \eqn{n \times v} matrix or data-frame for the multivariate case,
  where \eqn{v} is the number of variables.}
  \item{covariate}{a vector of length \eqn{n} or
    an \eqn{n \times p} matrix with covariate(s) values,
    where \eqn{p} is the number of covariates.}

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
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{read.geodata}} for reading data from an
  \emph{ASCII} file and \code{\link{list}} for general information on lists.
}

%\examples{}
\keyword{spatial}
\keyword{classes}
\keyword{manip}



