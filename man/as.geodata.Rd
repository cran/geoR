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
           covar.col = NULL, covar.names = "obj.names",
           realisations = NULL,
           na.action = c("ifany", "ifdata", "ifcovar", "none"),
           rep.data.action, rep.covar.action)
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
  \item{realisations}{optional. A vector indicating the realisation
    number or a number indicating a column in \code{obj} with the
    realisation indicator variable.
    See \code{DETAILS} below.  }
  \item{na.action}{string defining action to be taken in the presence of
    \code{NA}'s. The default option \code{"ifany"} excludes all points
    for which there are \code{NA}'s in the data or covariates. The option \code{"ifdata"} excludes points
    for which there are \code{NA}'s in the data. The default option \code{"ifcovar"} excludes all points
    for which there are \code{NA}'s in the covariates. The option
    \code{"none"} do not exclude points.  }
  \item{rep.data.action}{a string or a function. Defines action to be taken when there is more than
    one data at the same location. The default option \code{"none"}
    keeps the repeated locations, if any. The option \code{"first"} retains
    only the first data recorded at each location. Alternativelly a
    function can be passed and it will be used. For intance if
    \code{mean} is provided, the function will compute and return the average of
    the data at coincident locations.
    The non-default
    options will eliminate the repeated locations.  }
  \item{rep.covar.action}{idem to \code{rep.data.locations}, to be applied to
    the covariates, if any. Defaults to
    the same option set for \code{rep.data.locations}. }
}

\details{
  Objects of the class \code{"geodata"} contain data for
  geostatistical analysis using the package \bold{geoR}.
  Storing data in this format facilitates the usage of the functions in \bold{geoR}.
  However, conversion of objects to this class is not obligatory
  to carry out the analysis.

  \bold{Realisations}\cr
  Tipically geostatistical data correspond to a unique realisation of
  the spatial process.
  However, sometimes different "realisations" are possible.
  For instance, if data are collected in the same area at different
  points in time and independence between time points is assumed,
  each time can be considered a different "replicate" or "realisation"
  of the same process. The argument \code{realisations} takes a vector
  indication the replication number and can be passed to other \code{geoR}
  functions as, for instance, \code{\link{likfit}}.

  The data format is similar to the usual \code{geodata} format in
  \code{geoR}.
  Suppose there are realisations (times) \eqn{1, \ldots, J}{1, ..., J}
  and for each realisations \eqn{n_1, ..., n_j} observations are available.
  The coordinates for different realisations
  should be combined in a single \eqn{n \times 2}{nx2} object,
  where \eqn{n=n_1 + \ldots + n_J}{n=n_1 + ... + n_J}.
  Similarly for the data vector and covariates (if any).
}

\value{
  An object of the \code{\link{class}} \code{"geodata"} which is a list
  with two obligatory components (coords and data)
  and other optional components:
  \item{coords}{an \eqn{n \times 2}  matrix where \eqn{n} is the number
    of spatial locations.}
  \item{data}{a vector of length \eqn{n}, for the univariate case or,
    an \eqn{n \times v} matrix or data-frame for the multivariate case,
    where \eqn{v} is the number of variables.}
  \item{covariates}{a vector of length \eqn{n} or
    an \eqn{n \times p} matrix with covariate(s) values,
    where \eqn{p} is the number of covariates. Only returned if
    covariates are provided.  }
  \item{realisations}{a vector on size \eqn{n} with the replication
    number. Only returned if argument realisations is provided.  } 
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
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



