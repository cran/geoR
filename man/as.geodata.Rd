\name{as.geodata}
\alias{as.geodata}
\alias{geodata}
\alias{as.geodata.default}
\alias{as.geodata.SpatialPointsDataFrame}
\alias{as.data.frame.geodata}
\alias{as.geodata.geodata.frame}
\alias{is.geodata}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Converts an Object to the Class "geodata"}

\description{
  The default method converts a matrix or a data-frame
  to an object of the
  \code{\link{class}}  \code{"geodata"}.\cr
  Objects of the class \code{"geodata"} are lists with two obligatory
  components: \code{coords} and \code{data}.
  Optional components are allowed and a typical example is a vector or
  matrix with covariate(s) values.
}

\usage{
as.geodata(obj, \dots)

\method{as.geodata}{default}(obj, coords.col = 1:2, data.col = 3, data.names = NULL, 
                   covar.col = NULL, covar.names = "obj.names",
                   units.m.col = NULL, realisations = NULL,
                   na.action = c("ifany", "ifdata", "ifcovar", "none"),
                   rep.data.action, rep.covar.action, rep.units.action,
                   ...)

\method{as.data.frame}{geodata}(x, \dots, borders = TRUE)

\method{as.geodata}{geodata.frame}(obj, \dots)

\method{as.geodata}{SpatialPointsDataFrame}(obj, data.col = 1, \dots)

is.geodata(x)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{a matrix or data-frame where each line corresponds to one
    spatial location. It should contain values of 2D coordinates,
    data and, optionally, covariate(s) value(s) at the locations.
    A method for \code{\link[sp]{SpatialPointsDataFrame}} is also provided.
    It can also take an output of the function \code{\link{grf}}, see
    DETAILS below. }
  \item{coords.col}{a vector with the column numbers corresponding to the
    spatial coordinates.  }
  \item{data.col}{a scalar or vector with column number(s) corresponding
    to the data.  }
  \item{data.names}{optional. A string or vector of strings with names for
    the data columns. Only valid if there is more than one column of
    data. By default, takes the names from the original object.  }
  \item{covar.col}{optional. A scalar or numeric vector
    with the column number(s) corresponding to the covariate(s).
    Alternativelly can be a character vector with the names of the
    covariates. }
  \item{covar.names}{optional. A string or vector of strings with the name(s) of
    the covariates. By default take the names from the original object.  }
  \item{units.m.col}{optional. A scalar 
    with the column number corresponding to the offset variable.
    Alternativelly can be a character vector with the name of the
    offset. This option is particularly relevant when using the package
    \pkg{geoRglm}. All values must be greater then zero. }
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
    function can be passed and it will be used. For instance if
    \code{mean} is provided, the function will compute and return the average of
    the data at coincident locations. The non-default
    options will eliminate the repeated locations.  }
  \item{rep.covar.action}{idem to \code{rep.data.locations}, to be applied to
    the covariates, if any. Defaults to
    the same option set for \code{rep.data.locations}. }
  \item{x}{an object which is tested for the class \code{geodata}.}
  \item{rep.units.action}{a string or a function.
    Defines action to be taken on the element \code{units.m}, if present
    when there is more than
    one data at the same location. The default option is the same
    value set for \code{rep.data.action}.  }
  \item{borders}{logical. If TRUE the element borders in the
    \code{geodata}
  object is set as an attribute of the data-frame.  }
  \item{\dots}{values to be passed for the methods.  }
}

\details{
  Objects of the class \code{"geodata"} contain data for
  geostatistical analysis using the package \pkg{geoR}.
  Storing data in this format facilitates the usage of the functions in \pkg{geoR}.
  However, conversion of objects to this class is not obligatory
  to carry out the analysis.
  
  \code{NA}'s are not allowed in the coordinates. By default the
  respective rows will not be included in the output.

  \bold{Realisations}\cr
  Tipically geostatistical data correspond to a unique realisation of
  the spatial process.
  However, sometimes different "realisations" are possible.
  For instance, if data are collected in the same area at different
  points in time and independence between time points is assumed,
  each time can be considered a different "replicate" or "realisation"
  of the same process. The argument \code{realisations} takes a vector
  indication the replication number and can be passed to other \pkg{geoR}
  functions as, for instance, \code{\link{likfit}}.

  The data format is similar to the usual \code{geodata} format in
  \pkg{geoR}.
  Suppose there are realisations (times) \eqn{1, \ldots, J}{1, ..., J}
  and for each realisations \eqn{n_1, ..., n_j} observations are available.
  The coordinates for different realisations
  should be combined in a single \eqn{n \times 2}{nx2} object,
  where \eqn{n=n_1 + \ldots + n_J}{n=n_1 + ... + n_J}.
  Similarly for the data vector and covariates (if any).

  \bold{grf objects}\cr
  If an object of the class \code{grf} is provided the functions just
  extracts the elements \code{coords} and \code{data} of this object. 
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
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{read.geodata}} for reading data from an
  \emph{ASCII} file and \code{\link{list}} for general information on lists.
}

\examples{
\dontrun{
## converting the data-set "topo" from the package MASS (VR's bundle)
## to the geodata format:
if(require(MASS)){
topo
topogeo <- as.geodata(topo)
names(topogeo)
topogeo
}
}
}
\keyword{spatial}
\keyword{classes}
\keyword{manip}
