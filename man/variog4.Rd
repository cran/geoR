\name{variog4}
\alias{variog4}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computes Directional Variograms}
\description{
  Computes directional variograms for 4 directions provided by the user.
}
\usage{
variog4(geodata, coords = geodata$coords, data = geodata$data, 
        uvec = "default", breaks = "default", trend = "cte", lambda = 1,
        option = c("bin", "cloud", "smooth"),
        estimator.type = c("classical", "modulus"), 
        nugget.tolerance, max.dist, pairs.min = 2,
        bin.cloud = FALSE, direction = c(0, pi/4, pi/2, 3*pi/4), tolerance = pi/8,
        unit.angle = c("radians", "degrees"), messages, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{geodata}{a list containing element \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \pkg{geoR} data-set.
    If not provided the arguments
    \code{coords} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing
    coordinates of the \eqn{n} data locations in each row.
    Defaults to \code{geodata$coords}, if provided.}
  \item{data}{a vector or matrix with data values.
    If a matrix is provided, each column is regarded as one variable or realization.
    Defaults to \code{geodata$data}, if provided.}
  \item{uvec}{a vector with values to define the variogram binning. For
    further details see documentation for \code{\link{variog}}.  } 
  \item{breaks}{a vector with values to define the variogram binning. For
    further details see documentation for \code{\link{variog}}.  } 
  \item{trend}{specifies the mean part of the model.
    The options are:
    \code{"cte"} (constant mean),
    \code{"1st"} (a first order polynomial
    on the coordinates), \code{"2nd"} (a second order polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend). Defaults to \code{"cte"}.  }
  \item{lambda}{values of the Box-Cox transformation parameter.
    Defaults to \eqn{1} (no
    transformation). If another value is provided the variogram is
    computed after transforming the 
    data. A case of particular interest is \eqn{\lambda = 0}{lambda = 0}
    which corresponds to log-transformation.  }
  \item{option}{defines the output type: the options \code{"bin"} returns values of
    binned variogram, \code{"cloud"} returns the variogram cloud and
    \code{"smooth"} returns the kernel smoothed variogram.
    Defaults to \code{"bin"}.}
  \item{estimator.type}{\code{"classical"} computes the classical method of
    moments estimator.  \code{"modulus"} returns the variogram
    estimator suggested by Hawkins and Cressie (see Cressie, 1993, pg 75).
    Defaults to \code{"classical"}.  }
  \item{nugget.tolerance}{a numeric value. Points which are separated
    by a distance less than this value
    are considered co-located. Defaults to zero.  }
  \item{max.dist}{a numerical value defining the maximum distance for
    the variogram. Pairs of locations
    separated for distance larger than this value are ignored for the
    variogram calculation. Defaults to the maximum distance among the pairs
    of data locations.  }
  \item{pairs.min}{a integer number defining the minimum numbers of
    pairs for the bins.
    For \code{option = "bin"},
    bins with number of pairs smaller than this
    value are ignored. Defaults to \code{NULL}.  }
  \item{bin.cloud}{logical. If \code{TRUE} and
    \code{option = "bin"} the cloud values for each class are
    included in the output. Defaults to \code{FALSE}.  }
  \item{direction}{a vector with values of 4 angles, indicating the
    directions for which the variograms will be computed. Default
    corresponds to \code{c(0, 45, 90, 135)} (degrees). }
  \item{tolerance}{numerical value for the tolerance angle, when
    computing directional variograms. The value must be in the interval
    \eqn{[0, 90]} degrees.  Defaults to \eqn{\pi/8}{pi/8}.  }
  \item{unit.angle}{defines the unit for the specification of angles in
    the two previous arguments. Options are \code{"degrees"} and \code{"radians"}. }
  \item{messages}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{arguments to be passed to the function \code{\link{ksmooth}}, if
    \code{option = "smooth"}.  }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  The output is an object of the class \code{variog4},
  a list with five components.
  The first four elements are estimated variograms for the directions
  provided and the last is the omnidirectional variogram. 
  Each individual component is an object of the class \code{variogram},
  an output of the function \code{\link{variog}}.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{variog}} for variogram calculations and
  \code{\link{plot.variog4}} for plotting results}

\examples{
var4 <- variog4(s100, max.dist=1)
plot(var4)
}
\keyword{spatial}%-- one or more ...
