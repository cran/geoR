\name{variog}

\alias{variog}
\alias{rfm.bin}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Compute Empirical Variograms}

\description{
  Computes sample (empirical) variograms with options for  the \emph{classical} or \emph{robust}
  estimator. Output can be returned as a \code{binned variogram}, a \code{
    variogram cloud} or a \code{smoothed variogram}. Data
  transformation (Box-Cox) is allowed. Trends fitted by ordinary least
  squares can be removed. In this case variograms are computed using the
  residuals.
}

\usage{
variog(geodata, coords=geodata$coords, data=geodata$data, 
       uvec = "default", trend = "cte", lambda = 1, 
       option = c("bin", "cloud", "smooth"), 
       estimator.type = c("classical", "modulus"), 
       nugget.tolerance = 0, max.dist = NULL, pairs.min = 2, 
       bin.cloud = FALSE, direction = "omnidirectional",
       tolerance = pi/8, unit.angle = c("radians", "degrees"),
       messages.screen = TRUE, ...) 
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing element \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data-set.
    If not provided the arguments
    \code{coords} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing
    coordinates of the \eqn{n} data locations in each row.
    Defaults to \code{geodata$coords}, if provided.}
  \item{data}{a vector or matrix with data values.
    If a matrix is provided, each column is regarded as one variable or realization.
    Defaults to \code{geodata$data}, if provided.}
  \item{uvec}{a vector with values defining the variogram binning. Only
    used when
    \code{option = "bin"}. The values of \code{uvec} defines the mid-points of the bins.\cr
    If \eqn{uvec[1] > 0} the first bin is: \eqn{0 < u <= uvec[2] - 0.5*(uvec[2] - uvec[1])}.\cr
    If \eqn{uvec[1] = 0} first bin is: \eqn{0 < u <= 0.5*uvec[1]} and
    \eqn{uvec[1]} is replaced by the midpoint of this interval.}
  \item{trend}{specifies the mean part of the model.
    The options are:
    \code{"cte"} (constant mean),
    \code{"1st"} (a first degree polynomial
    on the coordinates), \code{"2nd"} (a second degree polynomial
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
  \item{direction}{a numerical value for the directional (azimuth) angle. This
    used to specify directional variograms. Default defines the
    omnidirectional variogram. The value must be in the interval
    \eqn{[0, 180]} degrees.  }
  \item{tolerance}{numerical value for the tolerance angle, when
    computing directional variograms. The value must be in the interval
    \eqn{[0, 90]} degrees.  Defaults to \eqn{\pi/8}{pi/8}.  }
  \item{unit.angle}{defines the unit for the specification of angles in
    the two previous arguments. Options are \code{"degrees"} and \code{"radians"}. }
  \item{messages.screen}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{arguments to be passed to the function \code{\link{ksmooth}}, if
    \code{option = "smooth"}.  }
}

\details{
  Variograms are widely used in geostatistical analysis for exploratory
  purposes, to estimate covariance parameters and/or to compare theoretical
  and fitted models against sample variograms.

  The two estimators currently implemented are:
  \itemize{
    \item \emph{classical} estimator:
    \deqn{\gamma(h) = \frac{1}{N_h} \sum_{i=1}^{N_h}[Y(x_{(i+h)}) -
      Y(x_i)]^2}{gamma(h) = (1/N_h) sum (Y(x_i+h) - Y(x_i))^2}
    \item Hawkins and Cressie's \emph{modulus} estimator
    \deqn{\gamma(h) = [\frac{1}{N_h} \sum_{i=1}^{N_h}|Y(x_{i+h}) -
      Y(x_i)|^{\frac{1}{2}}]^4}{gamma(h) = ((1/N_h) sum |Y(x_(i+h)) - Y(x_i)|^(1/2))^4}
}
}

\value{
  An object of the \code{\link{class}} \code{variogram} which is a
  list with the following components:
  \item{u }{a vector with distances.  }
  \item{v }{a vector with estimated variogram values at distances given
    in \code{u}.  }
  \item{n }{number of pairs in each bin, if
    \code{option = "bin"}.  }
  \item{var.mark }{variance of the data.  }
  \item{output.type }{echoes the \code{option} argument.  }
  \item{estimator.type }{echoes the type of estimator used.  }
  \item{direction }{direction for which the variogram was computed.  }
  \item{tolerance }{tolerance angle for directional variogram.  }
  \item{uvec }{lags provided in the function call.  }
  \item{call }{the function call.  }
  
}

\references{
  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.
  
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{variog4}} for more on computation of
  directional variograms, 
  \code{\link{variog.model.env}} and \code{\link{variog.mc.env}} for
  variogram envelopes,
  \code{\link{variofit}}  for variogram based
  parameter estimation and
  \code{\link{plot.variogram}} for graphical output.
}

\examples{
# Loading data:
if(is.R()) data(s100)
#
# computing variograms:
#
# binned variogram
vario.b <- variog(s100, max.dist=1)
# variogram cloud
vario.c <- variog(s100, max.dist=1, op="cloud")
#binned variogram and stores the cloud
vario.bc <- variog(s100, max.dist=1, bin.cloud=TRUE)
# smoothed variogram
vario.s <- variog(s100, max.dist=1, op="sm", band=0.2)
#
#
# plotting the variograms:
par(mfrow=c(2,2))
plot(vario.b, main="binned variogram") 
plot(vario.c, main="variogram cloud")
plot(vario.bc, bin.cloud=TRUE, main="clouds for binned variogram")  
plot(vario.s, main="smoothed variogram") 

# computing a directional variogram
vario.0 <- variog(s100, max.dist=1, dir=0, tol=pi/8)
plot(vario.b, type="l", lty=2)
lines(vario.0)
legend(0, 1.2, legend=c("omnidirectional", expression(0 * degree)), lty=c(2,1))
#
}

\keyword{spatial}
\keyword{smooth}





