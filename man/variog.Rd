\name{variog}

\alias{variog}
\alias{.rfm.bin}
\alias{.define.bins}

\concept{variogram}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Compute Empirical Variograms}

\description{
  Computes sample (empirical) variograms with options for  the \emph{classical} or \emph{robust}
  estimators. Output can be returned as a \code{binned variogram}, a \code{
    variogram cloud} or a \code{smoothed variogram}. Data
  transformation (Box-Cox) is allowed.
  ``Trends'' can be specified and are fitted by ordinary least
  squares in which case the variograms are computed using the
  residuals.
}

\usage{
variog(geodata, coords = geodata$coords, data = geodata$data, 
       uvec = "default", breaks = "default",
       trend = "cte", lambda = 1,
       option = c("bin", "cloud", "smooth"),
       estimator.type = c("classical", "modulus"), 
       nugget.tolerance, max.dist, pairs.min = 2,
       bin.cloud = FALSE, direction = "omnidirectional", tolerance = pi/8,
       unit.angle = c("radians","degrees"), angles = FALSE, messages, \dots) 
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
  \item{uvec}{a vector with values used to define the variogram binning. Only
    used when \code{option = "bin"}.
    See \code{DETAILS} below for more details on how to specify the bins.
  }
  \item{breaks}{a vector with values to define the variogram binning. Only
    used when \code{option = "bin"}.
    See \code{DETAILS} below for more details on how to specify the bins.
  }
  \item{trend}{specifies the mean part of the model.
    See documentation of \code{\link{trend.spatial}} for further details.
    Defaults to \code{"cte"}.  }
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
    variogram calculation. If not provided defaults takes the maximum
    distance among all pairs of data locations.  }
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
    \eqn{[0, \pi]}{[0, pi]} radians (\eqn{[0, 180]} degrees).  }
  \item{tolerance}{numerical value for the tolerance angle, when
    computing directional variograms. The value must be in the interval
    \eqn{[0, \pi/2]}{[0, pi/2]} radians (\eqn{[0, 90]} degrees).  Defaults to \eqn{\pi/8}{pi/8}.  }
  \item{unit.angle}{defines the unit for the specification of angles in
    the two previous arguments. Options are \code{"radians"} and
    \code{"degrees"}, with default to \code{"radians"}. }
  \item{angles}{Logical with default to \code{FALSE}. If \code{TRUE} the
    function also returns the angles between the pairs of points
    (unimplemented). }
  \item{messages}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{arguments to be passed to the function \code{\link{ksmooth}}, if
    \code{option = "smooth"}.  }
}

\details{
  Variograms are widely used in geostatistical analysis for exploratory
  purposes, to estimate covariance parameters and/or to compare theoretical
  and fitted models against sample variograms.

  \bold{Estimators}\cr
  The two estimators currently implemented are:
  \itemize{
    \item \emph{classical} (method of moments) estimator:
    \deqn{\gamma(h) = \frac{1}{2N_h} \sum_{i=1}^{N_h}[Y(x_{i+h}) -
      Y(x_i)]^2}{gamma(h) = (1/2N_h) sum [Y(x_i+h) - Y(x_i)]^2}

    \item Hawkins and Cressie's \emph{modulus} estimator
    \deqn{\gamma(h) = \frac{[\frac{1}{N_h}\sum_{i=1}^{N_h}|Y(x_{i+h}) -
	Y(x_i)|^{\frac{1}{2}}]^4}{0.914 + \frac{0.988}{N_h}}}{\gamma(h)
      = ([(1/N_h)\sum_{i=1}^{N_h} |Y(x_{i+h}) -
	Y(x_i)|^(1/2)]^4)/(0.914 + (0.988/N_h))}
      }

  \bold{Defining the bins}

  \emph{The defaults}\cr
  If arguments \code{breaks} and \code{uvec} are not provided, the
  binning is defined as follows:\cr
  \enumerate{
    \item read the argument \code{max.dist}. If not provided it is set
    to the maximum distance between the pairs of points.
    \item the center of the bins are initially defined by the sequence \code{u = seq(0,
      max.dist, l = 13)}.
    \item the interval spanned by each bin is given by the mid-points
    between the centers of the bins.
%    \item The center of the first bin is set to \code{u[1] = u[2]/4}.
  }

  If an vector is passed to the argument \code{breaks} its elements are
  taken as the limits of the bins (classes of distance) and the argument \code{uvec} is ignored.

  \emph{Variations on the default}\cr
  The default definition of the bins is different for some particular
  cases.
  \enumerate{
    \item if there are  coincident data locations
    the bins follows the default above but
    one more bin is added at the origin (distance zero) for
    the collocated points.
    
    \item if the argument \code{nugget.tolerance} is provided the
    separation distance between all pairs
    in the interval \eqn{[0, nugget.tolerance]} are set to zero.
    The first bin distance is set to zero (\code{u[1] = 0}).
    The remaining bins follows the default.

    \item if a scalar is provided to the argument \code{uvec} the
    default number of bins is defined by this number.
%    \code{u = seq(0, max.dist, l = uvec)}

    \item if a vector is provided to the argument \code{uvec},
    its elements are taken as central points of the bins.
  }
  
  %    The values of \code{uvec} defines the mid-points of the bins.\cr
  %    If \eqn{uvec[1] > 0} the first bin is: \eqn{0 < u <= uvec[2] - 0.5*(uvec[2] - uvec[1])}.\cr
  %    If \eqn{uvec[1] = 0} first bin is: \eqn{0 < u <= 0.5*uvec[1]} and
  %    \eqn{uvec[1]} is replaced by the midpoint of this interval.
  
}

\value{
  An object of the \code{\link{class}} \code{variogram} which is a
  list with the following components:
  \item{u }{a vector with distances.  }
  \item{v }{a vector with estimated variogram values at distances given
    in \code{u}.  }
  \item{n }{number of pairs in each bin, if
    \code{option = "bin"}.  }
  \item{sd }{standard deviation of the values in each bin.  }
  \item{bins.lim }{limits defining the interval spanned by each bin.  }
  \item{ind.bin }{a logical vector indicating whether the number of
    pairs in each bin is greater or equal to the value in the argument
    \code{pairs.min}.  }
  \item{var.mark }{variance of the data.  }
  \item{beta.ols }{parameters of the mean part of the model fitted by
    ordinary least squares.  }
  \item{output.type }{echoes the \code{option} argument.  }
  \item{max.dist}{maximum distance between pairs allowed in the
    variogram calculations.}
  \item{estimator.type }{echoes the type of estimator used.  }
  \item{n.data }{number of data.  }
  \item{lambda}{value of the transformation parameter.  }
  \item{trend}{trend specification.  }
  \item{nugget.tolerance}{value of the nugget tolerance argument.  }
  \item{direction }{direction for which the variogram was computed.  }
  \item{tolerance }{tolerance angle for directional variogram.  }
  \item{uvec }{lags provided in the function call.  }
  \item{call }{the function call.  }
}

\references{
  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.
  
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
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
legend("topleft", legend=c("omnidirectional", expression(0 * degree)), lty=c(2,1))
}

\keyword{spatial}
\keyword{smooth}
\keyword{robust}





