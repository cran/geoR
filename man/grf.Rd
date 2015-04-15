\name{grf}

\alias{grf}
\alias{geoR2RF}
\alias{.grf.aux1}
\alias{grfclass}
\alias{lines.grf}


\concept{simulation of Gaussian random fields}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Simulation of Gaussian Random Fields}

\description{
  \code{grf()} generates (unconditional)
  simulations of Gaussian random fields for
  given covariance parameters.
  \code{geoR2RF} converts model specification used by \pkg{geoR}
  to the correponding one in \pkg{RandomFields}.
}


\usage{
grf(n, grid = "irreg", nx, ny, xlims = c(0, 1), ylims = c(0, 1),
    borders, nsim = 1, cov.model = "matern",
    cov.pars = stop("missing covariance parameters sigmasq and phi"), 
    kappa = 0.5, nugget = 0, lambda = 1, aniso.pars,
    mean = 0, method, RF=TRUE, messages)

geoR2RF(cov.model, cov.pars, nugget = 0, kappa, aniso.pars)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{n}{number of points (spatial locations) in each simulations.  }
  \item{grid}{optional. An \eqn{n \times 2}{n x 2} matrix with coordinates of the
    simulated data.  }
  \item{nx}{optional. Number of points in the X direction.  }
  \item{ny}{optional. Number of points in the Y direction.  }
  \item{xlims}{optional. Limits of the area in the X direction. Defaults
  to \eqn{[0,1]}.  }
  \item{ylims}{optional. Limits of the area in the Y direction.  Defaults
  to \eqn{[0,1]}.  }
  \item{borders}{optional. Typically a two coluns matrix especifying a
    polygon. See DETAILS below.}
  \item{nsim}{Number of simulations. Defaults to 1.}
  \item{cov.model}{correlation function. See \code{\link{cov.spatial}} for
    further details. Defaults to the 
    \emph{exponential} model.  }
  \item{cov.pars}{a vector with 2 elements or an \eqn{n \times 2}{n x 2}
    matrix with values of  the covariance parameters
    \eqn{\sigma^2}{sigma^2} (partial sill) and \eqn{\phi}{phi} (range
    parameter). If a vector, the elements are the values of
    \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.
    If a matrix, corresponding to a model with several structures, the
    values of \eqn{\sigma^2}{sigma^2} 
    are in the first column and the values of \eqn{\phi}{phi} are in the second.  }
  \item{kappa}{additional smoothness parameter required only for the
    following correlation 
    functions: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} 
     and  \code{"gneiting.matern"}. More details on the documentation for the
    function \code{\link{cov.spatial}}.  }
  \item{nugget}{the value of the nugget effect parameter \eqn{\tau^2}{tau^2}.}
  \item{lambda}{value of the Box-Cox transformation parameter. The value \eqn{\lambda
      = 1}{lambda = 1} corresponds to no transformation, the default.
    For any other value of \eqn{\lambda}{lambda} Gaussian data is
    simulated and then transformed.  }
  \item{aniso.pars}{geometric anisotropy parameters. By default an
    isotropic field is assumed and this argument is ignored.
    If a vector with 2 values is provided, with values for the
    anisotropy angle \eqn{\psi_A}{psi_A} (in
    radians) and 
    anisotropy ratio \eqn{\psi_A}{psi_A}, the coordinates
    are transformed, 
    the simulation is performed on the isotropic (transformed) space
    and then the coordinates are back-transformed such that the resulting
    field is anisotropic. Coordinates transformation is performed 
    by the function \code{\link{coords.aniso}}.  }
  \item{mean}{a numerical vector, scalar or the same length of the
    data to be simulated. Defaults to zero.}
  \item{method}{simulation method with options for
    \code{"cholesky"}, \code{"svd"}, \code{"eigen"}, \code{"RF"}.
    Defaults to the \emph{Cholesky}
    decomposition.  See section \code{DETAILS} below.  }
  \item{RF}{logical, with defaults to TRUE, indicating whether
    the algorithm should try to use the function
    \code{\link[RandomFields]{GaussRF}} from the package
    \pkg{RandomFields} in case of \code{method} is missing
    and the number of points is greater than 500. }
  \item{messages}{logical,  indicating
    whether or not status messages are printed on the screen (or output device)
    while the function is running. Defaults to \code{TRUE}.  } 
}


\details{
  For the methods \code{"cholesky"}, \code{"svd"} and \code{"eigen"} the
  simulation consists of multiplying a vector of standardized
  normal deviates by a square root of the covariance matrix.
  The square root of a matrix is not uniquely defined. These 
  three methods differs in the way they compute the
  square root of the (positive definite) covariance matrix.

  The previously available
  method \code{method = "circular.embedding"} is no longer available
  in \pkg{geoR}. For simulations in a fine grid and/or with a large number
  of points use the package \pkg{RandomFields}.

  The option \code{"RF"} calls internally  the function
  \code{\link[RandomFields]{GaussRF}} from the package \pkg{RandomFields}.

  The argument \code{borders}, if provides takes a
  polygon data set following argument \code{poly}
  for the \pkg{splancs}' function \code{\link[splancs]{csr}}, in case of
  \code{grid="reg"} or \code{\link[splancs]{gridpts}}, in case of
  \code{grid="irreg"}. For the latter the simulation will have 
  \emph{approximately} \dQuote{n} points.

}


\value{
  \code{grf} returns a list with the components:
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix with the coordinates of the
    simulated data.  }
  \item{data}{a vector (if \code{nsim = 1}) or a matrix with the
    simulated values. For the latter each column corresponds to one
    simulation.  }
  \item{cov.model}{a string with the name of the correlation function.  }
  \item{nugget}{the value of the nugget parameter.}
  \item{cov.pars}{a vector with the values of \eqn{\sigma^2}{sigma^2}
    and \eqn{\phi}{phi}, respectively.  }
  \item{kappa}{value of the parameter \eqn{\kappa}{kappa}.  }
  \item{lambda}{value of the Box-Cox transformation parameter
    \eqn{\lambda}{lambda}.  }
  \item{aniso.pars}{a vector with values of the anisotropy parameters, if
    provided in the function call.  }
  \item{method}{a string with the name of the simulation method used.}
  \item{sim.dim}{a string "1d" or "2d" indicating the spatial dimension of the
    simulation.}
  \item{.Random.seed}{the random seed by the time the function was
    called.  }
  \item{messages}{messages produced by the function describing the
    simulation.  }
  \item{call}{the function call.  }

  \code{geoR2grf} returns a list with the components:
  \item{model}{\pkg{RandomFields} name of the correlation model}
  \item{param}{\pkg{RandomFields} parameter vector}
}


\references{
  Wood, A.T.A. and Chan, G. (1994) Simulation of stationary Gaussian
  process in \eqn{[0,1]^d}{[0,1]^d}.
  \emph{Journal of Computatinal and Graphical Statistics}, \bold{3}, 409--432.
  
  Schlather, M. (1999) \emph{Introduction to positive definite functions
    and to unconditional simulation of random fields}. Tech. Report
  ST--99--10, Dept Maths and Stats, Lancaster University.

  Schlather, M. (2001) \emph{Simulation and Analysis of Random Fields}. R-News \bold{1} (2), p. 18-20.
  
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}


\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br},  \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


\seealso{\code{\link{plot.grf}} and \code{\link{image.grf}} for
  graphical output,
  \code{\link{coords.aniso}} for anisotropy coordinates transformation
  and, \code{\link{chol}},
  \code{\link{svd}} and \code{\link{eigen}} for methods of matrix
  decomposition and \code{\link[RandomFields]{GaussRF}} function
  in the package \pkg{RandomFields}.
}

\examples{
sim1 <- grf(100, cov.pars = c(1, .25))
# a display of simulated locations and values
points(sim1)   
# empirical and theoretical variograms
plot(sim1)
## alternative way
plot(variog(sim1, max.dist=1))
lines.variomodel(sim1)
#
# a "smallish" simulation
sim2 <- grf(441, grid = "reg", cov.pars = c(1, .25)) 
image(sim2)
##
## 1-D simulations using the same seed and different noise/signal ratios
##
set.seed(234)
sim11 <- grf(100, ny=1, cov.pars=c(1, 0.25), nug=0)
set.seed(234)
sim12 <- grf(100, ny=1, cov.pars=c(0.75, 0.25), nug=0.25)
set.seed(234)
sim13 <- grf(100, ny=1, cov.pars=c(0.5, 0.25), nug=0.5)
##
par.ori <- par(no.readonly = TRUE)
par(mfrow=c(3,1), mar=c(3,3,.5,.5))
yl <- range(c(sim11$data, sim12$data, sim13$data))
image(sim11, type="l", ylim=yl)
image(sim12, type="l", ylim=yl)
image(sim13, type="l", ylim=yl)
par(par.ori)

## simulating within borders
data(parana)
pr1 <- grf(100, cov.pars=c(200, 40), borders=parana$borders, mean=500)
points(pr1)
pr1 <- grf(100, grid="reg", cov.pars=c(200, 40), borders=parana$borders)
points(pr1)
pr1 <- grf(100, grid="reg", nx=10, ny=5, cov.pars=c(200, 40), borders=parana$borders)
points(pr1)
}

\keyword{spatial}
\keyword{datagen}

