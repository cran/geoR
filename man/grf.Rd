\name{grf}

\alias{grf}
\alias{grf.aux2}
\alias{grf.aux1}
\alias{grfclass}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Simulation of Gaussian Random Fields}

\description{
  Generates simulations of Gaussian random fields for
  given covariance parameters.
}


\usage{
grf(n, grid = "irreg", nx, ny, xlims = c(0, 1), ylims = c(0, 1),
    nsim = 1, cov.model = "matern",
    cov.pars = stop("cov. parameters (sigmasq and phi) needed"), 
    kappa = 0.5, nugget = 0, lambda = 1, aniso.pars,
    method = c("cholesky", "svd", "eigen", "circular.embedding"),
    messages.screen = TRUE)
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
  \item{method}{simulation method. Defaults to the \emph{Cholesky}
    decomposition.  See section \code{DETAILS} below.  }
  \item{messages.screen}{logical,  indicating
    whether or not status messages are printed on the screen (or output device)
    while the function is running. Defaults to \code{TRUE}.  } 
}


\details{
  For the methods \code{"cholesky"}, \code{"svd"} and \code{"eigen"} the
  simulation consists of multiplying a vector of standardized
  normal deviates by a square root of the covariance matrix.
  The square root of a matrix is not uniquely defined. The 
  three available methods differs in the way they compute the
  square root of the (positive definite) covariance matrix.

  For \code{method = "circular.embedding"} the algorithm implements
  the method described by Wood & Chan (1994) which is based on Fourier
  transforms.
  Only regular and equally spaced grids can be generated using this method.
  \cr
  The code for the \code{"circular.embedding"} method
  was provided by Martin Schlather, University of Bayreuth\cr
  (\url{http://btgyn8.geo.uni-bayreuth.de/~martin/}).
 
  \bold{WARNING:} The code for the \code{"circular.embedding"} method is
  no longer being maintained. Martin has released a
  package called \code{RandomFields} (available on \code{CRAN})
  for simulation of random fields.
  We strongly recommend the use of this package for simulations on fine
  grids with large number of locations. 

}


\value{
  A list with the components:
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
  \item{.Random.seed}{the random seed at the time the function was
    called.  }
  \item{messages}{messages produced by the function describing the
    simulation.  }
  \item{call}{the function call.  }
}


\references{
  Wood, A.T.A. and Chan, G. (1994) Simulation of stationary Gaussian
  process in \eqn{[0,1]^d}{[0,1]^d}.
  \emph{Journal of Computatinal and Graphical Statistics}, \bold{3}, 409--432.
  
  Schlather, M. (1999) \emph{Introduction to positive definite functions
    and to unconditional simulation of random fields}. Tech. Report
  ST--99--10, Dept Maths and Stats, Lancaster University.
  
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
}


\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br},  \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


\seealso{\code{\link{plot.grf}} and \code{\link{image.grf}} for
  graphical output,
  \code{\link{coords.aniso}} for anisotropy coordinates transformation
  and, \code{\link{chol}},
  \code{\link{svd}} and \code{\link{eigen}} for methods of matrix decomposition.
}

\examples{
#
sim1 <- grf(100, cov.pars = c(1, .25))
# a display of simulated locations and values
points(sim1)   
# empirical and theoretical variograms
plot(sim1)
#
# a "smallish" simulation
sim2 <- grf(441, grid = "reg", cov.pars = c(1, .25)) 
image(sim2, xlab="X Coord", ylab="Y Coord")
#
# a "bigger" one
sim3 <- grf(40401, grid = "reg", cov.pars = c(10, .2), met = "circ") 
image(sim3, xlab="X Coord", ylab="Y Coord")
##
## 1-D simulations using the same seed and different noise/signal ratios
##
sim11 <- grf(100, ny=1, cov.pars=c(1, 0.25), nug=0)
.Random.seed <- sim11$.Random.seed
sim12 <- grf(100, ny=1, cov.pars=c(0.75, 0.25), nug=0.25)
.Random.seed <- sim11$.Random.seed
sim13 <- grf(100, ny=1, cov.pars=c(0.5, 0.25), nug=0.5)
##
par.ori <- par(no.readonly = TRUE)
par(mfrow=c(3,1), mar=c(3,3,.5,.5))
yl <- range(c(sim11$data, sim12$data, sim13$data))
image(sim11, type="l", ylim=yl)
image(sim12, type="l", ylim=yl)
image(sim13, type="l", ylim=yl)
par(par.ori)
}

\keyword{spatial}
\keyword{datagen}

