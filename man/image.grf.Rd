\name{image.grf}

\alias{image.grf}
\alias{persp.grf}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Image or Perspective Plot of Simulated Gaussian Random Field}

\description{
  Plots an image or perspective plot with a realization of a Gaussian
  random field, simulated using the function \code{\link{grf}}.  
}

\usage{
\method{image}{grf}(x, sim.number = 1, xlim, ylim, x.leg, y.leg, \dots)
\dots)

\method{persp}{grf}(x, sim.number = 1, \dots)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object of the class \code{grf}, typically an output of the function
    \code{\link{grf}}.  }
  \item{sim.number}{simulation number.
    Indicates the number of the simulation top be plotted.
    Only valid if the object contains more than one simulation.
    Defaults to 1.  }
  \item{xlim, ylim} {the ranges to be encompassed by the x and y axes.  }
  \item{x.leg, y.leg}{limits for the legend in the horizontal and vertical directions.  }
%  \item{cex.leg}{controls the size of the numerical values in the legend.  }
%  \item{vertical}{logical, defining whether the legend is drawn in
%    horizontal (default) or vertical position.  }
  \item{\dots}{further arguments to be passed to the functions
    \code{\link{image}} or \code{\link{persp}}.  }
}

%\details{
%}

\value{
  An image or perspective plot is produced on the current graphics
  device.
  No values are returned.
}


\references{
  Further information about the package \pkg{geoR} can be found at:\cr
  \url{http://www.est.ufpr.br/geoR}.
}


\author{Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{grf}} for simulation of Gaussian random fields,
  \code{\link{image}} and \code{\link{persp}} for the generic plotting functions.}

\examples{
# generating 4 simulations of a Gaussian random field
sim <- grf(441, grid="reg", cov.pars=c(1, .25), nsim=4)
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp = c(2,1,0))
for (i in 1:4)
  image(sim, sim.n=i)
par(op)
}

\keyword{spatial}
\keyword{dplot}
