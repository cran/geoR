\name{lines.variomodel.grf}

\alias{lines.variomodel.grf}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Lines with True Variogram for Simulated Data}

\description{
  This functions adds to the graphics device a line with the theoretical
  (true) variogram used when generating simulations with
  the function \code{\link{grf}}.
}

\usage{
\method{lines.variomodel}{grf}(x, max.dist = max(dist(obj$coords)),
length = 100, lwd = 2, \dots)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object from the class \code{grf} typically an output of
    the function \code{\link{grf}}.  }
  \item{max.dist}{maximum distance to compute and plot the true
    variogram. Defaults to the maximum distance between two data locations.  }
  \item{length}{number of points used to compute and draw the variogram line.}
  \item{lwd}{width of the line.  }
  \item{\dots}{further arguments to be passed to the function
    \code{\link{curve}}.  }
}

%\details{
%
%}

\value{
  A line with the true variogram model is added to the current
  plot on the graphics device. No values are returned.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{lines.variomodel}},
  \code{\link{grf}}, \code{\link{plot.grf}}, \code{\link{curve}}.}

\examples{
sim <- grf(100, cov.pars=c(1, .25)) # simulates data
plot(variog(sim, max.dist=1))       # plot empirical variogram
lines(sim, max.dist=1)              # plot true variogram
}

\keyword{spatial}
\keyword{aplot}%-- one or more ...