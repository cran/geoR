\name{lines.krige.bayes}

\alias{lines.krige.bayes}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Add a (Bayesian) Estimated Variogram to a Plot}

\description{
  Adds a estimated variogram to a plot with an empirical variogram.
  The estimate is a chosen summary (mean, mode or mean) of the
  posterior distribution returned by the function \code{\link{krige.bayes}}.
}

\usage{
lines.krige.bayes(obj, max.dist, length = 100,
             summary.posterior = c("mode", "median", "mean"), ...)
}


\arguments{
  \item{obj}{an object of the class \code{krige.bayes}, typically output of the function
    \code{\link{krige.bayes}}.  }
  \item{max.dist}{maximum distance for the x-axis.  }
  \item{length}{number of points to compute the variogram values.  }
  \item{summary.posterior}{specify which summary of the posterior
    distribution should be used as parameter estimation. Options are
    \code{"mean"}, \code{"median"} or
  \code{"mode"}.  }
  \item{\dots}{arguments to be passed to the function
    \code{\link{lines}}.  }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  A line with the estimated variogram plot is added to the plot in the
  current graphics device. No values are returned.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{krige.bayes}} and \code{\link{lines}}.}

\examples{
#See examples in the documentation of the function krige.bayes().
}

\keyword{spatial}
\keyword{aplot}

