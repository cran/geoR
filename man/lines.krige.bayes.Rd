\name{lines.krige.bayes}

\alias{lines.krige.bayes}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Add a (Bayesian) Estimated Variogram to a Plot}

\description{
  Adds a variogram model to a plot with an empirical variogram.
  The estimate is a chosen summary (mean, mode or mean) of the
  posterior distribution returned by the function \code{\link{krige.bayes}}.
}

\usage{
lines.krige.bayes(x, max.dist,
                  summary.posterior = c("mode", "median", "mean"), ...)
}


\arguments{
  \item{x}{an object of the class \code{krige.bayes}, typically output
    of the function \code{\link{krige.bayes}}.  }
  \item{max.dist}{maximum distance for the x-axis.  }
  \item{summary.posterior}{specify which summary of the posterior
    distribution should be used as parameter estimation. Options are
    \code{"mean"}, \code{"median"} or
  \code{"mode"}. See \code{DETAILS} below. }
  \item{\dots}{arguments to be passed to the function
    \code{\link{lines}}.  }
}

\details{
  The functions \code{\link{krige.bayes}} produces samples of the
  posterior distribution of the parameters \eqn{(\sigma^2, \phi and
    \tau^{2}_{rel})}{(sigma^2, phi and tausq.rel)}.
  For the option \code{mode}, the mode of \eqn{(\phi and
    \tau^{2}_{rel})}{(phi and tausq.rel)} is used with the conditional
  (on thses two values) mode of \eqn{\sigma^2}{sigma^2}.
  For the options \code{mean} and \code{median} these summaries are
  computed from the samples of the posterior.
}

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

