\name{lines.variomodel.krige.bayes}

\alias{lines.variomodel.krige.bayes}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Adds a Bayesian Estimate of the Variogram to a Plot}

\description{
  Adds a Bayesian estimate of the variogram model to a plot typically with an empirical variogram.
  The estimate is a chosen summary (mean, mode or mean) of the
  posterior distribution returned by the function \code{\link{krige.bayes}}.
}

\usage{
\method{lines.variomodel}{krige.bayes}(x, summary.posterior, max.dist, uvec,
                 posterior = c("variogram", "parameters"),  \dots)
}


\arguments{
  \item{x}{an object of the class \code{krige.bayes}, typically an output
    of the function \code{\link{krige.bayes}}.  }
  \item{summary.posterior}{specify which summary of the posterior
    distribution should be used as the parameter estimate.
    Options are \code{"mean"}, \code{"median"} or
  \code{"mode"}. See \code{DETAILS} below. }
  \item{max.dist}{numerical, the maximum distance for the x-axis.  }
  \item{uvec}{a numerical vector with support points to compute the
    variogram values. Only used if \code{posterior = "variogram"}.
    Defaults to \code{seq(0, max.dist, length = 51)}.
  }
  \item{posterior}{indicates whether the the variogram line is
    based on the posterior of the variogram function (default)
    or the posterior of the model parameters.  }
  \item{\dots}{arguments passed to the functions \code{\link{lines}} or \code{\link{curve}}.  }
}

\details{
  The function \code{\link{krige.bayes}} returns samples from the
  posterior distribution of the parameters \eqn{(\sigma^2, \phi,
    \tau^{2}_{rel})}{(sigma^2, phi, tausq.rel)}.
  
  This function allows for two basic options to draw a line with a
  summary of the variogram function.
  \describe{
    \item{"variogram": }{for each sample of the parameters the variogram
      function is computed at the support points defined in the
      argument \code{uvec}. Then a function provided by the user in the
      argument \code{summary.posterior} is used to compute a summary of
      the values obtained at each support point.  
    }
    \item{"parameters": }{in this case summaries of the posterior
      distribution of the model parameters as \dQuote{plugged-in} in the
      variogram function.
      One of the options \code{"mode"} (default) ,\code{"median"} or \code{"mean"}
      can be provided in the argument \code{summary.posterior}.
      The option \code{mode}, uses the mode of \eqn{(\phi,
	\tau^{2}_{rel})}{(phi, tau^2.rel)} and the mode of
      of \eqn{\sigma^2}{sigma^2} conditional on the modes of the former parameters.
      For the options \code{mean} and \code{median} these summaries are
      computed from the samples of the posterior.
    }
  }
}

\value{
  A line with the estimated variogram plot is added to the plot in the
  current graphics device. No values are returned.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR/}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{lines.variomodel}}, \code{\link{krige.bayes}} and \code{\link{lines}}.}

\examples{
#See examples in the documentation of the function krige.bayes().
}

\keyword{spatial}
\keyword{aplot}

