\name{statistics.predictive}

\alias{statistics.predictive}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summary statistics from predictive distributions}
\description{
  Computes summaries based on simulations of predictive distribution
  returned by
  \code{\link{krige.bayes}} and \code{\link{krige.conv}}.

}
\usage{
statistics.predictive(simuls, mean.var = TRUE, quantile, threshold,
                      sim.means, sim.vars)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{simuls}{object with simulations from the predictive distribution}
  \item{mean.var}{Logical. Indicates whether or not to compute mean and variances
    of the simulations at each location.  }
  \item{quantile}{defines quantile estimator. See
    documentation for \code{\link{output.control}} .  }
  \item{threshold}{defines probability estimator. See
    documentation for \code{\link{output.control}}.  }
  \item{sim.means}{Logical. Indicates whether or not to compute the mean
    of of the conditional simulations.  }
  \item{sim.vars}{Logical. Indicates whether or not to compute the variances
    of the conditional simulations.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  A list with one ore more of the following components.
  \item{mean }{mean at each prediction location.  }
  \item{variance}{variance at each prediction location.  }
  \item{quantiles}{quantiles, at each prediction location.}  
  \item{probabilities}{probabilities, at each prediction location, of been below
    the provided threshold.  }
  \item{sim.means }{vector with means of each conditional simulation.  }
  \item{sim.vars }{vector with variances of each conditional simulation.  }

}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~


%\examples{
%}
\keyword{spatial}% at least one, from doc/KEYWORDS
