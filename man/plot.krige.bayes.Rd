\name{plot.krige.bayes}

\alias{plot.krige.bayes}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plots Prior and/or Posterior Distributions}

\description{
  Produces
  plots the priors and posteriors distribuitions for the paramters
  \code{phi} and \code{tausq.rel} based on results returned by
  \code{\link{krige.bayes}}.
}

\usage{
\method{plot}{krige.bayes}(x, phi.dist = TRUE, tausq.rel.dist = TRUE, add = FALSE,
                 type=c("bars", "h", "l", "b", "o", "p"), thin, \dots)

}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object of the class \code{krige.bayes}, with an output of
    the funtions \code{\link{krige.bayes}}. }
  \item{phi.dist}{logical indicating whether or not plot the
    distributions for this parameter.  }
  \item{tausq.rel.dist}{logical indicating whether or not plot the
    distributions for this parameter.  }
  \item{add}{logical. If \code{TRUE} plots is added to current one.  }
  \item{type}{indicates the type of plot. Option \code{"bars"} uses the
    function \code{\link[graphics]{barplot}} and the others uses
    \code{\link[graphics]{matplot}}. }
  \item{thin}{a numerical vector defining the thining for values of the
    parameters \code{phi} and \code{tausq.rel} respectively.
    This improves visualisation when there are many values
    in the discrete distribution of the parameters.  }
  \item{\dots}{further arguments for the plotting function.  }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  For \code{plot.krige.bayes} a plot is produced or added to the current
  graphics device. No values are returned.
  
}


%\references{ ~put references to the literature/web site here ~ }

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
\code{\link{krige.bayes}}, \code{\link{barplot}}, \code{\link{matplot}}.
}

\examples{
## See documentation for krige.bayes
}

\keyword{spatial}
\keyword{dplot}
\keyword{aplot}
