\name{hist.krige.bayes}

\alias{hist.krige.bayes}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plots Sample from Posterior Distributions}

\description{
  Plots histograms and/or density estimation with samples
  from the posterior distribution of the model parameters
}

\usage{
\method{hist}{krige.bayes}(x, pars, density.est = TRUE, histogram = TRUE, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object of the class \code{krige.bayes}, with an output of
    the funtions \code{\link{krige.bayes}}. }
  \item{pars}{a vector with the names of one or more of the model parameters.
    Defaults to all model parameters. Setting to -1 excludes the intercept.  }
  \item{density.est}{logical indication whether a line with the density 
    estimation should be added to the plot.  }
  \item{histogram}{logical indicating whether the histogram is included in the 
    plot.  }
  \item{\dots}{further arguments for the plotting functions and or for the 
    density estimation.  }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  Produces a plot in the currently graphics device.\cr
  Returns a \code{\link{invisible}} list with the components:
  
  \item{histogram}{with the output of the function \code{\link{hist}} for 
    each parameter}
  \item{density.estimation}{with the output of the function 
    \code{\link{density}} for each parameter}
  
}


%\references{ ~put references to the literature/web site here ~ }

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{
  \code{\link{krige.bayes}}, \code{\link{hist}}, \code{\link{density}}.
}

\examples{
## See documentation for krige.bayes()
}

\keyword{spatial}
\keyword{dplot}
