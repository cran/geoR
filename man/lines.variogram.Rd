\name{lines.variogram}

\alias{lines.variogram}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Line with a Empirical Variogram}

\description{
  A sample (empirical) variogram computed using the
  function \code{\link{variog}} is added to the current plot.
}

\usage{
lines.variogram(obj, max.dist = max(obj$u), type = "o",
                scaled = F, ...)
}

\arguments{
  \item{obj}{an object of the class \code{"variogram"}, typically an
    output from the function \code{\link{variog}}.  }
  \item{max.dist}{maximum distance for the x-axis. The default is
    the maximum distance for which the sample variogram was computed.  }
  \item{type}{type of line for the empirical variogram. The default is
    \code{"o"} (dots and lines). See documentation for
    \code{\link{lines}} for further details.  }
  \item{scaled}{logical. If \code{TRUE} the variogram values are divided by the
    sample variance. This allows comparison between variograms of
    different variables.  }
  \item{\dots}{other arguments to be passed to the function
    \code{\link{lines}}.  }
}

%\details{
%  Adds a sample variogram(s) to a plot.
%  Together with \code{\link{lines.variogram}} can be
%  used to compare sample variograms for different variables.
%}


\value{
  A line with the empirical variogram is added to the plot in the current
  graphics device. No values are returned.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{variog}}, \code{\link{lines.variogram}},
  \code{\link{lines.variomodel}}, \code{\link{variog.model.env}},
  \code{\link{variog.mc.env}}, \code{\link{plot.grf}}, \code{\link{lines}}.}

%\examples{
%  TO DO: ADD AN EXAMPLE WITH 2 VARIABLES AND SCALED VARIOGRAMS
%}

\keyword{spatial}
\keyword{aplot}


