\name{lines.variogram}

\alias{lines.variogram}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Line with a Empirical Variogram}

\description{
  A sample (empirical) variogram computed using the
  function \code{\link{variog}} is added to the current plot.
}

\usage{
\method{lines}{variogram}(x, max.dist, type = "o",  scaled = FALSE,
         pts.range.cex, \dots)
}

\arguments{
  \item{x}{an object of the class \code{"variogram"}, typically an
    output from the function \code{\link{variog}}.  }
  \item{max.dist}{maximum distance for the x-axis. By default takes
    the maximum distance for which the sample variogram was computed.  }
  \item{type}{type of line for the empirical variogram. The default is
    \code{"o"} (dots and lines). See documentation for
    \code{\link{lines}} for further details.  }
  \item{scaled}{logical. If \code{TRUE} the variogram values are divided by the
    sample variance. This allows comparison between variograms of
    different variables.  }
  \item{pts.range.cex}{optional. A two elements vector with maximum and
    minimum values for the caracter expansion factor \code{cex}. If
    provided the point sizes in binned variogram are proportional to the
    number of pairs of points used to compute each bin.  } 
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
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
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


