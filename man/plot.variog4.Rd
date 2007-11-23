\name{plot.variog4}
\alias{plot.variog4}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Directional Variograms}

\description{
  This function plot directional variograms computed by the function
  \code{\link{variog4}}. The omnidirectional variogram can be also included
  in the plot.
}

\usage{
\method{plot}{variog4}(x, omnidirectional=FALSE, same.plot=TRUE, legend = TRUE, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object of the class \code{variog4}, typically an output of the function \code{\link{variog4}}. }
  \item{omnidirectional}{logical. Indicates whether the omnidirectional
    variogram is included in the plot. }
  \item{same.plot}{logical. Indicates whether the directional variograms
    are plotted in the same or separated plots.  }
  \item{legend}{logical indicating whether legends are automatically
    included in the plots.  }
  \item{\dots}{further arguments to be passed to the function
    \code{\link{plot}}. Typical arguments are \code{col}, \code{lty},
    \code{lwd}. For \code{same.plot = TRUE} the arguments are passed to
    the function 
    \code{\link{matplot}} which is used to produce the plot.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  A plot is produced on the output device. No values returned.
}

\references{
  Further information about the \pkg{geoR} package can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{variog4}} for variogram calculations and
  \code{\link{matplot}} for multiple lines plotting.  }

\examples{
s100.v4 <- variog4(s100, max.dist=1)
# Plotting variograms for the four directions
plot(s100.v4)
# changing plot options
plot(s100.v4, lwd=2)
plot(s100.v4, lty=1, col=c("darkorange", "darkblue", "darkgreen","darkviolet"))
plot(s100.v4, lty=1, lwd=2)
# including the omnidirectional variogram
plot(s100.v4, omni=TRUE)
# variograms on different plots
plot(s100.v4, omni=TRUE, same=FALSE)
}
\keyword{spatial}
\keyword{dplot}

