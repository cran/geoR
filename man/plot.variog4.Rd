\name{plot.variog4}
\alias{plot.variog4}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Directional Variograms}

\description{
  This function plot directional variograms computed in four different directions.
}

\usage{
plot.variog4(obj, omnidirectional=FALSE, same.plot=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{obj}{an object of the class \code{variog4}, typically an output of the function \code{\link{variog4}}. }
  \item{omnidirectional}{logical. Indicates whether the omnidirectional
    variogram is included in the plot. }
  \item{same.plot}{logical. Indicates whether the directional variograms
  are plotted in the same or separated plots.  }
  \item{\dots}{further arguments to be passed to the function
    \code{\link{plot}}.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  A plot is produced on the output device. No values returned.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{variog4}} for variogram calculations and
  \code{\link{plot.variog4}} for plottig results.  }

\examples{
if(is.R()) data(s100)
plot(variog4(s100, max.dist=1))
plot(variog4(s100, max.dist=1), omni=TRUE)
plot(variog4(s100, max.dist=1), omni=TRUE, same=F)
}
\keyword{spatial}
\keyword{dplot}