\name{summary.variofit}

\alias{summary.variofit}
\alias{print.summary.variofit}
\alias{print.variofit}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summarize Results of Variogram Estimation}

\description{
  This function prints a summary of the parameter estimation results given
  by \code{\link{variofit}}.
}

\usage{
\method{summary}{variofit}(object, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{object}{an object of the class \code{"variomodel"}
    typically an output of
    \code{\link{variofit}}.}
  \item{\dots}{other arguments to be passed to the function
    \code{\link{print}} or \code{\link{summary}}.  }

}

%\details{
%}

\value{
  Prints a summary of the estimation results on the screen or
  other output device. 
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

\seealso{The functions \code{\link{variofit}} for
  variogram based estimation. For likelihood based parameter estimation see  \code{\link{likfit}}. }

\examples{
s100.vario <- variog(s100, max.dist=1)
wls <- variofit(s100.vario, ini=c(.5, .5), fix.nugget = TRUE)
wls
summary(wls)
}

\keyword{spatial}












