\name{summary.variomodel}

\alias{summary.variomodel}
\alias{print.summary.variomodel}
\alias{print.variomodel}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summarize Results of Variogram Estimation}

\description{
  This function prints a summary of the parameter estimation results given
  by \code{\link{olsfit}} or \code{\link{wlsfit}}.
}

\usage{
summary.variomodel(obj)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{an object of the class \code{"variomodel"}
    typically an output of
     \code{\link{olsfit}} or
    \code{\link{wlsfit}}.}
}

%\details{
%}

\value{
  Prints a summary of the estimation results on the screen or
  other output device. 
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

\seealso{The functions \code{\link{olsfit}} and \code{\link{wlsfit}} for
  variogram based estimation. For likelihood based parameter estimation see  \code{\link{likfit}}. }

\examples{
if(is.R()) data(s100)
s100.vario <- variog(s100, max.dist=1)
wls <- wlsfit(s100.vario, ini=c(.5, .5), fix.nugget = TRUE)
wls
summary(wls)
}

\keyword{spatial}












