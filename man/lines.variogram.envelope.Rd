\name{lines.variogram.envelope}

\alias{lines.variogram.envelope}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Adds Envelopes Lines to a Variogram Plot}

\description{
  Variogram envelopes computed by \code{\link{variog.model.env}} or
  \code{\link{variog.mc.env}} are added to the current variogram plot.
}

\usage{
\method{lines}{variogram.envelope}(x, lty = 3, \dots)
}

\arguments{
  \item{x}{an object of the class \code{"variogram.envelope"},
    typically an output of the functions \code{\link{variog.model.env}} or
    \code{\link{variog.mc.env}}.  } 
  \item{lty}{line type. Defaults to 3.  }
  \item{\dots}{arguments to be passed to the function \code{\link{lines}}.  }
}

%\details{
%}

\value{
  Lines defining the variogram envelope
  are added to the plotin the current graphics device.
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

\seealso{\code{\link{variog}} for variogram computation,
  \code{\link{variog.model.env}} and \code{\link{variog.mc.env}} for
  computation of variogram envelopes, and \code{\link{lines}} for the
  generic function.  }

\examples{
s100.vario <- variog(s100, max.dist = 1)
s100.ml <- likfit(s100, ini=c(.5, .5))
s100.mod.env <- variog.model.env(s100, obj.variog = s100.vario,
   model = s100.ml) 
s100.mc.env <- variog.mc.env(s100, obj.variog = s100.vario)
plot(s100.vario)
lines(s100.mod.env)
lines(s100.mc.env, lwd=2)

}

\keyword{spatial}
\keyword{aplot}


