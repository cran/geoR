\name{lines.variomodel.likGRF}

\alias{lines.variomodel.likGRF}

\title{Adds a Variogram Line to a Variogram Plot}

\description{
  This function adds a fitted variogram based on the estimates of the
  model parameters returned by the function \code{\link{likfit}}
  to a current variogram plot.
}


\usage{
\method{lines.variomodel}{likGRF}(x, max.dist, scaled = FALSE, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object of the class \code{likGRF}
    which is a list containing
    information about the fitted model parameters, typically an output
    of the function \code{\link{likfit}}.  }
  \item{max.dist}{maximum distance (x-axis) to compute and draw the line
    representing the  variogram model.
    The default is the distance given by \code{obj$max.dist}.  }
  \item{scaled}{logical. If \code{TRUE} the total sill in the plot is
    equals to \eqn{1}.  }
  \item{\dots}{arguments to be passed to the function
    \code{\link{curve}}.  }
}

\details{
  Adds variogram model(s)  to a plot.
  In conjuction with \code{\link{plot.variogram}} can be
  used to compare sample variograms against fitted models returned by
  \code{\link{likfit}}.
}

\value{
  A line with a variogram model is added to a plot on the current
  graphics device. No values are returned.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.}

%\note{ ~~further notes~~ }

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{lines.variomodel}},
  \code{\link{lines.variomodel.variofit}},
  \code{\link{plot.variogram}}, \code{\link{lines.variogram}},
  \code{\link{variofit}}, \code{\link{likfit}}, \code{\link{curve}}. 
}

\examples{
# compute and plot empirical variogram
vario <- variog(s100, max.dist = 1)
plot(vario)
# estimate parameters
vario.ml <- likfit(s100, ini = c(1, .3), fix.nugget = TRUE)
# adds fitted model to the plot  
lines(vario.ml)
}
\keyword{spatial}
\keyword{aplot}

