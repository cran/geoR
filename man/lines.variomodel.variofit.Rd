\name{lines.variomodel.variofit}

\alias{lines.variomodel.variofit}

\title{Adds a Line with a Fitted Variogram Model to a Variogram Plot}

\description{
  This function adds a line with the variogram model
  fitted by the function
  \code{\link{variofit}} to a current variogram plot.
}

\usage{
\method{lines.variomodel}{variofit}(x, max.dist, scaled = FALSE, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object of the class \code{variofit}
    which is a list containing information about
    the fitted model parameters, typically an output of the function
    \code{\link{variofit}}.  }
  \item{max.dist}{maximum distance (x-axis) to compute and draw the line
    representing the variogram model.
    The default is the distance given by \code{x$max.dist}.  }
  \item{scaled}{logical. If \code{TRUE} the total sill in the plot is
    set to \eqn{1}.  }
  \item{\dots}{arguments to be passed to the function
    \code{\link{curve}}.  }
}

\details{
  Adds fitted variogram model to a plot.
  In conjuction with \code{\link{plot.variogram}} can be
  used to compare empirical variograms against fitted models returned by
  \code{\link{variofit}}.
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
  \code{\link{lines.variomodel.likGRF}},
  \code{\link{plot.variogram}}, \code{\link{lines.variogram}},
  \code{\link{variofit}}, \code{\link{likfit}}, \code{\link{curve}}. 
}

\examples{
# compute and plot empirical variogram
vario <- variog(s100, max.dist = 1)
plot(vario)
# estimate parameters
vario.wls <- variofit(vario, ini = c(1, .3), fix.nugget = TRUE)
# adds fitted model to the plot  
lines(vario.wls)
}
\keyword{spatial}
\keyword{aplot}

