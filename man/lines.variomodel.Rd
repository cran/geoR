\name{lines.variomodel}

\alias{lines.variomodel}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Line with a  Variogram Model}

\description{
  This function adds a theoretical and/or fitted variogram to the
  current plot. The variogram model to be
  added is typically an output of a variogram estimation function.
  Alternatively  a
  model cab be specified by the user. 
}


\usage{
lines.variomodel(obj, max.dist = obj$max.dist, length = 100, ...)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{an object of the class \code{variomodel} which is a list containing
    information about the model parameters.  }
  \item{max.dist}{maximum distance (x-axis) to compute and draw the line
    representing the  variogram model.
    The default is the distance given by \code{obj$max.dist}.  }
  \item{length}{number of points between \code{0} and \code{max.dist} to
    compute the values of the variogram model.  }
  \item{\dots}{arguments to be passed to the function
    \code{\link{lines}}.  }
}

\details{
  Allows theoretical and/or fitted variogram(s) be added to a plot.
  Together with \code{\link{plot.variogram}} can be
  used to compare sample variograms against fitted models returned by
  \code{\link{olsfit}}, \code{\link{wlsfit}} and/or \code{\link{likfit}}.
}

\value{
  A line with a variogram model is added to a plot on the current
  graphics device. No values are returned.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.}

%\note{ ~~further notes~~ }

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{plot.variogram}}, \code{\link{lines.variogram}},
  \code{\link{olsfit}}, \code{\link{wlsfit}}, \code{\link{likfit}}, \code{\link{lines}}. 
}

\examples{
sim <- grf(100, cov.pars=c(1, .3))
# generates simulated data
vario <- variog(sim)
# computes empirical variogram
plot(vario)
# plots the empirical variogram
vario.wls <- wlsfit(vario, ini=c(1, .3), fix.nugget = T)
# estimate parameters
lines(vario.wls)
# adds fitted model to the plot  

}

\keyword{spatial}
\keyword{aplot}
