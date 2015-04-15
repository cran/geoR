\name{plot.variogram}

\alias{plot.variogram}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Empirical Variogram}

\description{
  Plots sample (empirical) variogram computed using the
  function \code{\link{variog}}.
}

\usage{
\method{plot}{variogram}(x, max.dist, vario.col = "all", scaled = FALSE,
               var.lines = FALSE, envelope.obj = NULL,
               pts.range.cex, bin.cloud = FALSE,  \dots)
}


\arguments{
  \item{x}{an object of the class \code{"variogram"}, typically an
    output of the function \code{\link{variog}}.  }
  \item{max.dist}{maximum distance for the x-axis. The default is
    the maximum distance for which the sample variogram was computed.}
  \item{vario.col}{only used if \code{obj} has information on more than
    one empirical variogram. The default \code{"all"} indicates that
    variograms of all variables should be plotted.  Alternativelly a
    numerical vector can be used to select  variables. }
  \item{scaled}{If \code{TRUE} the variogram values are divided by the
    sample variance. This allows comparison of variograms 
    of variables measured in different scales.  }
  \item{var.lines}{If \code{TRUE} a horizontal line is drawn at the value
    of the variance of the data (if \code{scaled = F}) or at 1 (if
    \code{scaled = T}).  }
  \item{envelope.obj}{adds a variogram envelope computed by
    the function \code{\link{variog.model.env}} or
    \code{\link{variog.mc.env}}.  }
  \item{pts.range.cex}{optional. A two elements vector with maximum and
    minimum values for the caracter expansion factor \code{cex}. If
    provided the point sizes in binned variogram are proportional to the
    number of pairs of points used to compute each bin.  } 
  \item{bin.cloud}{logical. If \code{TRUE}
    and the sample variogram
    was computed with the option \code{bin.cloud = TRUE}, box-plots
    of values at each bin are plotted instead of the empirical
    variograms.  }
  \item{\dots}{other arguments to be passed to the function
    \code{\link{plot}} or \code{\link{matplot}}  }

}

\details{
  This function plots empirical  variograms.
  Toghether with  \code{\link{lines.variogram}}
  can be used to compare sample variograms of different variables
  and 
  to compare variogram models against the
  empirical variogram.
  
  It uses the function \code{\link{matplot}} when plotting variograms
  for more them one variable.
}

\value{
  Produces a plot with the sample variogram on the current graphics
  device.
  No values are returned.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}


\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}
}

%\note{ ~~further notes~~ }
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~
%
\seealso{
  \code{\link{variog}} for variogram calculations,
  \code{\link{lines.variogram}} and \code{\link{lines.variomodel}} for
  adding lines to the current plot,
  \code{\link{variog.model.env}} and \code{\link{variog.mc.env}} for
  variogram envelops computation, \code{\link{matplot}} for multiple
  lines plot 
  and \code{\link{plot}} for generic plot function.
}

\examples{
op <- par(no.readonly = TRUE)
sim <- grf(100, cov.pars=c(1, .2)) # simulates data
vario <- variog(sim, max.dist=1)   # computes sample variogram
par(mfrow=c(2,2))
plot(vario)                     # the sample variogram
plot(vario, scaled = TRUE)      # the scaled sample variogram
plot(vario, max.dist = 1)       # limiting the maximum distance
plot(vario, pts.range = c(1,3)) # points sizes proportional to number of pairs
par(op)
}

\keyword{spatial}
\keyword{dplot}
