\name{plot.variogram}

\alias{plot.variogram}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Empirical Variogram}

\description{
  Plots sample (empirical) variogram computed using the
  function \code{\link{variog}}.
}

\usage{
plot.variogram(obj, max.dist = max(obj$u), ylim = "default",
               type = "b", scaled = FALSE, var.lines = FALSE,
               envelope.obj = NULL, bin.cloud = FALSE, ...)
}


\arguments{
  \item{obj}{an object of the class \code{"variogram"}, typically an
    output of the function \code{\link{variog}}.  }
  \item{max.dist}{maximum distance for the x-axis. The default is
    the maximum distance for which the sample variogram was computed.}
  \item{ylim}{limits for the variogram values in the y-axis. The default
    is from 0 to the maximum value in \code{obj$v}.  }
  \item{type}{type of line for the empirical variogram. The default is
    \code{"b"} (dots and lines).  For further details see documentation
    for \code{\link{lines}}.  }
  \item{scaled}{If \code{TRUE} the variogram values are divided by the
    sample variance. This allows comparison between variograms from
    different variables.  }
  \item{var.lines}{If \code{TRUE} a horizontal line is drawn at the value
    of the variance of the data (if \code{scaled = F}) or at 1 (if
    \code{scaled = T}).  }
  \item{envelope.obj}{adds a variogram envelope computed by
    the function \code{\link{variog.model.env}} or
    \code{\link{variog.mc.env}}.  }
  \item{bin.cloud}{logical. If \code{TRUE}
    and the sample variogram
    was computed with the option \code{keep.cloud = TRUE}, box-plots
    of values at each bin are plotted instead of the empirical
    variograms.  }
  \item{\dots}{other arguments to be passed to the function
    \code{\link{plot}}.  }

}

\details{
  This function allows visualization of sample variogram and
  together with \code{\link{lines.variogram}} can be
  used to compare sample variograms for different variables.
  Furthermore, together with \code{\link{lines.variomodel}} can be used
  to compare theoretical and/or fitted variogram models against the
  empirical variogram.
}

\value{
  Produces a plot with the sample variogram on the current graphics
  device.
  No values are returned.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}


\author{Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
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
  variogram envelops computation,
  and \code{\link{plot}} for generic plot function.
}

\examples{
sim <- grf(100, cov.pars=c(1, .2))# generates simulated data
vario <- variog(sim, max.dist=1)  # computes sample variogram
par(mfrow=c(1,2))
plot(vario)                       # the sample variogram
plot(vario, scaled=T)             # the scaled sample variogram

}

\keyword{spatial}
\keyword{dplot}
