\name{image.kriging}

\alias{image.kriging}
\alias{persp.kriging}
\alias{prepare.graph.kriging}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Image or Perspective Plot with Kriging Results}

\description{
  Plots image or perspective plots with results of the kriging
  calculations.
}

\usage{
\method{image}{kriging}(x, locations, borders, values = obj$predict,
              coords.data, x.leg, y.leg, cex.leg = 0.75,
              vertical = FALSE, ...)

persp.kriging(x, locations, borders, values = obj$predict, ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object of the  class \code{kriging}, typically
    with the output of the functions \code{\link{krige.conv}} or \code{\link{ksline}}.  }
  \item{locations}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates of the prediction locations, which should define a
    regular grid
    in order to be plotted by \code{\link{image}} or
    \code{\link{persp}}. By default does not need to be provided and
    evaluates the attribute "prediction.locations" from the input object. }
  \item{borders}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates defining the borders of a region inside the grid defined
    by \code{locations}. Elements in the
    argument \code{values} are assigned to locations internal to the
    borders and \code{NA}'s to the external ones.
  }
  \item{values}{a vector with values to be plotted. Defaults to \code{obj$predict}.  }
  \item{coords.data}{optional. If an \eqn{n \times 2}{n x 2}
    matrix with the data coordinates is
    provided, points indicating the data locations are included in the
    plot.  }
  \item{x.leg}{limits for the legend in the \eqn{x} direction.  }
  \item{y.leg}{limits for the legend in the \eqn{y} direction.  }
  \item{cex.leg}{controls the size of the numerical values in the legend.  }
  \item{vertical}{logical, defining whether the legend is drawn in
    horizontal (default) or vertical position.  }
  \item{\dots}{further arguments to be passed to the functions
    \code{\link{image}} or \code{\link{persp}}.  } 
  
}

%\details{
%}

\value{
  An image or perspective plot is produced o the current graphics
  device.
  No values are returned.
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

 
\seealso{\code{\link{krige.conv}} and \code{\link{ksline}} for kriging
  calculations. Documentation for
  \code{\link{image}} and  \code{\link{persp}} contain basic information
  on the plotting functions.}

\examples{
if(is.R()) data(s100) 
loci <- expand.grid(seq(0,1,l=31), seq(0,1,l=31))
kc <- krige.conv(s100, loc=loci,
                 krige=krige.control(cov.pars=c(1, .25)))
par(mfrow=c(1,2))
image.kriging(kc, loc=loci)
image.kriging(kc, loc=loci, val=kc$krige.var)

}

\keyword{spatial}
\keyword{dplot}
