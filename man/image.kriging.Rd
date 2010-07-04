\name{image.kriging}

\alias{image.kriging}
\alias{contour.kriging}
\alias{persp.kriging}
\alias{.prepare.graph.kriging}
\alias{plot.1d}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Image or Perspective Plot with Kriging Results}

\description{
  Plots image or perspective plots with results of the kriging
  calculations.
}

\usage{
\method{image}{kriging}(x, locations, borders, values = x$predict,
              coords.data, x.leg, y.leg, \dots)
\method{contour}{kriging}(x, locations, borders, values = x$predict,
              coords.data, filled=FALSE, \dots)
\method{persp}{kriging}(x, locations, borders, values = x$predict, \dots)
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
  \item{x.leg, y.leg}{limits for the legend in the horizontal and vertical directions.  }
  \item{filled}{logical. If \code{FALSE} the function
    \code{\link{contour}} is used otherwise
    \code{\link{filled.contour}}. Defaults to \code{FALSE}.  }
  \item{\dots}{further arguments to be passed to the functions
    \code{\link{image}}, \code{\link{contour}}, \code{\link{filled.contour}}, \code{\link{persp}} or
    \code{\link{legend.krige}}.
    For instance, the
    argument \code{zlim} can be used to set the the minimum and maximum
    `z' values for which colors should be plotted.
    See documentation for those function for
    possible arguments.
  }   
}

\details{
  \code{plot1d} and \code{prepare.graph.kriging} are auxiliary functions
  called by the others.
}

\value{
  An image or perspective plot is produced o the current graphics
  device.
  No values are returned.
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

 
\seealso{\code{\link{krige.conv}} and \code{\link{ksline}} for kriging
  calculations. Documentation for
  \code{\link{image}}, \code{\link{contour}}, \code{\link{filled.contour}} and  \code{\link{persp}} contain basic information
  on the plotting functions.}

\examples{
loci <- expand.grid(seq(0,1,l=51), seq(0,1,l=51))
kc <- krige.conv(s100, loc=loci,
                 krige=krige.control(cov.pars=c(1, .25)))
image(kc)
contour(kc)
image(kc)
contour(kc, add=TRUE, nlev=21)
persp(kc, theta=20, phi=20)
contour(kc, filled=TRUE)
contour(kc, filled=TRUE, color=terrain.colors)
contour(kc, filled=TRUE, col=gray(seq(1,0,l=21)))
# adding data locations
image(kc, coords.data=s100$coords)
contour(kc,filled=TRUE,coords.data=s100$coords,color=terrain.colors)
#
# now dealing with borders
#
bor <- matrix(c(.4,.1,.3,.9,.9,.7,.9,.7,.3,.2,.5,.8),
              ncol=2)
# plotting just inside borders
image(kc, borders=bor)
contour(kc, borders=bor)
image(kc, borders=bor)
contour(kc, borders=bor, add=TRUE)
contour(kc, borders=bor, filled=TRUE, color=terrain.colors)
# kriging just inside borders
kc1 <- krige.conv(s100, loc=loci,
                 krige=krige.control(cov.pars=c(1, .25)),
                 borders=bor)
image(kc1)
contour(kc1)
# avoiding the borders 
image(kc1, borders=NULL)
contour(kc1, borders=NULL)

op <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(3,3,0,0), mgp=c(1.5, .8,0))
image(kc)
image(kc, val=sqrt(kc$krige.var))

# different ways to add the legends and pass arguments:
image(kc, ylim=c(-0.2, 1), x.leg=c(0,1), y.leg=c(-0.2, -0.1))
image(kc, val=kc$krige.var, ylim=c(-0.2, 1))
legend.krige(y.leg=c(-0.2,-0.1), x.leg=c(0,1), val=sqrt(kc$krige.var))

image(kc, ylim=c(-0.2, 1), x.leg=c(0,1), y.leg=c(-0.2, -0.1), cex=1.5)
image(kc, ylim=c(-0.2, 1), x.leg=c(0,1), y.leg=c(-0.2, -0.1), offset.leg=0.5)

image(kc, xlim=c(0, 1.2))
legend.krige(x.leg=c(1.05,1.1), y.leg=c(0,1), kc$pred, vert=TRUE)
image(kc, xlim=c(0, 1.2))
legend.krige(x.leg=c(1.05,1.1), y.leg=c(0,1),kc$pred, vert=TRUE, off=1.5, cex=1.5)

par(op)
}

\keyword{spatial}
\keyword{dplot}
