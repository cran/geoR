\name{plot.geodata}

\alias{plot.geodata}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Exploratory Geostatistical Plots}

\description{
  This function produces a \eqn{2 \times 2}{2 x 2} display
  with the following plots: 
  the first indicates the spatial locations assign different
  colors to data
  in different quartiles,
  the next two shows data against the \emph{X} and
  \emph{Y} coordinates and the last is an histogram of the data values or optionally,
  a 3-D plot with spatial locations and associated data values.
}

\usage{
\method{plot}{geodata}(x, coords=x$coords, data = x$data,
             borders, trend="cte", lambda = 1, col.data = 1,
             weights.divide = "units.m", lowess = FALSE, scatter3d = FALSE,
             density = TRUE, rug = TRUE, qt.col, \dots)
}

%- maybe also `usage' for other objects do cumented here.

\arguments{
  \item{x}{a list containing elements \code{coords} and
    \code{data} described next. Typically an object of the class
    \code{"geodata"} - a \pkg{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{ an \eqn{n \times 2}{n x 2} matrix containing in each row Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
  \item{borders}{If an \eqn{n \times 2}{nx2} matrix or data-frame with
    the borders of the area is provided, the borders are included in the
    first plot. By default it searches for a element named "borders"
    in  the geodata object.  }
  \item{trend}{specifies the mean part of the model. The options are:
    \code{"cte"} (constant mean - default option), \code{"1st"} (a first order polynomial
    on the coordinates), \code{"2nd"} (a second order polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend).
    If provided the trend is "removed" using the function
    \code{\link[stats]{lm}} and the residuals are plotted.  }
  \item{lambda}{value of the Box-Cox transformation parameter. Two particular cases
    are \eqn{\lambda = 1}{lambda = 1} which corresponds to no
    transformation and  \eqn{\lambda = 0}{lambda = 0} corresponding to
    the log-transformation.  }
  \item{col.data}{indicates the column number for the data
    to be plotted. Only valid if more than one data-set is available
    i.e., if the argument \code{data} is a matrix.  }
  \item{weights.divide}{if a vector of weights with the same length as
    the data is provided each  data is
    divided by the corresponding element in this vector.
    Defaults divides the data by the element \code{units.m} in the
    data object, if present, otherwise no action is taken and original
    data is used.
    The usage of \code{units.m} is common for data objects
    to be analysed using the package \pkg{geoRglm}.  }
  \item{lowess}{logical. Indicates whether the function
    \code{\link{lowess}} should be used in the plots of the data against
  the coordinates.  }
  \item{scatter3d}{logical. If \code{TRUE} the last plot is
    produced by
    \code{\link[scatterplot3d]{scatterplot3d}} showing
    a 3d plot with data locations and corresponding values.  }
  \item{density}{logical. If \code{TRUE} (default) a line with density
    estimation is added to the histogram.  }
  \item{rug}{logical. If \code{TRUE} a rug plot is added to the
    histogram. }
  \item{qt.col}{colors for the quartiles in the first plot. If missing
    defaults to blue, green, yellow and red.}
  \item{\dots}{further arguments to be passed to the function
    \code{\link{hist}} or
    \code{\link[scatterplot3d]{scatterplot3d}}.  }
}

%\details{
%  By default, this function requires the package
%  \pkg{scatterplot3d} in order to produce a 3-D plot with data
%  locations and coordinates. If this package is not available an
%  histogram of the data replaces the 3-D plot.
%}

\value{
  A plot is produced on the graphics device. No values are returned. 
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

\seealso{\code{\link{points.geodata}},
  \code{\link[scatterplot3d]{scatterplot3d}}, \code{\link{lowess}}, 
  \code{\link{density}}, \code{\link{rug}}. }

\examples{
require(geoR)
plot(s100)
plot(s100, scatter3d=TRUE)
plot(s100, qt.col=1)

plot(ca20, bor=borders)          # original data
plot(ca20, trend=~altitude+area) # residuals from an external trend
plot(ca20, trend='1st')          # residuals from a polynomial trend

plot(sic.100, bor=sic.borders)           # original data
plot(sic.100, bor=sic.borders, lambda=0) # logarithm of the data
}

\keyword{spatial}
\keyword{dplot}%-- one or more ...
