\name{plot.geodata}

\alias{plot.geodata}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Exploratory Geostatistical Plots}

\description{
  This function produces a \eqn{2 \times 2}{2 x 2} graphics display
  with the following plots: 
  the first indicates spatial locations,
  the second is a  3-D plot with spatial locations and associated data values
  (or a histogram of the data), and the last two shows data against the \emph{X} and
  \emph{Y} coordinates.   
}

\usage{
plot.geodata(geodata, coords=geodata$coords, data = geodata$data,
     borders = NULL, trend="cte", lambda = 1, col.data = 1,
     weights.divide = FALSE, window.new = FALSE, ...)
}

%- maybe also `usage' for other objects do cumented here.

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{ an \eqn{n \times 2}{n x 2} matrix containing in each row Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
  \item{borders}{If an \eqn{n \times 2}{nx2} matrix or data-frame with
    the borders of the area is provided, the borders are included in the
    first plot.  }
  \item{trend}{specifies the mean part of the model. The options are:
    \code{"cte"} (constant mean - default option), \code{"1st"} (a first degree polynomial
    on the coordinates), \code{"2nd"} (a second degree polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend).
    If provided the trend is "removed" using the function
    \code{\link[base]{lm}} and the residuals are plotted.  }
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
    Defaults to \code{NULL}.  }
  \item{window.new}{logical. If \code{TRUE} a new graphic device is openned,
    otherwise the current one is used. Defaults to \code{FALSE}.  }
  \item{\dots}{further arguments to be passed to the function
    \code{\link[scatterplot3d]{scatterplot3d}}.  }
}

\details{
  By default, this function requires the package
  \code{scatterplot3d} in order to produce a 3-D plot with data
  locations and coordinates. If this package is not available an
  histogram of the data replaces the 3-D plot.
}

\value{
  A plot is produced on the graphics device. No values are returned. 
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

\seealso{\code{\link{points.geodata}}, \code{\link[scatterplot3d]{scatterplot3d}}. }

\examples{
if(is.R()) data(s100)
plot(s100)
}

\keyword{spatial}
\keyword{dplot}%-- one or more ...
