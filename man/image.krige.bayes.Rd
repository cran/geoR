\name{image.krige.bayes}

\alias{image.krige.bayes}
\alias{contour.krige.bayes}
\alias{persp.krige.bayes}
\alias{.prepare.graph.krige.bayes}

\title{Plots Results of the Predictive Distribution}

\description{
  This function produces an image or perspective plot of a selected
  element
  of the predictive distribution
  returned by the function \code{\link{krige.bayes}}.
}

\usage{
\method{image}{krige.bayes}(x, locations, borders,
                  values.to.plot=c("mean", "variance",
                            "mean.simulations", "variance.simulations",
                            "quantiles", "probabilities", "simulation"),
                  number.col, coords.data, x.leg, y.leg, messages, \dots)
\method{contour}{krige.bayes}(x, locations, borders, 
                  values.to.plot = c("mean", "variance",
                       "mean.simulations", "variance.simulations",
                       "quantiles", "probabilities", "simulation"),
                  filled=FALSE, number.col, coords.data,
                  x.leg, y.leg, messages, \dots)
\method{persp}{krige.bayes}(x, locations, borders,
                  values.to.plot=c("mean", "variance",
                       "mean.simulations", "variance.simulations",
                       "quantiles", "probabilities", "simulation"),
                  number.col, messages, \dots)
}

\arguments{
  \item{x}{an object of the class \code{krige.bayes}, typically an
    output of the function \code{\link{krige.bayes}}.  }
  \item{locations}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates of the prediction locations, which should define a
    regular grid
    in order to be plotted by \code{\link{image}} or
    \code{\link{persp}}. By default does not need to be provided and
    evaluates the attribute "prediction.locations" from the input object.  }
  \item{borders}{an \eqn{n \times 2}{n x 2} matrix with the
    coordinates defining the borders of a region inside the grid defined
    by \code{locations}. Elements in the
    argument \code{values} are assigned to locations internal to the
    borders  and \code{NA}'s to the external ones.  }
  \item{values.to.plot}{select the element of the predictive
    distribution to be plotted. See DETAILS below.}
  \item{filled}{logical. If \code{FALSE} the function
    \code{\link{contour}} is used otherwise
    \code{\link{filled.contour}}. Defaults to \code{FALSE}.  }
  \item{number.col}{Specifies the number of the column to be plotted.
    Only used if previous argument is set to one of \code{"quantiles"}, 
    \code{"probabilities"}  or \code{"simulation"}.  }
  \item{coords.data}{optional. If an \eqn{n \times 2}{n x 2}
    matrix with the data coordinates is
    provided, points indicating the data locations are included in the
    plot.  }
  \item{x.leg, y.leg}{limits for the legend in the horizontal and vertical directions.  }
%  \item{cex.leg}{controls the size of the numerical values in the legend.  }
%  \item{vertical}{logical, defining whether the legend is drawn in
%    horizontal (default) or vertical position.  }
  \item{messages}{logical, if TRUE status messages are printed while
    running the function.  }
  \item{\dots}{extra arguments to be passed to the plotting function
    \code{\link{image}} or \code{\link{persp}}.  }
}

\details{
  The function \code{\link{krige.bayes}} returns
  summaries and other results about the predictive distributions.
  The argument \code{values.to.plot} specifies which result will be
  plotted. It can be passed to the function in two different forms:
  \itemize{
    \item a vector with the object containing the values to be plotted, or
    \item one of the following options: \code{"moments.mean"},
    \code{"moments.variance"}, 
    \code{"mean.simulations"},
    \code{"variance.simulations"},
    \code{"quantiles"},
    \code{"probability"} or
    \code{"simulation"}.  }
  For the last three options, if the results are stored in matrices,
  a column number must be provided using the argument \code{number.col}.

  The documentation for the function \code{\link{krige.bayes}} provides
  further details about these options.
}

\value{
  An \code{\link{image}} or \code{\link{persp}} plot is produced on the
  current graphics device. No values are returned. 
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br},
  \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{krige.bayes}} for Bayesian Kriging computations and, \code{\link{image}} and  \code{\link{persp}} for the generic plotting functions.}

\examples{
#See examples in the documentation for the function krige.bayes().
}

\keyword{spatial}

