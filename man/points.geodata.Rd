\name{points.geodata}

\alias{points.geodata}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plots Spatial Locations and Data Values}

\description{
  This function produces a plot with 
  points indicating the  data locations. Arguments can control the
  points sizes, patterns and colors. These can be set to be proportional
  to data values, ranks  or quantiles. Alternatively, points can be added
  to the current  plot.
}

\usage{
points.geodata(x, coords=x$coords, data=x$data, data.col = 1, borders = NULL,
               pt.divide=c("data.proportional","rank.proportional",
                           "quintiles", "quartiles", "deciles", "equal"),
               lambda = 1, trend = "cte", weights.divide = NULL, cex.min,
               cex.max, pch.seq, col.seq, add.to.plot = FALSE,
               round.quantiles = FALSE, graph.pars = FALSE, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a list containing elements \code{coords} and
    \code{data} described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing
    coordinates of the \eqn{n} data locations in each row.
    Defaults to \code{geodata$coords}.  }
  \item{data}{a vector or matrix with data values.
    If a matrix is provided each column is regarded as one variable or realization.
    Defaults to \code{geodata$data}.  }
  \item{data.col}{the number of the data column. Only used if
    \code{data} is a matrix with columns corresponding to different
    variables or simulations.  }
  \item{borders}{If an \eqn{n \times 2}{nx2}
    matrix or data-frame with the coordinates of the borders of the
    regions is provided, the borders are added to the plot.  }
  \item{pt.divide}{defines the division of the points in categories.
    See \code{DETAILS} below for the available options.
    Defaults to \code{pt.divide = "data.proportional"}.   }
  \item{trend}{specifies the mean part of the model. The options are:
    \code{"cte"} (constant mean - default option), \code{"1st"} (a first order polynomial
    on the coordinates), \code{"2nd"} (a second order polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend).
    If provided the trend is "removed" using the function
    \code{\link[base]{lm}} and the residuals are plotted.  }
  \item{lambda}{value of the Box-Cox transformation parameter. Two particular cases
    are \eqn{\lambda = 1}{lambda = 1} which corresponds to no
    transformation and  \eqn{\lambda = 0}{lambda = 0} corresponding to
    the log-transformation.  }
  \item{weights.divide}{if a vector of weights with the same length as
    the data is provided each  data is
    divided by the corresponding element in this vector.
    Defaults to \code{NULL}.  }
  \item{cex.min}{minimum value for the graphical parameter
    \code{cex}. This value defines the size of the point corresponding the minimum
    of the data. Defaults to 0.5.  }
  \item{cex.max}{maximum value for the graphical parameter
    \code{cex}. This value defines the size of the point corresponding the maximum
    of the data. If \code{pt.divide = "equal"} it is used to set
    the value for the
    graphical parameter \code{cex}. Defaults to 1.5.  }
  \item{pch.seq}{number(s) defining the graphical parameter \code{pch}.  }
  \item{col.seq}{number(s) defining the colors in the graphical parameter
    \code{col}.  }
  \item{add.to.plot}{logical. If \code{TRUE} the points are added
    to the current plot
    otherwise a new plot is created. Defaults to \code{FALSE}.  }
  \item{round.quantiles}{logical. Defines whether or not the values
    of the quantiles should be rounded. Defaults to \code{FALSE}.  }
  \item{graph.pars}{logical. If \code{TRUE} the graphics
    parameters used to produce the plots are returned.
    Defaults to \code{FALSE}.  }
  \item{\dots}{further arguments to be passed to the function
    \code{\link{plot}}, if \code{add.to.plot = FALSE}; or to the function
    \code{\link{points}}, if \code{add.to.plot = TRUE}.  }
}

\details{
  The points can be devided in categories and have different sizes
  and/or colours according to the argument
  \code{pt.divide}. The options are:
  
  \describe{
    \item{"data.proportional"}{sizes proportional to the data values.  }
    \item{"rank.proportional"}{sizes proportional to the rank of the
      data.  }
    \item{"quintiles"}{five different sizes according to the
      quintiles of the data.  }
    \item{"quartiles"}{four different sizes according to the
      quartiles of the data.  } 
    \item{"deciles"}{ten different sizes according to the
      deciles of the data. }
    \item{"equal"}{all points with the same size.  }
    \item{a numerical vector with quantiles}{the values in the
      vector will be used by the function \code{\link{cut}} as break
      points to divide the data in classes. }
  }

  For cases where points have different sizes the arguments
  \code{cex.min} and \code{cex.max} set the minimum and the maximum
  point sizes. Additionally,
  \code{pch.seq} can set different patterns for the points and
  \code{col.seq} can be used to define colors.
  For example, different colors
  can be used for quartiles, quintiles and deciles while a sequence of
  gray tones (or a color sequence) can be used
  for point sizes proportional to the data or their ranks.
  For more details see the section \code{EXAMPLES}.
}

\value{
  A plot is created or points are added to the current graphics device.\cr
  By default no value is returned. However, if \code{graph.pars = TRUE}
  a list with graphical parameters used to produce the plot is returned.
  According to the input options, the list has
  some or all of the
  following components:

  \item{quantiles}{the values of the quantiles used to divide the data.  }
  \item{cex}{the values of the graphics expansion parameter \code{cex}.
  } 
  \item{col}{the values of the graphics color parameter \code{col}.  }
  \item{pch}{the values of the graphics pattern parameter \code{pch}.  }

}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ }

\seealso{\code{\link{plot.geodata}} for another display of the data and
  \code{\link{points}} and \code{\link{plot}} for information on the
  generic \R functions. The documentation of
  \code{\link{par}} provides details on graphical parameters.
  For  color schemes in \R see \code{\link{gray}} and
  \code{\link{rainbow}}.  }

\examples{
if(is.R()) data(s100)
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(mgp = c(2,1,0))
points(s100, xlab="Coord X", ylab="Coord Y")
points(s100, xlab="Coord X", ylab="Coord Y", pt.divide="rank.prop")
points(s100, xlab="Coord X", ylab="Coord Y", cex.max=1.7,
               col=gray(seq(1, 0.1, l=100)), pt.divide="equal")
               # the function gray() works only for R
points(s100, pt.divide="quintile", xlab="Coord X",
               ylab="Coord Y")
par(op)
}

\keyword{spatial}
\keyword{dplot}
\keyword{aplot}
%-- one or more ...



