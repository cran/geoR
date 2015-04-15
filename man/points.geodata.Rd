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
\method{points}{geodata}(x, coords=x$coords, data=x$data, data.col = 1, borders,
               pt.divide=c("data.proportional","rank.proportional",
                           "quintiles", "quartiles", "deciles", "equal"),
               lambda = 1, trend = "cte", abs.residuals = FALSE,
               weights.divide = "units.m", cex.min, cex.max, cex.var,
               pch.seq, col.seq, add.to.plot = FALSE,
               x.leg, y.leg = NULL, dig.leg = 2, 
               round.quantiles = FALSE, permute = FALSE, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a list containing elements \code{coords} and
    \code{data} described next. Typically an object of the class
    \code{"geodata"} - a \pkg{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  
}
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
    regions is provided, the borders are added to the plot.
    By default it searches for a element named "borders" in  the geodata
    object.  }
  \item{pt.divide}{defines the division of the points in categories.
    See \code{DETAILS} below for the available options.
    Defaults to \code{pt.divide = "data.proportional"}.   }
  \item{trend}{specifies the mean part of the model. The options are:
    \code{"cte"} (constant mean - default option), \code{"1st"} (a first order polynomial
    on the coordinates), \code{"2nd"} (a second order polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend).
    If provided the trend is "removed" using the function
    \code{\link[stats]{lm}} and the residuals are plotted.  }
  \item{abs.residuals}{logical. If \code{TRUE} and the value passed to
    the
    argument \code{trend} is different from \code{"cte"} the point sizes
    are proportional to absolute values of the residuals.  }
  \item{lambda}{value of the Box-Cox transformation parameter. Two particular cases
    are \eqn{\lambda = 1}{lambda = 1} which corresponds to no
    transformation and  \eqn{\lambda = 0}{lambda = 0} corresponding to
    the log-transformation.  }
  \item{weights.divide}{if a vector of weights with the same length as
    the data is provided each  data is
    divided by the corresponding element in this vector.
    Defaults divides the data by the element \code{units.m} in the
    data object, if present, otherwise no action is taken and original
    data is used.
    The usage of \code{units.m} is common for data objects
    to be analysed using the package \pkg{geoRglm}. }
  \item{cex.min}{minimum value for the graphical parameter
    \code{cex}. This value defines the size of the point corresponding the minimum
    of the data. Defaults to 0.5.  }
  \item{cex.max}{maximum value for the graphical parameter
    \code{cex}. This value defines the size of the point corresponding the maximum
    of the data. If \code{pt.divide = "equal"} it is used to set
    the value for the
    graphical parameter \code{cex}. Defaults to 1.5.  }
  \item{cex.var}{a numeric vector with the values of a variable defining
    the size of the points. Particularly useful for displaying 2
    variables at once.  }
  \item{pch.seq}{number(s) defining the graphical parameter \code{pch}.  }
  \item{col.seq}{number(s) defining the colors in the graphical parameter
    \code{col}.  }
  \item{add.to.plot}{logical. If \code{TRUE} the points are added
    to the current plot or image otherwise a display is open. Defaults to \code{FALSE}.  }
  \item{x.leg, y.leg}{\code{x} and \code{y} location of the legend as
    documented in \code{\link{legend}}.  }
  \item{dig.leg}{the desired number of digits after the decimal
    point. Printing values in the legend uses \code{\link{formatC}} with
    argument \code{format = "f"}. }
  \item{round.quantiles}{logical. Defines whether or not the values
    of the quantiles should be rounded. Defaults to \code{FALSE}.  }
  \item{permute}{logical indication whether the data values should be
    randomly re-alocatted to the coordinates. See \code{DETAILS} below.}
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
    \item{a scalar}{defines a number of quantiles, the number provided
      defines the  number of different points sizes and colors.  }
    \item{a numerical vector with quantiles and length > 1}{the values in the
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

  The argument \code{cex.var} allows for displaying 2 variables
  at once. In this case one variable defines the backgroung colour
  of the points and the other defines the points size.
  
  The argument \code{permute} if set to \code{TRUE}
  randomly realocates the data in the coordinates.
  This may be used to
  contrast the spatial pattern of original data against another
  situation where there is no spatial dependence (when setting
  \code{permute = TRUE}). If a \code{trend} is provided the residuals
  (and not the original data) are permuted.
}

\value{
  A plot is created or points are added to the current graphics device.\cr
  A list with graphical parameters used to produce the plot is returned invisibily.
  According to the input options, the list has some or all of the
  following components:

  \item{quantiles}{the values of the quantiles used to divide the data.  }
  \item{cex}{the values of the graphics expansion parameter \code{cex}.
  } 
  \item{col}{the values of the graphics color parameter \code{col}.  }
  \item{pch}{the values of the graphics pattern parameter \code{pch}.  }

}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
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
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp = c(2,1,0))
points(s100, xlab="Coord X", ylab="Coord Y")
points(s100, xlab="Coord X", ylab="Coord Y", pt.divide="rank.prop")
points(s100, xlab="Coord X", ylab="Coord Y", cex.max=1.7,
               col=gray(seq(1, 0.1, l=100)), pt.divide="equal")
points(s100, pt.divide="quintile", xlab="Coord X", ylab="Coord Y")
par(op)

points(ca20, pt.div='quartile', x.leg=4900, y.leg=5850, bor=borders)

par(mfrow=c(1,2), mar=c(3,3,1,1), mgp = c(2,1,0))
points(s100, main="Original data")
points(s100, permute=TRUE, main="Permuting locations")

## Now an example using 2 variable, 1 defining the
## gray scale and the other the points size
points.geodata(coords=camg[,1:2], data=camg[,3], col="gray",
               cex.var=camg[,5])
points.geodata(coords=camg[,1:2], data=camg[,3], col="gray",
               cex.var=camg[,5], pt.div="quint")
}
\keyword{spatial}
\keyword{dplot}
\keyword{aplot}



