\name{plot.proflik}

\alias{plot.proflik}
\alias{.proflik.plot.aux1}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plots Profile Likelihoods}

\description{
  This function produces plots of the profile likelihoods computed by
  the function \code{\link{proflik}}.
}

\usage{
\method{plot}{proflik}(x, pages = c("user", "one", "two"), uni.only, bi.only,
             type.bi = c("contour", "persp"), conf.int = c(0.90, 0.95),
             yaxis.lims = c("conf.int", "as.computed"),
             by.col = TRUE, log.scale = FALSE, use.splines = TRUE,
             par.mar.persp = c(0, 0, 0, 0), ask = FALSE, \dots)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{an object of the class \code{proflik}, typically an
    output of the function \code{\link{proflik}}.  }
  \item{pages}{specify how the plots will be arranged in the
    graphics device. The default option, \code{"user"}, uses the current
    graphical parameters. The option \code{"one"} places all the
    profiles in the same page and the option \code{"two"} places the
    univariate profiles in one page and the bivariate profiles in
    a second page.  }
  \item{uni.only}{only 1-D profiles are plotted even if the object
    contains results about the 2-D profiles.  }
  \item{bi.only}{only 2-D profile are plotted even if the object
    contains results about the 1-D profiles.  }
  \item{type.bi}{Type of plot for the 2-D profiles. Options are
    \code{"contour"} for  contour plot
    (the default) and \code{"persp"} for perspective plot.  }
  \item{conf.int}{a vector with numbers in the interval \eqn{[0,1]}
    specifying levels of the (approximated) confidence
    intervals. Defaults corresponds to the levels  90\% and 95\%.  }
  \item{yaxis.lims}{defines the lower limits for the y-axis in the 1-D
    plots. If \code{"conf.int"} the limit is determined by the level of
    the confidence interval (the default) otherwise will be determined
    by the smallest computed value.  }
  \item{by.col}{logical, If \code{TRUE} the plots are arranged by columns in a multiple graphics
    device.  }
  \item{log.scale}{plots the x-axis in the logarithmic scale. Defaults to
    \code{FALSE}.  }
  \item{use.splines}{logical. If \code{TRUE} (the default) the function
    \code{\link{spline}} is used to interpolate between the points
    computed by \code{proflik}.  }
  \item{par.mar.persp}{graphical parameters to be used with
    \code{\link{persp}} plots. For more details see \code{\link{par}}.  }
  \item{ask}{logical. Defines whether or not the user is prompted before
    each plot is produced.  }
  \item{\dots}{additional arguments to be passed to the functions
    \code{\link{plot}}, \code{\link{contour}} and/or \code{\link{persp}}.  }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  Produces plots with the profile likelihoods on the current graphics device.
  No values are returned.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{proflik}} for computation of the profile
  likelihoods. For the generic plotting functions see
  \code{\link{plot}},  \code{\link{contour}}, \code{\link{persp}}.
  See \code{\link{spline}} for interpolation.  }

\examples{
# see examples in the documentation for the function proflik()
}
\keyword{spatial}
\keyword{dplot}

