\name{pred_grid}

\alias{pred_grid}

\title{Generates a 2D Prediction Grid }
\description{
  This function facilitates the generation of a 2D prediction grid for
  geostatistical kriging.
}
\usage{
pred_grid(coords, y.coords = NULL, \dots, y.by = NULL,
          y.length.out = NULL, y.along.with = NULL)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{coords}{a list, matrix or data-frame with xy-coordinates
    of prediction points or a vector with x-coordinates.  }
  \item{y.coords}{a vector with y-coordinates. Needed if
    argument \code{coords} provides only x-coordinates.  }
  \item{\dots}{arguments \code{by} or \code{length.out}
    to be passed to the function \code{\link{rep}}.
    These arguments are used for the x-coordinates and are
    default optioons for y-coordinates.  }
  \item{y.by}{Optional. \code{by} argument for \code{\link{rep}} to be
    used with the y-coordinates. }
  \item{y.length.out}{Optional. \code{length.out} argument for
    \code{\link{rep}} to be used with the y-coordinates. }
  \item{y.along.with}{Optional. \code{along.with} argument for \code{\link{rep}} to be
      used with the y-coordinates. }
}

%      \details{
%  ~~ If necessary, more details than the description above ~~
}
\value{
 An two column data-frame which is on output of \code{\link{expand.grid}}. 
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 

\seealso{See \code{\link{seq}} and \code{\link{expand.grid}} which are
  used internally and  \code{\link{locations.inside}} and
  \code{\link{polygrid}} to select points inside a border.
}

\examples{
data(ca20)
loc0 <- pred_grid(ca20$borders, by=20)
points(ca20, borders=borders)
points(loc0, pch="+")
points(locations.inside(loc0, ca20$border), pch="+", col=2)
}

\keyword{spatial}
% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
