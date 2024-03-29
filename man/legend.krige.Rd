\name{legend.krige}

\alias{legend.krige}

\title{Add a legend to a image with kriging results}

\description{
  This function allows adds a legend to an image plot generated by
  \code{\link{image.kriging}} or \code{\link{image.krige.bayes}}.
  It can be called internally by these functions or directly by the user.
}

\usage{
legend.krige(x.leg, y.leg, values, scale.vals,
             vertical = FALSE, offset.leg = 1, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x.leg}{limits for the legend in the \eqn{x} direction.  }
  \item{y.leg}{limits for the legend in the \eqn{y} direction.  }
  \item{values}{values plotted in the image.  }
  \item{scale.vals}{optional. Values to appear in the legend.
    If not provided the function \code{\link{pretty}} is used to define
    the values.  }
  \item{vertical}{If \code{TRUE} the legend is drawn in the vertical
    direction. Defaults to \code{FALSE}.  }
  \item{offset.leg}{numeric value controlling the distance between the
    legend text and the legend box. }
  \item{\dots}{further arguments to be passed to the function
    \code{\link{text}}.  }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  A legend is added to the current plot. No values are returned.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR/}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{image.kriging}}, \code{\link{image.krige.bayes}}.  }

\examples{
# See examples in the documentation for image.kriging
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{aplot}% __ONLY ONE__ keyword per line
