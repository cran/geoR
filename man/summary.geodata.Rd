\name{summary.geodata}

\alias{summary.geodata}
\alias{print.summary.geodata}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summaries for geodata object}
\description{
  Sumarises each of the main elements of an object of the class \code{geodata}.
}
\usage{
\method{summary}{geodata}(object, lambda =1, add.to.data = 0,
                by.realisations=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of the class \code{geodata}.  }
  \item{lambda}{value of the Box-Cox transformation parameter. Two particular cases
    are \eqn{\lambda = 1}{lambda = 1} which corresponds to no
    transformation and  \eqn{\lambda = 0}{lambda = 0} corresponding to
    the log-transformation.  }
  \item{add.to.data}{scalar, Constant value to be added to the data
    values.
    Only used if a value different from 1 is passed to the argument \code{lambda}.}
  \item{by.realisations}{logical. Indicates whether the summary must be performed separatly for each realisation, if the \code{geodata} object contains the element \code{realisations}. Defaults to \code{TRUE}.  }
\item{\dots}{further arguments to be passed to the function
    \code{\link{summary.default}}.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  A list with components
  \item{coords.summary }{a matrix with minimum and maximum values for
    the coordinates.  }
    \item{distances.summary}{minimum and maximum distances between pairs of points.  }
  \item{borders.summary }{a matrix with minimum and maximum values for
    the coordinates. Only returned if there is an element \code{borders}
  in the \code{geodata} object.  }
  \item{data.summary }{summary statistics (min, max, quartiles and mean)
    for the data.  }
  \item{units.m.summary }{summary statistics (min, max, quartiles and mean)
    for the offset variable. Only returned if there is an element \code{units.m}
  in the \code{geodata} object. }
  \item{covariate.summary }{summary statistics (min, max, quartiles and mean)
    for the covariate(s). Only returned if there is an element \code{covariate}
    in the \code{geodata} object. }
    \item{others}{names of other elements if present in the \code{geodata} object.  }
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
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{summary}}, \code{\link{as.geodata}}. }

\examples{
summary(s100)

summary(ca20)
}
\keyword{univar}
\keyword{spatial}
