\name{summary.geodata}

\alias{summary.geodata}
\alias{print.summary.geodata}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Summaries for geodata object}
\description{
  Sumarises each of the main elements of an object of the class \code{geodata}.
}
\usage{
summary.geodata(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of the class \code{geodata}.  }
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
}
\references{
  Further information on the package \bold{geoR} can be found at:\cr
  \url{http://www.est.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{summary}}, \code{\link{as.geodata}}. }

\examples{
data(s100)
summary(s100)

data(ca20)
summary(ca20)
}
\keyword{univar}
\keyword{spatial}
