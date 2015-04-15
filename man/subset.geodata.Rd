\name{subset.geodata}
\alias{subset.geodata}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Method for subsetting geodata objects}
\description{
  Subsets a object of the class \code{geodata}
  by transforming it to a data-frame, using \code{subset}
  and back transforming to a \code{geodata} object.
}
\usage{
\method{subset}{geodata}(x, ..., other = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of the class \code{geodata}.  }
  \item{\dots}{arguments to be passed to
    \code{\link{subset.data.frame}}.  }
  \item{other}{logical. If \code{TRUE} non-standard
    \code{geodata} elements
    of the original \code{geodata} object are copied to the
    resulting object.  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list which is an object of the class \code{geodata}.  
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ %
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{subset}} for the generic function and methods and
  \code{\link{as.geodata}} for more information on geodata objects.  }
\examples{
subset(ca20, data > 70)
subset(ca20, area == 1)
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{manip}% __ONLY ONE__ keyword per line
