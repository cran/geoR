\name{locations.inside}

\alias{locations.inside}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{select prediction locations inside borders}

\description{
  Selects the prediction locations located inside a polygon defining
  borders of a region where prediction is aimed.
  Typically internally called by \code{\link{krige.bayes}}, \code{\link{krige.conv}}, \code{\link{ksline}}. 
}

\usage{
locations.inside(locations, borders)
}

\arguments{
  \item{locations}{a two columns matrix or dqata frame with coordinates
    of the prediction locations.  }
  \item{borders}{a two column matrix or data-frame with coordinates of a
    polygon defining the borders of the region.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  A two columns matrix with coordinates of points inside the borders.  
}

%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

%\seealso{}

%\examples{
%}
\keyword{spatial}% at least one, from doc/KEYWORDS
