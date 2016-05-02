\name{locations.inside}

\alias{locations.inside}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Select prediction locations inside borders}

\description{
  Selects the prediction locations located inside a polygon defining
  borders of a region where prediction is aimed.
  Typically internally called by \pkg{geoR} functions
  \code{\link{krige.bayes}}, \code{\link{krige.conv}},
  \code{\link{ksline}}. 
}

\usage{
locations.inside(locations, borders, as.is = TRUE, \ldots)
}

\arguments{
  \item{locations}{a two columns matrix or dqata frame with coordinates
    of the prediction locations.  }
  \item{borders}{a two column matrix or data-frame with coordinates of a
    polygon defining the borders of the region.  }
  \item{as.is}{logical defining if the returned object of of the same
    type (list, data-frame or matrix) as the provided in
    \code{locations}. If \code{FALSE} the function returns a matrix.  }
  \item{\ldots}{arguments to be passed to the internal function
    \code{.geoR_pip} and currently not used. }
}
%\details{
%  This function works basically as a wrapper for \pkg{splancs} functions
%  \code{\link[splancs]{pip}} which in turn uses
%  \code{\link[splancs]{inout}}.
%}

\value{
  A two columns matrix, data-frame or a list with 2 elements with coordinates of points inside the borders.  
}

%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link[sp]{over}},\code{\link[sp]{coordinates}},
    \code{\link[sp]{SpatialPoints}}.
}

\examples{
gr <- pred_grid(parana$borders, by=20)
plot(gr, asp=1, pch="+")
polygon(parana$borders)
gr.in <- locations.inside(gr, parana$borders)
points(gr.in, col=2, pch="+")
}
\keyword{spatial}% at least one, from doc/KEYWORDS
