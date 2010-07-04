\name{set.coords.lims}

\alias{set.coords.lims}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Sets Limits to Scale Plots}
\description{
  This is an function typically called by functions in the package
  \pkg{geoR} to set limits for the axis when plotting spatial data.
}

\usage{
set.coords.lims(coords, borders = coords, xlim, ylim, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{coords}{an \eqn{n \times 2}{nx2} matrix with coordinates.  }
  \item{borders}{an \eqn{n \times 2}{nx2} matrix with coordinates.  }
  \item{xlim, ylim}{the ranges to be encompassed by the x and y axes.  }
  \item{\dots}{not used, included just for internal handling.  }
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  A \eqn{2 \times 2} matrix with limits for the axis.
}

%\references{ ~put references to the literature/web site here ~ }
\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ }

 
%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

%\examples{}
\keyword{spatial}% at least one, from doc/KEYWORDS

