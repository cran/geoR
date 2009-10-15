\name{summary.likGRF}

\alias{summary.likGRF}
\alias{print.summary.likGRF}
\alias{print.likGRF}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Summarizes Parameter Estimation Results for Gaussian Random Fields}

\description{
  Summarizes results returned by the function \code{\link{likfit}}.\cr
  Functions are \emph{methods} for \code{\link[base]{summary}} and
  \code{\link[base]{print}} for the classes \code{likGRF} and \code{summary.likGRF}.
}

\usage{
\method{summary}{likGRF}(object, \dots)
\method{print}{likGRF}(x, digits = max(3, getOption("digits") - 3), \dots)
\method{print}{summary.likGRF}(x, digits = max(3, getOption("digits") - 3), \dots)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of \code{\link[base]{class}} \emph{likGRF}, typically
    a result of a call to \code{\link{likfit}}.  }
  \item{x}{an object of \code{\link[base]{class}} \emph{likGRF} or
    \code{\link[base]{class}} \emph{summary.likGRF}, typically
    resulting from a call to \code{\link{likfit}}.  }
  \item{digits}{the number of significant digits to use when printing. }
  \item{\dots}{extra arguments for \code{\link[base]{print}}.  }
}

\details{
  A detailed summary of a object of the class \code{likGRF} is produced by
  by \code{summary.likGRF} and printed by \code{print.summary.likGRF}.
  This includes model specification with values of fixed and estimated parameters. 
  A simplified summary of the parameter estimation is printed by
  \code{print.likGRF}. 
}

\value{
  \code{print.likGRF} prints the parameter estimates and the value of the
  maximized likelihood.\cr
  \code{summary.likGRF} returns a list with main results of a call to
  \code{\link{likfit}}.\cr
  \code{print.summary.likGRF} prints these results on the screen (or other
  output device) in a "nice" format.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{likfit}}, \code{\link[base]{print}},
  \code{\link[base]{summary}}.
}

\examples{
# See examples for the function likfit()
}

\keyword{spatial}
\keyword{print}

