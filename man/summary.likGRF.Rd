\name{summary.likGRF}

\alias{summary.likGRF}
\alias{print.summary.likGRF}
\alias{print.likGRF}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Summarizes Parameter Estimation Results for Gaussian Random Fields}

\description{
  Summarizes results returned by the function \code{\link{likfit}}.
  Functions are \emph{methods} for classes \code{likGRF} and
  \code{summary.likGRF}.
}

\usage{
summary(obj, ...)
print(obj, digits = "default", ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{obj}{an object of \code{\link{class}} \emph{likGRF}, typically
    a result of a call to \code{\link{likfit}}.  }
  \item{\dots}{extra arguments for \code{\link{print}}. A commonly used
    argument is \code{digits} which
    specifies the number of digits for the numerical output.
    The default is given by \code{options()$digits}.  }
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
  output device) in a "nice" way.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{likfit}}, \code{\link{print}},
  \code{\link{summary}}.
}

\examples{
# See examples for the function likfit()
}

\keyword{spatial}
\keyword{print}

