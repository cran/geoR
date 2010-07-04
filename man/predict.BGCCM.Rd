\name{predict.BGCCM}
\alias{predict.BGCCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction for the bivariate Gaussian common component
geostatistical model}
\description{
  Performs prediction for the bivariate Gaussian common component
  geostatistical model
}
\usage{
\method{predict}{BGCCM}(object, locations, borders,
              variable.to.predict = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{on object of the class \code{BGCCMfit}, which is
    an output of \code{\link{likfitBGCCM}}.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data-frame with
    the 2-D coordinates of the \eqn{N} prediction locations,
    or a list for which
    the first two components are used. Input is internally checked by the
    function \code{check.locations}.  }
  \item{borders}{optional. If missing, by default reads the element
    \code{borders} of the \code{geodata} object of the variable to be
    predicted. Ignored if set to \code{NULL}.
    If a two column matrix defining a polygon is
    provided the prediction is performed only at locations inside this polygon. }
  \item{variable.to.predict}{scalar with options for values or 2
    indicating which variable is to be predicted. }
  \item{\dots}{not yet used.  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list of the class \code{BGCCMpred} with components:
  \item{predicted}{predicted values.  }
  \item{krige.var}{prediction variances.  }
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\section{Warning  }{ This is a new function and still in draft format and pretty much untested.}

\seealso{\code{\link{likfitBGCCM}}}
\examples{
# see http://www.leg.ufpr.br/geoR/tutorials/CCM.R
}
\keyword{spatial}% at least one, from doc/KEYWORDS
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

