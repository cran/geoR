\name{varcovBGCCM}

\alias{varcovBGCCM}
\alias{.dist12}
\alias{.cov012.model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Covariance matrix for the bivariate Gaussian
  common component geostatistical model}
\description{
  Covariance matrix for the bivariate Gaussian
  common component geostatistical model or its inverse,
  and optionally the determinant of the matrix.
}
\usage{
varcovBGCCM(dists.obj, cov0.pars, cov1.pars, cov2.pars,
            cov0.model = "matern", cov1.model = "matern",
            cov2.model = "matern", kappa0 = 0.5, kappa1 = 0.5,
            kappa2 = 0.5, scaled = TRUE, inv = FALSE, det = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dists.obj}{ ~~Describe \code{dists.obj} here~~ }
  \item{cov0.pars}{ ~~Describe \code{cov0.pars} here~~ }
  \item{cov1.pars}{ ~~Describe \code{cov1.pars} here~~ }
  \item{cov2.pars}{ ~~Describe \code{cov2.pars} here~~ }
  \item{cov0.model}{ ~~Describe \code{cov0.model} here~~ }
  \item{cov1.model}{ ~~Describe \code{cov1.model} here~~ }
  \item{cov2.model}{ ~~Describe \code{cov2.model} here~~ }
  \item{kappa0}{ ~~Describe \code{kappa0} here~~ }
  \item{kappa1}{ ~~Describe \code{kappa1} here~~ }
  \item{kappa2}{ ~~Describe \code{kappa2} here~~ }
  \item{scaled}{ ~~Describe \code{scaled} here~~ }
  \item{inv}{logical. If \code{TRUE} the inverse of the covariance
    matrix
    is returned instead.  }
  \item{det}{logical. Optional, if \code{TRUE} the logarithm of the
  detarminant of the covariance matrix is returned as an attribute.  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A matrix which is the covariance matrix for the
  bivariate Gaussian
  common component geostatistical model or its inverse if
  \code{inv=TRUE}.
  If \code{det=T} the logarithm of the determinant
  of the matrix is also returned as an attribute
  named \code{logdetS}. 
}
%\references{ ~put references to the literature/web site here ~ }

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\section{Warning  }{This is a new function and still in draft format and pretty much untested.}

\seealso{\code{\link{cov.spatial}}, \code{\link{varcov.spatial}}}
\examples{
# see http://www.est.ufpr.br/geoR/tutorials/CCM.R
}
\keyword{spatial}% at least one, from doc/KEYWORDS
