\name{practicalRange}
\alias{practicalRange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pratical range for correlation functions}
\description{
  Computes practical ranges for the correlation functions implemented in
  the geoR package
}
\usage{
practicalRange(cov.model, phi, kappa = 0.5, correlation = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cov.model}{correlation model as documented in
    \code{\link{cov.spatial}}.  }
  \item{phi}{correlation parameter as documented in \code{\link{cov.spatial}}. }
  \item{kappa}{additional correlation parameter as documented in
    \code{\link{cov.spatial}}.  }
  \item{correlation}{correlation threshold for asymptotic
    models. Defaults to 0.05.  }
  \item{\dots}{arguments to be passed to \code{\link{optimise}}. }
}
%\details{
% 
%}
\value{
  A scalar with the value of the practical range.
}
%\references{ ~put references to the literature/web site here ~ }
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}
\seealso{\code{\link{cov.spatial}}}
\examples{
practicalRange("exp", phi=10)
practicalRange("sph", phi=10)
practicalRange("gaus", phi=10)
practicalRange("matern", phi=10, kappa=0.5)
practicalRange("matern", phi=10, kappa=1.5)
practicalRange("matern", phi=10, kappa=2.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
%\keyword{}% __ONLY ONE__ keyword per line
