\name{sample.posterior}

\alias{sample.posterior}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Samples from the posterior distribution}

\description{
  Sample quadruples \eqn{(\beta, \sigma^2, \phi, \tau^2_{rel})}{(beta,
    sigma^2, phi, tau^2.rel)} from the posterior
  distribution returned by \code{\link{krige.bayes}}.
}

\usage{
sample.posterior(n, kb.obj)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{n}{number of samples}
  \item{kb.obj}{on object with an output of \code{\link{krige.bayes}}.  }
}

%\details{
%  
%}

\value{
  A \eqn{n \times 4}{n x 4} data-frame with samples from the posterior distribution of the model parameters.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{krige.bayes}} and \code{\link{sample.posterior}}.  }

%\examples{
%}

\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{distribution}% __ONLY ONE__ keyword per line
