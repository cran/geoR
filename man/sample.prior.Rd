\name{sample.prior}

\alias{sample.prior}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Sample the prior distribution}

\description{
  Sample quadruples \eqn{(\beta, \sigma^2, \phi, \tau^2_{rel})}{(beta,
    sigma^2, phi, tau^2.rel)} from the prior distribution of parameters specifying a Gaussian
  random field.
  Typically the prior is specified in the same manner as when calling \code{\link{krige.bayes}}.
}

\usage{
sample.prior(n, kb.obj=NULL, prior=prior.control())
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{n}{number of samples}
  \item{kb.obj}{on object with an output of \code{\link{krige.bayes}}.  }
  \item{prior}{an call to \code{\link{prior.control}}. Unnecessary if
    \code{kb.obj} is provided.  }
}
%\details{
%  
%}

\value{
  A \eqn{p+3 \times 4}{p+3 x 4} data-frame with a sample of the prior
  distribution of model parameters, where \eqn{p}{p} is the length of
  the mean parameter \eqn{\beta}{beta}.
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

\examples{
sample.prior(50, 
             prior=prior.control(beta.prior = "normal", beta = .5, beta.var.std=0.1, 
                                 sigmasq.prior="sc", sigmasq=1.2, df.sigmasq= 2, 
                                 phi.prior="rec", phi.discrete = seq(0,2, l=21)))
}


\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{distribution}% __ONLY ONE__ keyword per line
