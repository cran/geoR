\name{sample.prior}

\alias{sample.prior}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Samples from the prior}

\description{
  Sample quadruples \eqn{(\beta, \sigma^2, \phi, \tau^2_{rel})}{(beta,
    sigma^2, phi, tau^2.rel)} from the prior distribution of parameters specifying a Gaussian
  random field. Typically the prior is specified as for when calling \code{\link{krige.bayes}}.
}

\usage{
sample.prior(n, kb.obj=NULL, prior=prior.control())
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{n}{number of samples}
  \item{kb.obj}{on object with an output of \code{\link{krige.bayes}}.  }
  \item{prior}{an call to \code{\link{prior.control}}. Not needed if
    \code{kb.obj} is provided.  }
}
%\details{
%  
%}

\value{
  A data-frame with samples from the model parameters.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.est.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{krige.bayes}} and \code{\link{sample.posterior}}.  }

\examples{
sample.prior(50, prior=prior.control(beta.prior = "normal", beta = .5, beta.var.std=0.1, sigmasq.prior="sc", sigmasq=1.2, df.sigmasq= 2, phi.prior="rec", phi.discrete = seq(0,2, l=21)))
}


\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{distribution}% __ONLY ONE__ keyword per line
