\name{likfitBGCCM}
\alias{likfitBGCCM}
\alias{.negloglikBGCCM}
\alias{loglikBGCCM}
\alias{.naiveLL.BGCCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits the bivariate Gaussian
  common component geostatistical model}
\description{Computes maximum likelihood estimates of the 
  bivariate Gaussian
  common component geostatistical model.
}

\usage{
likfitBGCCM(geodata1, geodata2, ini.sigmasq, ini.phi,
            cov0.model="matern", cov1.model="matern", cov2.model="matern",
            kappa0=0.5, kappa1=0.5, kappa2=0.5,
            fc.min = c("optim", "nlminb"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geodata1}{an object of the class \code{geodata} with the first variable. }
  \item{geodata2}{an object of the class \code{geodata} with the second variable. }
  \item{ini.sigmasq}{optional, a vector with initial values for the
        variance parameters. If not provided default values are used.  }
  \item{ini.phi}{optional, a vector with initial values for the
    correlation range parameters. If not provided default values are
    used. }
  \item{cov0.model, cov1.model, cov2.model}{covariance model for each of
    the processes. See \code{\link{cov.spatial}} for details.  }
  \item{kappa0, kappa1, kappa2}{extra parameter for some covariance
    models.  }
  \item{fc.min}{a string indication which function
  should be used to minimise the negative of the log-likelihood. }
\item{\dots}{further arguments to be passed to \code{\link{optim}}
  or \code{\link{nlminb}}.  }
}
%\details{
%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A list with model fitting information to which
  the class \code{BGCCM} is assigned.
  \item{mu}{a 2 elements vector with mean estimates.}
  \item{sigmasq}{a 4 elements vector with variance estimates.}
  \item{phi}{a 3 elements vector with estimated correlation
    parameters values.  }
  \item{loglik}{a scalar. Maximised value of the log-likelihood.  }
  \item{optim}{results returned by \code{\link{optim}} or
    \code{\link{nlminb}}. }
  \item{\dots}{and other information related to the model fitting. }
}
%\references{ ~put references to the literature/web site here ~ }
\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ 
%
% ~Make other sections like Warning with \section{Warning }{....} ~
%}

\section{Warning  }{ This is a new function and still in draft format and pretty much untested.}

\seealso{\code{\link{optim}}, \code{\link{nlminb}},
  \code{\link{varcovBGCCM}},
  \code{\link{as.geodata}}, \code{\link{likfit}}. }
\examples{
# see http://www.leg.ufpr.br/geoR/tutorials/CCM.R
}
\keyword{spatial}% at least one, from doc/KEYWORDS
\keyword{models}
