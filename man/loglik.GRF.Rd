\name{loglik.GRF}

\alias{loglik.GRF}
%- Also NEED an `\alias' for EACH other topic documented here.

\concept{Likelihood function for a Gaussian Random Field}

\title{Log-Likelihood for a Gaussian Random Field}

\description{
  This function computes the value of the log-likelihood for a
  Gaussian random field.
}

\usage{
loglik.GRF(geodata, coords = geodata$coords, data = geodata$data,
           obj.model = NULL, cov.model = "exp", cov.pars, nugget = 0,
           kappa = 0.5, lambda = 1, psiR = 1, psiA = 0,
           trend = "cte", method.lik = "ML", compute.dists = TRUE,
           realisations = NULL)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \pkg{geoR}
    data-set.
    If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
  \item{obj.model}{a object of the class \code{variomodel} with a fitted
    model. Tipically an output of 
    \code{\link{likfit}} or \code{\link{variofit}}. }
  \item{cov.model}{a string specifying the model for the correlation
    function. For further details see 
    documentation for \code{\link{cov.spatial}}. }
  \item{cov.pars}{a vector with 2 elements with values of  the covariance parameters
    \eqn{\sigma^2}{sigma^2} (partial sill) and \eqn{\phi}{phi} (range
    parameter).  }
  \item{nugget}{value of the nugget parameter. Defaults to \eqn{0}.  }
  \item{kappa}{value of the smoothness parameter. Defaults to
    \eqn{0.5}.  }
  \item{lambda}{value of the Box-Cox tranformation parameter. Defaults
    to \eqn{1}.  }
  \item{psiR}{value of the anisotropy ratio parameter. Defaults to
    \eqn{1}, corresponding to isotropy.  }
  \item{psiA}{value (in radians) of the anisotropy rotation parameter. Defaults to zero.}
  \item{trend}{specifies the mean part of the model.
    The options are:
    \code{"cte"} (constant mean),
    \code{"1st"} (a first order polynomial
    on the coordinates), \code{"2nd"} (a second order polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend). Defaults to \code{"cte"}.  }
  \item{method.lik}{options are \code{"ML"} for likelihood and \code{"REML"} for
    restricted likelihood. Defaults to \code{"ML"}.  }
  \item{compute.dists}{for internal use with other function. Don't
    change the default unless you know what you are doing.  }
  \item{realisations}{optional. A vector indicating replication number
    for each data. For more details see \code{\link{as.geodata}}. }
}

\details{
  The expression log-likelihood is:
  \deqn{l(\theta) = -\frac{n}{2} \log (2\pi) - \frac{1}{2} \log |\Sigma|
    - \frac{1}{2} (y -
  F\beta)' \Sigma^{-1} (y - F\beta),}{l(theta) = -(n/2) * log(2 * pi) - 0.5 * log|V| - 0.5 * (y - F b)' V^{-1} (y - F
    b),}
  where \eqn{n} is the size of the data vector \eqn{y}, \eqn{\beta}{b}
  is the mean (vector) parameter with dimention \eqn{p}, \eqn{\Sigma}{V} is the covariance
  matrix and \eqn{F} is the matrix with the values of the covariates (a
  vector of \eqn{1}'s if the mean is constant.

  The expression restricted log-likelihood is:
  \deqn{rl(\theta) =  -\frac{n-p}{2} \log (2\pi) + \frac{1}{2} \log |F' F|
    - \frac{1}{2} \log |\Sigma| - \frac{1}{2} \log |F' \Sigma F| - \frac{1}{2} (y -
  F\beta)' \Sigma^{-1} (y - F\beta).}{rl(\theta) =  -((n-p)/2) *
  log (2 * pi) + 0.5 * log |F'F| - 0.5 * log |V|  - 0.5 * log |F'VF| -
  0.5 * (y - Fb)' V^(-1) (y - Fb).}
  
}

\value{
  The numerical value of the log-likelihood.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}


\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br},  \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{likfit}} for likelihood-based parameter estimation.}

\examples{
loglik.GRF(s100, cov.pars=c(0.8, .25), nugget=0.2)
loglik.GRF(s100, cov.pars=c(0.8, .25), nugget=0.2, met="RML")

## Computing the likelihood of a model fitted by ML
s100.ml <- likfit(s100, ini=c(1, .5))
s100.ml
s100.ml$loglik
loglik.GRF(s100, obj=s100.ml)

## Computing the likelihood of a variogram fitted model
s100.v <- variog(s100, max.dist=1)
s100.vf <- variofit(s100.v, ini=c(1, .5))
s100.vf
loglik.GRF(s100, obj=s100.vf)
}

\keyword{spatial}
