\name{loglik.GRF}
\alias{loglik.GRF}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Log-Likelihood for a Gaussian Random Field}
\description{
  This function computes the value of the log-likelihood for a
  realisation of a Gaussian random field.
}
\usage{
loglik.GRF(geodata, coords=geodata$coords, data=geodata$data,
           cov.model="exp", cov.pars, nugget=0, kappa=0.5,
           lambda=1, psiR=1, psiA=0,
           trend="cte", method="ML", compute.dists=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \bold{geoR}
    data-set.
    If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
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
    \code{"1st"} (a first degree polynomial
    on the coordinates), \code{"2nd"} (a second degree polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend). Defaults to \code{"cte"}.  }
  \item{method}{options are \code{"ML"} for likelihood and \code{"REML"} for
    restricted likelihood. Defaults to \code{"ML"}.  }
  \item{compute.dists}{for internal use with other function. Don't
    change the default unless you know what you are doing.  }
}
\details{
  The expression log-likelihood is:
  \deqn{l(\theta) = -\frac{n}{2} \log (2\pi) - \frac{1}{2} \log |\Sigma|
    - \frac{1}{2} (y -
  F\beta)\prime \Sigma^{-1} (y - F\beta),}{l(theta) = -(n/2) * log(2 * pi) - 0.5 * log|V| - 0.5 * (y - F b)' V^{-1} (y - F
    b),}
  where \eqn{n} is the size of the data vector \eqn{y}, \eqn{\beta}{b}
  is the mean (vector) parameter with dimention \eqn{p}, \eqn{\Sigma}{V} is the covariance
  matrix and \eqn{F} is the matrix with the values of the covariates (a
  vector of \eqn{1}'s if the mean is constant.

  The expression restricted log-likelihood is:
  \deqn{rl(\theta) =  -\frac{n-p}{2} \log (2\pi) + \frac{1}{2} \log |F\prime F|
    - \frac{1}{2} \log |\Sigma| - \frac{1}{2} \log |F\prime \Sigma F| - \frac{1}{2} (y -
  F\beta)\prime \Sigma^{-1} (y - F\beta).}{rl(\theta) =  -((n-p)/2) *
  log (2 * pi) + 0.5 * log |F'F| - 0.5 * log |V|  - 0.5 * log |F'VF| -
  0.5 * (y - Fb)' V^(-1) (y - Fb).}
  
}

\value{
  The numerical value of the log-likelihood.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}


\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br},  \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{likfit}} for likelihood-based parameter estimation.}

\examples{
if(is.R()) data(s100)
loglik.GRF(s100, cov.pars=c(0.8, .25), nugget=0.2)
loglik.GRF(s100, cov.pars=c(0.8, .25), nugget=0.2, met="RML")
}

\keyword{spatial}
