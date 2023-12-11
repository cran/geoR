\name{boxcox}

\alias{rboxcox}
\alias{dboxcox}

\concept{Box-Cox transformation} 

%- Also NEED an `\alias' for EACH other topic documented here.

\title{The Box-Cox Transformed Normal Distribution}

\description{Functions related with the Box-Cox family of transformations.
  Density and random generation for the Box-Cox transformed normal
  distribution with mean
  equal to \code{mean} and standard deviation equal to \code{sd}, \emph{in the normal scale}.
}

\usage{
rboxcox(n, lambda, lambda2 = NULL, mean = 0, sd = 1)

dboxcox(x, lambda, lambda2 = NULL, mean = 0, sd = 1)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{lambda}{numerical value(s) for the transformation parameter
    \eqn{\lambda}{lambda}.  }
  \item{lambda2}{logical or numerical value(s) of the additional transformation
    (see DETAILS below). Defaults to \code{NULL}.  }
  \item{n}{number of observations to be generated.  }
  \item{x}{a vector of quantiles (\code{dboxcox}) or an output of
    \code{boxcoxfit} (\code{print}, \code{plot}, \code{lines}).  }
  \item{mean}{a vector of mean values at the normal scale.  }
  \item{sd}{a vector of standard deviations at the normal scale.  }
}

\details{
  Denote \eqn{Y} the variable at the original scale and \eqn{Y'} the
  transformed variable. The Box-Cox transformation is defined by:
  
  \deqn{Y' = \left\{ \begin{array}{ll}
    log(Y)
    \mbox{ , if $\lambda = 0$} \cr
    \frac{Y^\lambda - 1}{\lambda} \mbox{ ,  otherwise}
    \end{array} \right.}{Y' = log(Y) if lambda = 0 , Y' = ((Y ^ lambda) - 1)/lambda otherwise }. 
  
  An additional shifting parameter \eqn{\lambda_2}{lambda2} can be
  included in which case the transformation is given by: 
  
  \deqn{Y' = \left\{
    \begin{array}{ll}
    log(Y + \lambda_2)
    \mbox{ , $\lambda = 0$ } \cr
    \frac{(Y + \lambda_2)^\lambda - 1}{\lambda} \mbox{ ,  otherwise}
    \end{array} \right.}{Y' = log(Y + lambda2) if lambda = 0 ,
    Y' = (((Y + lambda2) ^ lambda) - 1)/lambda otherwise }.

 The function \code{rboxcox} samples \eqn{Y'}{Y'} from the normal distribution using
 the function \code{\link{rnorm}} and backtransform the values according to the
 equations above to obtain values of \eqn{Y}{Y}.
 If necessary the back-transformation truncates the values such that
 \eqn{Y' \geq \frac{1}{\lambda}}{Y' >= -1/lambda} results in
 \eqn{Y = 0}{Y = 0} in the original scale.
 Increasing the value of the mean and/or reducing the variance might help to avoid truncation.
}

\value{
  The functions returns the following results:
  \item{rboxcox}{a vector of random deviates.  }
  \item{dboxcox}{a vector of densities.  }
}

\references{
  Box, G.E.P. and Cox, D.R.(1964) An analysis of transformations. JRSS B
  \bold{26}:211--246.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{The parameter estimation function is \code{\link{boxcoxfit}}.
  Other packages has BoxCox related functions such as \code{\link[MASS]{boxcox}} in the package \pkg{MASS} and
  the function \code{box.cox} in the package \samp{car}.
%  the function \code{\link[car]{box.cox}} in the package \pkg{car}.
}

\examples{
## Simulating data
simul <- rboxcox(100, lambda=0.5, mean=10, sd=2)
##
## Comparing models with different lambdas,
## zero  means and unit variances
curve(dboxcox(x, lambda=-1), 0, 8)
for(lambda in seq(-.5, 1.5, by=0.5))
  curve(dboxcox(x, lambda), 0, 8, add = TRUE)
}
\keyword{distribution}
