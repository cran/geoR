\name{BoxCox}

\alias{boxcox.fit}
\alias{boxcox.negloglik}
\alias{rboxcox}
\alias{dboxcox}
\alias{plot.boxcox.fit}
\alias{lines.boxcox.fit}
\alias{print.boxcox.fit}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Box-Cox transformation related functions}

\description{
  Estimate parameters and plot results , computes density and
  random generation for the Box-Cox normal transformed distribution.
}

\usage{
boxcox.fit(data, xmat, lambda, lambda2 = FALSE, add.to.data = 0,...)

print.boxcox.fit(x, ...)

plot.boxcox.fit(x, hist = TRUE, data = eval(x$call$data), ...)

lines.boxcox.fit(x, data = eval(x$call$data), ...)

rboxcox(n, lambda, lambda2 = NULL, mean = 0, sd = 1)

dboxcox(x, lambda, lambda2 = NULL, mean = 0, sd = 1)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{data}{a vector with data.  }
  \item{xmat}{a matrix of covariates. Defaults to \code{rep(1,
      length(y))}.  } 
  \item{lambda}{numerical value(s) of the transformation parameter
    \eqn{lambda}{lambda}. 
    Used as initial value for the parameter estimation function.
    If not provided default values are assumed. If multiple values are passed
    the one with highest likelihood is used as initial value. }
  \item{lambda2}{numerical value(s) of the additional transformation
    (see DETAILS below). Defaults to \code{FALSE}.
    If \code{TRUE} this parameter is also estimated and the initial
    value is set to the absolute value of the minimum data.
    A numerical value provided is used as initial
    value. As for \code{lambda} multiple values are allowed.  }
  \item{add.to.data}{a constant to be added to the data.  }
  \item{n}{number of observations to be generated.  }
  \item{x}{a vector of quantiles (\code{dboxcox}) or an output of
    \code{boxcox.fit} (\code{print}, \code{plot}, \code{lines}).  }
  \item{mean}{a vector of means at for the normal scale.  }
  \item{sd}{a vector of standard deviations at the normal scale.  }
  \item{hist}{logical indicating whether histograms are to be plotted.  }
  \item{\dots}{extra parameters to be passed to the minimization
    function \code{\link{optim}} (\code{boxcox.fit}), \code{\link{hist}}
    (\code{plot}) or \code{\link{curve}} (\code{lines}).  }

}


\details{
  Denote \eqn{Y} the variable at the original scale.
  The Box-Cox transformation is defined by:
  
  \deqn{Y\prime = \left\{ \begin{array}{ll}
    log(Y)
    \mbox{ , if $\lambda = 0$} \cr
    \frac{Y^\lambda - 1}{\lambda} \mbox{ ,  otherwise}
    \end{array} \right.}{Y' = log(Y) if lambda = 0 , Y' = ((Y ^ lambda) - 1)/lambda otherwise }. 
  
  An additional shifting parameter \eqn{lambda_2}{lambda2} can be
  included in which case the transformation is given by: 
  
  \deqn{Y\prime = \left\{ \begin{array}{ll}
    log(Y + \lambda_2)
    \mbox{ , $\lambda = 0$ } \cr
    \frac{(Y + lambda_2)^\lambda - 1}{\lambda} \mbox{ ,  otherwise}
    \end{array} \right.}{Y' = log(Y + lambda2) if lambda = 0 ,
    Y' = (((Y + lambda2) ^ lambda) - 1)/lambda otherwise }.

 The function \code{rboxcox} samples \eqn{Y\prime}{Y'} from the normal distribution using
 \code{\link{rnorm}} and backtransform the values according to the
 equations above. If necessary the back-transformation truncates the values such that
 \eqn{Y\prime \geq \frac{1}{\lambda}}{Y' >= -1/lambda}. This results in
 values of \eqn{Y = 0} in the original scale.
 Increasing the mean and/or reducing the variance might help to avoid truncation.
    
}

\value{

  The functions returns the following results:
  
  \itemize{
    \item{boxcox.fit}{a list with estimated  parameters and results on
      the numerical minimization.  }
    \item{print.boxcox.fit}{print estimated parameters. No values
      returned.  }
    \item{plot.boxcox.fit}{plots histogram of the data (optional) and
      the model.  No values returned.  This function is only valid if
      covariates are not included in \code{boxcox.fit}. }
    \item{lines.boxcox.fit}{adds a line with the fitted model to the
      current plot. No values returned.  This function is only valid if
      covariates are not included in \code{boxcox.fit}.  }
    \item{rboxcox}{a vector with simulated observations.  }
    \item{dboxcox}{a vector of densities.  }
  }
}

\references{
  Box, G.E.P. and Cox, D.R.(1964) An analysis of transformations. JRSS B 26 211-246. 
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{The minimization function \code{\link{optim}},
  \code{\link[MASS]{boxcox}} in the package \code{MASS} and
  \code{\link[car]{boxcox}} in the package \code{car}.
}

\examples{
## Simulating data
simul <- rboxcox(100, lambda=0.5, mean=10, sd=2)
## Finding the ML estimates
ml <- boxcox.fit(simul)
ml
## Ploting histogram and fitted model
plot(ml)
##
## Comparing models with different lambdas,
## zero  means and unit variances
curve(dboxcox(x, lambda=-1), 0, 8)
for(lambda in seq(-.5, 1.5, by=0.5))
  curve(dboxcox(x, lambda), 0, 8, add = TRUE)
##
## Another example, now estimating lambda2
##
simul <- rboxcox(100, lambda=0.5, mean=10, sd=2)
ml <- boxcox.fit(simul, lambda2 = TRUE)
ml
plot(ml)
##
## An example with a regression model
##
if(require(MASS)){
  data(trees)
  boxcox.fit(data = trees[,3], xmat = trees[,1:2])
  }

}

\keyword{distribution}
