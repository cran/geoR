\name{boxcox.fit}

\alias{boxcox.fit}
\alias{print.boxcox.fit}
\alias{plot.boxcox.fit}
\alias{lines.boxcox.fit}
\alias{negloglik.boxcox}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Parameter Estimation for the Box-Cox Transformation}

\description{Parameter estimation and plotting of the results for the
  Box-Cox transformed normal distribution.
}

\usage{
boxcox.fit(object, xmat, lambda, lambda2 = NULL, add.to.data = 0, \dots)

print.boxcox.fit(x, \dots)

plot.boxcox.fit(x, hist = TRUE, data = eval(x$call$object), \dots)

lines.boxcox.fit(x, data = eval(x$call$object), \dots)

}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{a vector with the data.  }
  \item{xmat}{a matrix with covariates values. Defaults to \code{rep(1,
      length(y))}.  } 
  \item{lambda}{numerical value(s) for the transformation parameter
    \eqn{\lambda}{lambda}. 
    Used as the initial value in the function for parameter estimation.
    If not provided default values are assumed. If multiple values are passed
    the one with highest likelihood is used as initial value. }
  \item{lambda2}{logical or numerical value(s) of the additional transformation
    (see DETAILS below). Defaults to \code{NULL}.
    If \code{TRUE} this parameter is also estimated and the initial
    value is set to the absolute value of the minimum data.
    A numerical value is provided it is used as the initial
    value. Multiple values are allowed as for \code{lambda}.  }
  \item{add.to.data}{a constant value to be added to the data.  }
  \item{x}{a list, typically an output of the function
    \code{boxcox.fit}.  }
  \item{hist}{logical indicating whether histograms should to be
    plotted.  }
  \item{data}{data values. }
  \item{\dots}{extra parameters to be passed to the minimization
    function \code{\link{optim}} (\code{boxcox.fit}), \code{\link{hist}}
    (\code{plot}) or \code{\link{curve}} (\code{lines}).  }
  
}

\value{
  The functions returns the following results:
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
}

\references{
  Box, G.E.P. and Cox, D.R.(1964) An analysis of transformations. JRSS B
  \bold{26}:211--246.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{rboxcox}} and  \code{\link{dboxcox}} for the
  expression and more on the
  Box-Cox transformation,
  the minimization function \code{\link{optim}},
  the function \code{\link[MASS]{boxcox}} in the package \pkg{MASS} and
  the function \code{\link[car]{boxcox}} in the package \pkg{car}.
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
  boxcox.fit(object = trees[,3], xmat = trees[,1:2])
  }
}

\keyword{regression}
\keyword{models}
\keyword{hplot}

