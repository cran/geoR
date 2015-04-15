\name{variog.model.env}

\alias{variog.model.env}
\alias{boot.variofit}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Envelops for Empirical Variograms Based on Model Parameters}

\description{
  Computes envelopes for a empirical variogram by simulating
  data for given model parameters.

  Computes bootstrap paremeter estimates
}

\usage{
variog.model.env(geodata, coords = geodata$coords, obj.variog,
                 model.pars, nsim = 99, save.sim = FALSE, messages)

boot.variofit(geodata, coords = geodata$coords, obj.variog,
              model.pars, nsim = 99, trace = FALSE, messages)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing element \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \pkg{geoR} data-set.
    If not provided the argument
    \code{coords} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing
    Euclidean 
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{obj.variog}{an object of the class \code{"variogram"}, typically
    an output of the function \code{\link{variog}}.  }
  \item{model.pars}{a list with model specification and parameter
    values. The input is typically an object of the class
    \code{variomodel} which is an output of \code{\link{likfit}},
    \code{\link{variofit}}.
    The required components of the list are:
    \itemize{
      \item \code{beta}, the mean parameter. Defaults to zero.
      \item \code{cov.model}, the covariance model. Defaults to
      "exponential".
      \item \code{cov.pars}, the covariance parameters
      \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}.
      \item \code{kappa}, the extra covariance parameters for some of the
      covariance models. Defaults to \eqn{0.5}.
      \item \code{nugget}, the error component variance. Defaults to zero.
      \item \code{estimator.type}, the type of variogram
      estimator. Options for \code{"classical"} and
      \code{"robust"}. Defaults to \code{obj.variog$estimator}.  }
  } 
  \item{nsim}{number of simulations used to compute the envelopes.
    Defaults to 99.  }
  \item{save.sim}{logical. Indicates whether or not the simulated data
    are included in the output. Defaults to \code{FALSE}.  }
  \item{trace}{logical. If \code{TRUE} the fitted values for the
    bootstrap parameter estimation are printend while the function is
    running.  }
  \item{messages}{logical. If \code{TRUE}, the default, status messages
    are printed while the function is running.  }
}

\details{
  The envelopes are computed assuming a (transformed) Gaussian random
  field model.
  Simulated values are generated at the data locations, given the model parameters.
  The empirical variogram is computed for each simulation
  using the same lags as for the original variogram of the data.
  The envelopes are computed by taking, at each lag, the maximum and minimum
  values of the variograms for the simulated data.
}

\value{
  An object of the \code{\link{class}} \code{"variogram.envelope"} which is a
  list with the components:
  \item{u}{a vector with distances. }
  \item{v.lower}{a vector with the minimum variogram values at each
    distance in \code{u}.  }
  \item{v.upper}{a vector with the maximum variogram values at each
    distance in \code{u}.  }
  \item{simulations}{a matrix with the simulated data.
    Only returned if \code{save.sim = TRUE}.  }
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br},  \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{variog.mc.env}} for envelops computed by
  using data permutation, 
  \code{\link{variog}} for variogram calculations,
  \code{\link{plot.variogram}} and \code{\link{variog.mc.env}} for
  graphical output. The functions
  \code{\link{likfit}}, \code{\link{variofit}}
  are used to estimate the model parameters.
}

\examples{
s100.ml <- likfit(s100, ini = c(0.5, 0.5), fix.nugget = TRUE)
s100.vario <- variog(s100, max.dist = 1)
s100.env <- variog.model.env(s100, obj.v = s100.vario,
                             model.pars = s100.ml)
plot(s100.vario, env = s100.env)
}

\keyword{spatial}%-- one or more ...

