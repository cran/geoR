\name{output.control}

\alias{output.control}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{ Defines output options for prediction functions }

\description{
  Auxiliary function defining output options for
  \code{\link{krige.bayes}} and \code{\link{krige.conv}}.
}

\usage{
output.control(n.posterior, n.predictive, moments, n.back.moments,
               simulations.predictive, mean.var, quantile,
               threshold, sim.means, sim.vars, signal, messages)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{n.posterior}{number of samples to be taken from the posterior
    distribution. Defaults to 1000. }
  \item{n.predictive}{number of samples to be taken from the
    predictive distribution. Default equals to
    \code{n.posterior}.  }
  \item{moments}{logical. Indicates whether the moments of the
    predictive distribution are returned. If \code{lambda = 1}
    there is no transformation/back-transformation.
    If \code{lambda = 0} or \code{lambda = 0.5} the moments are
    back-transformed by analytical expressions. For other cases the
    back-transformation is done by simulation. Defaults to \code{TRUE}. }
  \item{n.back.moments}{number of sample to back-transform moments by
    simulation. Defaults to 1000.  }
  \item{simulations.predictive}{logical. Defines whether to draw simulations
    from the predictive distribution.
    Only considered if prediction
    locations are provided in the argument \code{locations} of the main functions.
    Defaults to \code{FALSE} but changed to
    \code{TRUE} if an integer greater then zero is provided in the
    argument \code{n.predictive} and/or simulations are required in
    order to compute quantities required by other arguments such as
    threshold, quantiles and some values of the transformation parameter.
  } 
  \item{mean.var}{logical (optional). Indicates whether mean and
    variances of the simulations of the predictive distributions are
    computed and returned.   } 
  \item{quantile}{a (optional) numeric vector.
    If provided indicates whether quantiles of the
    simulations from the
    predictive distribution are computed and returned.
    If a vector with numbers in the interval
    \eqn{[0,1]}{[0,1]} is provided, the output includes the object
    \code{quantiles}, which contains values of corresponding
    estimated quantiles. For
    example, if \code{quantile = c(0.25, 0.50, 0.75)} the function
    returns the quartiles of the predictive distributions at each of the
    prediction locations.
    If \code{quantile = TRUE} default values \code{c(0.025, 0.5,
      0.975)} are assumed. 
    A measure of uncertainty of the predictions, an alternative to the kriging
    standard error, computed by \eqn{(quantile_0.975 - quantile_0.025)/4}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{threshold}{Optional. A numerical vector.
    If one or more values are provided, an object named
    \code{probabilities} is included in the output.
    This object contains, for each
    prediction location, the probability that the variable is less
    than or equal than the threshold provided by the user.
    Defaults to \code{FALSE}. }
  \item{sim.means}{logical (optional). Indicates whether mean
    of each of the conditional simulations of the predictive
    distribution should be computed and returned. Defaults to
  \code{TRUE}, if simulations from the predictive are required.  } 
  \item{sim.vars}{logical (optional). Indicates whether variance
    of each of the conditional simulations of the predictive
    distribution should be computed and returned. Defaults to \code{FALSE}.  } 
  \item{signal}{logical indicating whether the signal or the variable is
    to be predicted. Different defaults are set internally by
    functions calling \code{output.control}. See \code{DETAILS}
    below.  }
  \item{messages}{logical. Indicates
    whether or not status messages are printed on the output device
    while the function is running. Defaults to \code{TRUE}. }

}

\details{

  \bold{SIGNAL}

  This function is typically called by the \pkg{geoR}'s prediction functions
  \code{\link{krige.bayes}} and \code{\link{krige.conv}}
  defining the output.

  By default, \code{\link{krige.bayes}} sets \code{signal = TRUE}
  and \code{\link{krige.conv}} sets \code{signal = FALSE}.

  The underlying model
  \deqn{Y(x) = \mu + S(x) + \epsilon}{Y(x) = mu + S(x) + e}
  assumes that observations \eqn{Y(x)}{Y(x)} are noisy
  versions of a \emph{signal} \eqn{S(x)}{S(x)} and
  \eqn{Var(\epsilon)=\tau^2}{Var(e)=tau^2} is the nugget variance.

  If \eqn{\tau^2 = 0}{tau^2 = 0} the \eqn{Y} and \eqn{S} are
  indistiguishable.\cr

  If \eqn{\tau^2 > 0}{tau^2 > 0} and regarded as measurement error, the
  option \code{signal} defines whether the \eqn{S} (\code{signal =
    TRUE}) or the variable \eqn{Y} (\code{signal = FALSE}) is to be
  predicted.\cr
  For the latter the predictions will "honor" the data,
  i.e. predicted values will coincide with the data, at data locations.\cr
  For unsampled locations and untransformed data,
  the predicted values equals data
  regardless \code{signal = TRUE} or \code{FALSE}, however
  predictions variances will differ.

  The function \code{\link{krige.conv}} has an argument
  \code{micro.scale}. If \eqn{micro.scale > 0} the error term is
  divided as \eqn{\epsilon = \epsilon_{ms} + \epsilon_{me}}{e = e_{ms} +
    e_{me}} and the nugget variance is divided into two terms: \emph{micro-scale variance}
  and \emph{measurement error}.\cr
  If \code{signal = TRUE} the term \eqn{\epsilon_{ms}}{e_{ms}} is
  regarded as part of the signal and consequently the \emph{micro-scale variance} is added to
  the prediction variance.\cr
  If \code{signal = FALSE} the total error variance \eqn{\tau^2}{tau^2}
  is added to the prediction variance.
}

\value{
  A list with processed arguments to be passed to the main function.
}
%\references{ ~put references to the literature/web site here ~ }

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{The prediction functions \code{\link{krige.bayes}} and  \code{\link{krige.conv}}.}

%\examples{
%}

\keyword{spatial}% at least one, from doc/KEYWORDS

