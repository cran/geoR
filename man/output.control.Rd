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
               threshold, signal, messages.screen)
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
    predictive distributions are computed and returned.
    If a vector with numbers in the interval
    \eqn{[0,1]}{[0,1]} is provided, the output includes the object
    \code{quantiles}, which contains values of corresponding
    estimated quantiles. For
    example, if \code{quantile = c(0.25, 0.50, 0.75)} the function
    returns the quartiles of the distributions at each of the
    prediction locations.
    If \code{quantile = TRUE} default values equals to \code{c(0.025, 0.5,
      0.975)}, are assumed. 
    A measure of uncertainty of the predictions, analogous to the kriging
    standard error, can be computed by \eqn{(quantile0.975 - quantile0.025)/4}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{threshold}{a (optional) numerical vector.
    If one or more values are provided, object named
    \code{probabilities} is included in the output.
    This object contains, for each
    prediction location, the probability that the variable is less
    than or equal than the threshold provided by the user.
    Defaults to \code{FALSE}. }
  \item{signal}{logical. If \code{TRUE} the signal is predicted, otherwise the
    variable is predicted. If no transformation is performed the
    expectations are the same in both cases and
    the kriging variances are different when the nugget parameter
    is different from zero.  }  
  \item{messages.screen}{logical. Indicates
    whether or not status messages are printed on the output device
    while the function is running. Defaults to \code{TRUE}. }

}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  A list with processed arguments to be passed to the main function.
}
%\references{ ~put references to the literature/web site here ~ }

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}
%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{krige.bayes}} and  \code{\link{krige.conv}}}

%\examples{
%}

\keyword{spatial}% at least one, from doc/KEYWORDS

