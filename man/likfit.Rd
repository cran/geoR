\name{likfit}

\alias{likfit}
\alias{likfit.limits}
\alias{negloglik.GRF}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Likelihood Based Parameter Estimation for Gaussian Random Fields}

\description{
  \emph{Maximum likelihood} (ML) or
  \emph{restricted maximum likelihood} (REML) parameter estimation for (transformed) Gaussian
  random fields.
}

\usage{
likfit(geodata, coords = geodata$coords, data = geodata$data,
       trend = "cte", ini.cov.pars, fix.nugget = FALSE, nugget = 0,
       fix.kappa = TRUE, kappa = 0.5, fix.lambda = TRUE, lambda = 1,
       fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
       cov.model = "matern", realisations = NULL,
       method.lik = "ML", components = FALSE,
       nospatial = TRUE, limits = likfit.limits(),
       print.pars = FALSE, messages.screen = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \bold{geoR}
    data-set.
    If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix where each row has the 2-D
    coordinates of the \eqn{n} data locations.
    By default it takes the
    component \code{coords} of the argument \code{geodata}, if provided.
  }
  \item{data}{a vector with \emph{n} data values. By default it takes the
    component \code{data} of the argument \code{geodata}, if provided.  }
  \item{trend}{specifies the mean part of the model. See documentation
    of \code{\link{trend.spatial}} for further details.
    Defaults to \code{"cte"}.  }
  \item{ini.cov.pars}{initial values for the covariance parameters:
    \eqn{\sigma^2}{sigma^2} (partial sill) and \eqn{\phi}{phi} (range
    parameter).  Typically a vector with two components. However a
    matrix can be used to provide several initial values. See
    \code{DETAILS} below.}
  \item{fix.nugget}{logical, indicating whether the parameter
    \eqn{\tau^2}{tau^2} (nugget variance) should be regarded as fixed
    (\code{fix.nugget = TRUE}) or should be
    estimated (\code{fix.nugget = FALSE}). Defaults to
    \code{FALSE}.  }    
  \item{nugget}{value of the nugget parameter.
    Regarded as a
    fixed value if \code{fix.nugget = TRUE} otherwise
    as the initial value for the 
    minimization algorithm. 
    Defaults to zero.  }
  \item{fix.kappa}{logical, indicating whether the extra parameter
    \eqn{\kappa}{kappa} should be regarded as fixed
    (\code{fix.kappa = TRUE}) or should be 
    estimated (\code{fix.kappa = FALSE}). Defaults to
    \code{TRUE}.  }    
  \item{kappa}{value of the extra parameter \eqn{\kappa}{kappa}.
    Regarded as a fixed value if \code{fix.kappa = TRUE}
    otherwise as the initial value for the
    minimization algorithm. Defaults to
    \eqn{0.5}. This parameter is valid only if the covariance function is one
    of: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} or
    \code{"gneiting.matern"}.  For more details on covariance functions
    see documentation for \code{\link{cov.spatial}}.  }
  \item{fix.lambda}{logical, indicating whether the Box-Cox transformation parameter
    \eqn{\lambda}{lambda} should be regarded as fixed
    (\code{fix.lambda = TRUE}) or should be
    be estimated (\code{fix.lambda = FALSE}). Defaults to
    \code{TRUE}.  }    
  \item{lambda}{value of the Box-Cox transformation parameter
    \eqn{\lambda}{lambda}.
    Regarded as a fixed value if \code{fix.lambda = TRUE} otherwise 
    as the initial value for the
    minimization algorithm. Defaults to
    \eqn{1}. Two particular cases are \eqn{\lambda = 1}{lambda = 1}
    indicating no transformation and \eqn{\lambda = 0}{lambda
      = 0} indicating log-transformation.  }
  \item{fix.psiA}{logical, indicating whether the anisotropy angle parameter
    \eqn{\psi_R}{psi_R} should be regarded as fixed
    (\code{fix.psiA = TRUE}) or should 
    be estimated (\code{fix.psiA = FALSE}). Defaults to
    \code{TRUE}.  }    
  \item{psiA}{value (in radians) for the anisotropy angle parameter
    \eqn{\psi_A}{psi_A}. 
    Regarded as a fixed value if \code{fix.psiA = TRUE}
    otherwise  as the initial value for the
    minimization algorithm.
    Defaults to \eqn{0}. See \code{\link{coords.aniso}} for further
    details on anisotropy correction.  }
  \item{fix.psiR}{logical, indicating whether the anisotropy ratio parameter
    \eqn{\psi_R}{psi_R} should be regarded as fixed
    (\code{fix.psiR = TRUE}) or should be 
    estimated (\code{fix.psiR = FALSE}). Defaults to
    \code{TRUE}.  }    
  \item{psiR}{value, always greater than 1, for the anisotropy ratio  parameter
    \eqn{\psi_R}{psi_R}. 
    Regarded as a fixed value if \code{fix.psiR = TRUE}
    otherwise as the initial value for the
    minimization algorithm.
    Defaults to \eqn{1}. See \code{\link{coords.aniso}} for further
    details on anisotropy correction.  }
  \item{cov.model}{a string specifying the model for the correlation
    function. For further details see 
    documentation for \code{\link{cov.spatial}}.
    Defaults are equivalent to the \emph{exponential} model.  }
  \item{realisations}{optional. A vector indicating replication number
    for each data. For more details see \code{\link{as.geodata}}. }
  \item{method.lik}{options are \code{"ML"} for maximum likelihood and \code{"REML"} for
    restricted maximum likelihood. Defaults to \code{"ML"}.  }
  \item{components}{an \eqn{n \times 3}{n x 3} data-frame with fitted
    values for the three model components: trend, spatial and residuals.
    See the section \code{DETAILS} below for the model specification.}
  \item{nospatial}{logical. If \code{TRUE} parameter estimates for the
    model without spatial component are included in the output.  }
  \item{limits}{values defining lower and upper limits for the model
    parameters used in the numerical minimization.
    The auxiliary function \code{likfit.limits()} is called to set the
    limits.  }
  \item{print.pars}{logical. If \code{TRUE} the parameters and the value
  of the negative log-likelihood (up to a constant) are printed each
  time the function to be minimised is called.  }
  \item{messages.screen}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{additional parameters to be passed to the minimization
    function. Typically  arguments of the type  \code{control()} which controls the
    behavior of the minimization algorithm. For further details see documentation 
    for the minimization function \code{\link{optim}}.  }
}

\details{
  This function estimate the parameters of the Gaussian random field model,
  specified here by:
  \deqn{Y(x) = \mu(x) + S(x) + e}{Y(x) = mu(x) + S(x) + e}
  where
  \itemize{
    \item \eqn{x} defines a spatial location. Typically Euclidean
    coordinates on a plane. 
    \item \eqn{Y} is the variable been observed.
    \item \eqn{\mu(x) = X\beta}{mu(x) = X \%*\% beta} is the mean component of
    the model (trend). 
    \item \eqn{S(x)} is a stationary Gaussian process with variance \eqn{\sigma^2}{sigma^2}
    (partial sill) and a correlation function parametrized by
    \eqn{\phi}{phi} (the range parameter). Possible extra parameters
    for the correlation function are the smoothness parameter
    \eqn{\kappa}{kappa} and the anisotropy parameters \eqn{\phi_R}{phi_R} and
    \eqn{\phi_A}{phi_A} (anisotropy ratio and angle, respectively).  
    \item \eqn{e} is the error term with variance parameter  \eqn{\tau^2}{tau^2}
    (nugget variance).
  }
  
  The additional parameter \eqn{\lambda}{lambda} allows the Box-Cox
  transformation. If used \eqn{Y(x)} above is replaced by \eqn{g(Y(x))}
  such that
  \deqn{g(Y(x)) = \frac{Y^\lambda(x) - 1}{\lambda}.}{g(Y(x)) = ((Y^lambda(x)) -
    1)/lambda .}

  Two cases of particular interest are \eqn{\lambda = 1}{lambda = 1}
  indicating no transformation and \eqn{\lambda = 0}{lambda = 0}
  indicating log-transformation.
  
  Parameter estimation is performed numerically using the \R
  function \code{\link{optim}} to minimize the
  negative log-likelihood computed by \code{negloglik.GRF}.
  
  Lower and upper limits for parameter values can be
  individually specified using the function \code{likfit.limits()}.
  For example, including the following in the function call:\cr
  \code{limits = likfit.limits(phi=c(0, 10), lambda=c(-2.5, 2.5))}, \cr
  will change the limits for the parameters \eqn{\phi}{phi} and \eqn{\lambda}{lambda}.
  Default values are used if the argument \code{limits} is not provided.  
  
  If the \code{fix.lambda = FALSE} and \code{nospatial = FALSE} the
  Box-Cox parameter for the model without the spatial component is
  obtained numerically, with log-likelihood computed by the function
  \code{boxcox.ns}.

  \bold{Multiple initial values} can be specified providing a \eqn{n
    \time 2}{n x 2} matrix for the argument \code{ini.cov.pars} and/or
  providing a vector for the values of the remaining model parameters.
  In this case the log-likelihood is computed for all combinations of
  model parameters. The set with results in the maximum value of the
  log-likelihood is then used to start the minimisation algorithm.
  
}


\value{
  An object of the classes \code{"likGRF"} and \code{"variomodel"}.\cr
  The function \code{\link{summary.likGRF}} is used to print a summary
  of the fitted model.\cr 
  The object is a list with the following components:
  \item{cov.model}{a string with the name of the correlation
    function.  }
  \item{nugget}{value of the nugget parameter \eqn{\tau^2}{tau^2}.
    This is an estimate if \code{fix.nugget = FALSE} otherwise, a fixed
    value.  }
  \item{cov.pars}{a vector with the estimates of the parameters
    \eqn{\sigma^2}{sigma^2} 
    and \eqn{\phi}{phi}, respectively.  }
  \item{kappa}{value of the smoothness parameter. Valid only if 
    the correlation function is one of: \code{"matern"},
    \code{"powered.exponential"}, \code{"cauchy"}
    or \code{"gneiting.matern"}.  }
  \item{beta}{estimate of mean parameter \eqn{\beta}{beta}.
    This can be a scalar or vector depending on the
    trend (covariates) specified in the model.  }
  \item{beta.var}{estimated variance (or covariance matrix) for the mean
    parameter \eqn{\beta}{beta}.  }
  \item{lambda}{values of the Box-Cox transformation parameter. A fixed value if
    \code{fix.lambda = TRUE} otherwise the estimate value.  }
  \item{aniso.pars}{fixed values or estimates of the anisotropy
    parameters, according to the function call.  }
  \item{method.lik}{estimation method used, \code{"ML"} (maximum likelihood)
    or \code{"REML"} (restricted maximum likelihood).  }
  \item{loglik}{the value of the maximized likelihood.  }
  \item{npars}{number of estimated parameters.  }
  \item{AIC}{value of the Akaike information criteria.  }
  \item{BIC}{value of the Bayesian information criteria.  }
  \item{parameters.summary}{a data-frame with all model parameters,
    their status (estimated or fixed) and values.  }
  \item{info.minimisation}{results returned by the minimisation function.}
  \item{max.dist}{maximum distance between 2 data points. This
    information relevant for other functions which use outputs from
    \code{likfit}.  }
  \item{trend.matrix}{the trend (covariates) matrix \eqn{X}.  }
  \item{log.jacobian}{numerical value of the logarithm of the Jacobian of
    the transformation.  }
  \item{nospatial}{estimates for the model without the spatial
    component. }
  \item{call}{the function call.  }
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }
%
% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{summary.likGRF}} for summary of the results,
  \code{\link{plot.variogram}}, \code{\link{lines.variogram}} and
  \code{\link{lines.variomodel}} for graphical output,
  \code{\link{proflik}} for computing profile likelihoods,
  \code{\link{variofit}} and for other estimation methods, 
  and \code{\link{optim}} for the numerical minimization function.
}

\examples{
if(is.R()) data(s100)
ml <- likfit(s100, ini=c(0.5, 0.5), fix.nug = TRUE)
ml
summary(ml)
reml <- likfit(s100, ini=c(0.5, 0.5), fix.nug = TRUE, met = "REML")
summary(reml)
plot(variog(s100))
lines(ml)
lines(reml, lty = 2)

\testonly{
if(is.R()) data(s100)
ml <- likfit(s100, ini=c(0.5, 0.5), nug=0.2)
ml <- likfit(s100, ini=c(0.5, 0.5), fix.nug=TRUE, nug=0.2)
ml <- likfit(s100, data=exp(s100$data), ini=c(0.5, 0.5), nug=0.2, fix.lambda=F)
ml <- likfit(s100, ini=c(0.5, 0.5), nug=0.2, fix.psiR = F)
ml <- likfit(s100, ini=c(0.5, 0.5), nug=0.2, fix.psiA = F, fix.psiR = F)
ml <- likfit(s100, ini=c(0.5, 0.5), cov.model="matern", fix.kappa=TRUE, kappa=2)
ml <- likfit(s100, ini=c(0.5, 0.5), cov.model="matern", fix.kappa=F)

## Multiple realizations
ap1 <- grf(40, cov.pars=c(1, .3))
ap2 <- grf(30, cov.pars=c(1, .3))
ap3 <- grf(50, cov.pars=c(1, .3))
ap <- list(coords = rbind(ap1$coords, ap2$coords, ap3$coords), data =
c(ap1$data, ap2$data, ap3$data))
ap$realisations <- c(rep(1,40), rep(2,30), rep(3,50))
ap.fit <- likfit(ap, ini=c(.5, .5), reali=ap$real)

}

}

\keyword{spatial}
\keyword{models}%-- one or more ...
