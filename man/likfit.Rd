\name{likfit}

\alias{likfit}
\alias{logLik.likGRF}
\alias{resid.likGRF}
\alias{residuals.likGRF}
\alias{fitted.likGRF}
\alias{likfit}
\alias{likfit.limits}
\alias{.negloglik.GRF}
\alias{logLik.likGRF}
\alias{check.parameters.values}

\concept{variogram parameter estimation}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Likelihood Based Parameter Estimation for Gaussian Random Fields}

\description{
  \emph{Maximum likelihood} (ML) or \emph{restricted maximum likelihood} (REML)
  parameter estimation for (transformed) Gaussian random fields.
}

\usage{
likfit(geodata, coords = geodata$coords, data = geodata$data,
       trend = "cte", ini.cov.pars, fix.nugget = FALSE, nugget = 0,
       fix.kappa = TRUE, kappa = 0.5, fix.lambda = TRUE, lambda = 1,
       fix.psiA = TRUE, psiA = 0, fix.psiR = TRUE, psiR = 1, 
       cov.model, realisations, lik.method = "ML", components = TRUE,
       nospatial = TRUE, limits = pars.limits(),
       print.pars = FALSE, messages, \dots)

\method{fitted}{likGRF}(object, spatial = TRUE, \ldots)

\method{resid}{likGRF}(object, spatial = FALSE, \ldots)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"}.
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
    Regarded as a fixed value if \code{fix.nugget = TRUE} otherwise
    as the initial value for the minimisation algorithm. Defaults to zero.  }
  \item{fix.kappa}{logical, indicating whether the extra parameter
    \eqn{\kappa}{kappa} should be regarded as fixed
    (\code{fix.kappa = TRUE}) or should be 
    estimated (\code{fix.kappa = FALSE}). Defaults to \code{TRUE}.  }    
  \item{kappa}{value of the extra parameter \eqn{\kappa}{kappa}.
    Regarded as a fixed value if \code{fix.kappa = TRUE}
    otherwise as the initial value for the minimisation algorithm. Defaults to
    \eqn{0.5}. This parameter is valid only if the covariance function is one
    of: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} or
    \code{"gneiting.matern"}.  For more details on covariance functions
    see documentation for \code{\link{cov.spatial}}.  }
  \item{fix.lambda}{logical, indicating whether the Box-Cox transformation parameter
    \eqn{\lambda}{lambda} should be regarded as fixed
    (\code{fix.lambda = TRUE}) or should be
    be estimated (\code{fix.lambda = FALSE}). Defaults to \code{TRUE}.  }    
  \item{lambda}{value of the Box-Cox transformation parameter
    \eqn{\lambda}{lambda}.
    Regarded as a fixed value if \code{fix.lambda = TRUE} otherwise 
    as the initial value for the minimisation algorithm. Defaults to
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
    otherwise  as the initial value for the minimisation algorithm.
    Defaults to \eqn{0}. See \code{\link{coords.aniso}} for further
    details on anisotropy correction.  }
  \item{fix.psiR}{logical, indicating whether the anisotropy ratio parameter
    \eqn{\psi_R}{psi_R} should be regarded as fixed
    (\code{fix.psiR = TRUE}) or should be estimated (\code{fix.psiR = FALSE}). Defaults to
    \code{TRUE}.  }    
  \item{psiR}{value, always greater than 1, for the anisotropy ratio  parameter
    \eqn{\psi_R}{psi_R}. 
    Regarded as a fixed value if \code{fix.psiR = TRUE}
    otherwise as the initial value for the minimisation algorithm.
    Defaults to \eqn{1}. See \code{\link{coords.aniso}} for further
    details on anisotropy correction.  }
  \item{cov.model}{a string specifying the model for the correlation
    function. For further details see documentation for \code{\link{cov.spatial}}.
    Reads values from an \code{variomodel} object passed to \code{ini.cov.pars} if any, otherwise
    defaults to the \emph{exponential} model.  }
  \item{realisations}{optional. Logical or a vector indicating the number of replication
    for each datum. For further information see \code{DETAILS} below and
    documentation for \code{\link{as.geodata}}. }
  \item{lik.method}{(formely method.lik) options are \code{"ML"} for maximum likelihood and \code{"REML"} for
    restricted maximum likelihood. Defaults to \code{"ML"}.  }
  \item{components}{an \eqn{n \times 3}{n x 3} data-frame with fitted
    values for the three model components: trend, spatial and residuals.
    See the section \code{DETAILS} below for the model specification.  }
  \item{nospatial}{logical. If \code{TRUE} parameter estimates for the
    model without spatial component are included in the output.  }
  \item{limits}{values defining lower and upper limits for the model
    parameters used in the numerical minimisation.
    The auxiliary function \code{\link{pars.limits}} is called to set
    the limits.
    See also \bold{Limits} in DETAILS below.}
  \item{print.pars}{logical. If \code{TRUE} the parameters and the value
  of the negative log-likelihood (up to a constant) are printed each
  time the function to be minimised is called.  }
  \item{messages}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{additional parameters to be passed to the minimisation
    function. Typically  arguments of the type  \code{control()} which controls the
    behavior of the minimisation algorithm. For further details see documentation 
    for the minimisation function \code{\link{optim}}.  }
  \item{object}{an object with output of the function \code{likfit}.  }
  \item{spatial}{logical, determines whether the spatial component of
    the model in included in the output. The geostatistical model
    components are: \emph{trend}, \emph{spatial} and \emph{residulas}.
    See DETAILS.}
}

\details{
  This function estimate the parameters of the Gaussian random field model,
  specified as:
  \deqn{Y(x) = \mu(x) + S(x) + e}{Y(x) = mu(x) + S(x) + e}
  where
  \itemize{
    \item \eqn{x} defines a spatial location. Typically Euclidean
    coordinates on a plane. 
    \item \eqn{Y} is the variable been observed.
    \item \eqn{\mu(x) = X\beta}{mu(x) = X \%*\% beta} is the mean component of
    the model (trend). 
    \item \eqn{S(x)} is a stationary Gaussian process with variance \eqn{\sigma^2}{sigma^2}
    (partial sill) and a correlation function parametrized in its
    simplest form by \eqn{\phi}{phi} (the range parameter). Possible extra parameters
    for the correlation function are the smoothness parameter
    \eqn{\kappa}{kappa} and the anisotropy parameters \eqn{\phi_R}{phi_R} and
    \eqn{\phi_A}{phi_A} (anisotropy ratio and angle, respectively).  
    \item \eqn{e} is the error term with variance parameter  \eqn{\tau^2}{tau^2}
    (nugget variance).
  }
  
  The additional parameter \eqn{\lambda}{lambda} allows for the Box-Cox
  transformation of the response variable.
  If used (i.e. if \eqn{\lambda \neq 1}) \eqn{Y(x)} above is replaced by \eqn{g(Y(x))}
  such that
  \deqn{g(Y(x)) = \frac{Y^\lambda(x) - 1}{\lambda}.}{g(Y(x)) = ((Y^lambda(x)) -
    1)/lambda .}

  Two particular cases are \eqn{\lambda = 1}{lambda = 1}
  which indicates no transformation and \eqn{\lambda = 0}{lambda = 0}
  indicating the log-transformation.
  
  \bold{Numerical minimization}

  In general parameter estimation is performed numerically using the \R
  function \code{\link{optim}} to minimise the
  negative log-likelihood computed by the function \code{negloglik.GRF}.
  If the nugget, anisotropy (\eqn{\psi_A, \psi_R}{psiA, psiR}),
  smoothness (\eqn{\kappa}{kappa}) and transformation (\eqn{\lambda}{lambda}) parameters
  are held fixed then the numerical minimisation can be reduced to
  one-dimension and the function \code{\link{optimize}} is used instead
  of \code{optim}. In this case initial values are irrelevant.

  \bold{Limits}
  
  Lower and upper limits for parameter values can be
  individually specified using the function \code{link{pars.limits}}.
  For example, including the following in the function call:\cr
  \code{limits = pars.limits(phi=c(0, 10), lambda=c(-2.5, 2.5))}, \cr
  will change the limits for the parameters \eqn{\phi}{phi} and \eqn{\lambda}{lambda}.
  Default values are used if the argument \code{limits} is not provided.  

  There are internal reparametrisation depending on the options for
  parameters to be estimated.
  For instance for the common situation when \code{fix.nugget=FALSE} the
  minimisation is performed in a reduced
  parameter space using
  \eqn{\tau^2_{rel} = \frac{\tau^2}{\sigma^2}}{tau^2_{rel} = tau^2/sigma^2}.
  In this case values of \eqn{\sigma^2}{sigma^2} and \eqn{\beta}{beta}
  are then given by
  analytical expressions which are function of the two parameters
  remaining parameters and limits for these two parameters will be ignored.

  Since parameter values are found by numerical optimization using 
  the function  \code{\link{optim}}, 
  in given circunstances the algorithm may not converge to correct
  parameter values when called with default options and the user may
  need to pass extra options for the optimizer. For instance the
  function \code{optim} takes a \code{control} argument.
  The user should try different initial values and if the parameters have
  different orders of magnitude may need to use options to scale the parameters.
  Some possible workarounds in case of problems include:
  \itemize{
  \item rescale you data values (dividing by a constant, say)
  \item rescale your coordinates (subtracting values and/or dividing by constants)
  \item Use the  mechanism to pass \code{control()} options for the
  optimiser internally
  }

  \bold{Transformation}  
  If the \code{fix.lambda = FALSE} and \code{nospatial = FALSE} the
  Box-Cox parameter for the model without the spatial component is
  obtained numerically, with log-likelihood computed by the function
  \code{boxcox.ns}.

  \bold{Multiple initial values} can be specified providing a \eqn{n
    \time 2}{n x 2} matrix for the argument \code{ini.cov.pars} and/or
  providing a vector for the values of the remaining model parameters.
  In this case the log-likelihood is computed for all combinations of
  the model parameters. The parameter set which maximises the
  value of the log-likelihood is then used to start the
  minimisation algorithm.

  Alternatively the argument \code{ini.cov.pars} can take an object of
  the class \code{eyefit} or \code{variomodel}. This allows the usage 
  of an output of the functions \code{\link{eyefit}}, \code{\link{variofit}} or
  \code{\link{likfit}} be used as initial value.

  The argument \bold{realisations} allows sets of data \emph{assumed to be
  independent} replications of the same process.
  Data on different realisations may or may not be co-located.
  For instance, data collected at different times
  can be pooled together in the parameter estimation assuming
  time independence.
  The argument \code{realisations} takes a vector indicating the
  replication number (e.g. the times).
  If \code{realisations = TRUE} the code looks for an element
  named \code{realisations} in the \code{geodata} object.
  The log-likelihoods are computed for each replication and added together.
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
  \item{AIC}{value of the Akaike Information Criteria, \eqn{AIC=-2
      ln(L) + 2 p} where \eqn{L} is the maximised likelihood and p is
    the number of parameters in the model.  }
  \item{BIC}{value of the Bayesian Information Criteria,
    \eqn{BIC=-2ln(L) + p log(n)}, where \eqn{n} is the number of data,
      \eqn{L,p} as for AIC above.  }
  \item{parameters.summary}{a data-frame with all model parameters,
    their status (estimated or fixed) and values.  }
  \item{info.minimisation}{results returned by the minimisation function.}
  \item{max.dist}{maximum distance between 2 data points. This
    information relevant for other functions which use outputs from
    \code{likfit}.  }
  \item{trend}{the trend (covariates) matrix \eqn{X}.  }
  \item{log.jacobian}{numerical value of the logarithm of the Jacobian of
    the transformation.  }
  \item{nospatial}{estimates for the model without the spatial
    component. }
  \item{call}{the function call.  }
}

\references{
  Further information on the package
  \pkg{geoR}
  can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
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
  and \code{\link{optim}} for the numerical minimisation function.
}

\examples{
\dontrun{
ml <- likfit(s100, ini=c(0.5, 0.5), fix.nug = TRUE)
ml
summary(ml)
reml <- likfit(s100, ini=c(0.5, 0.5), fix.nug = TRUE, lik.met = "REML")
summary(reml)
plot(variog(s100))
lines(ml)
lines(reml, lty = 2)
}
}
\keyword{spatial}
\keyword{models}%-- one or more ...
