\name{likfit.old}

\alias{likfit.old}
\alias{proflik.phi}
\alias{proflik.ftau}
\alias{proflik.nug}
\alias{loglik.spatial}
\alias{likfit.nospatial}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Likelihood Estimation for Gaussian Random Fields}

\description{
  This function estimates the parameters of a (transformed) Gaussian
  random field using either \emph{maximum likelihood} (ML) or \emph{restricted maximum
    likelihood} (REML).
  Includes option for the estimation the Box-Cox transformation parameter.
  
  This is a earlier version of the function \code{\link{likfit}} and it
  will be made obsolete.
}

\usage{
likfit.old(geodata, coords=geodata$coords, data=geodata$data,
           trend = "cte", ini, fix.nugget = FALSE, nugget = 0, 
           cov.model = "matern", kappa = 0.5, fix.lambda = TRUE,
           lambda = 1, method = "ML", 
           predicted = FALSE, residuals = FALSE, 
           minimisation.function = c("optim","nlmP", "nlm"),
           automatic.refit = FALSE, range.limits,
           messages.screen = TRUE, ...) 
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing the elements \code{coords} and
    \code{data} as described next. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing in each row Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.}
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.}
  \item{trend}{specifies the mean part of the model. The options are:
    \code{"cte"} (constant mean - the default option), \code{"1st"} (a first degree polynomial
    on the coordinates), \code{"2nd"} (a second degree polynomial
    on the coordinates), or a formula of the type \code{~X} where \code{X}
    is a matrix with the covariates (external trend). }
  \item{ini}{initial values of the covariance parameters. These values
    are used to start the numerical minimizer. If the model to be
    fitted has the following 3 parameters: nugget (\eqn{\tau^2}{tau^2}), sill (\eqn{\sigma^2}{sigma^2})
    and range (\eqn{\phi}{phi}), \code{ini} should be a 3 elements vector with initial values for
    these three parameters, in this order.
    If the model to be fitted does not includes the 
    nugget parameter, \code{ini} is a vector with 2 components,
    the initial values for
    \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.
    If a matrix with several
    initial values on the rows is provided,
    a search for the best initial value is
    performed first, taking the one which maximizes the likelihood.  }
  \item{fix.nugget}{logical, indicating whether or not the nugget parameter
    must be included in the estimation. If \code{fix.nugget = T} the
    argument \code{ini} must be a vector of length 2 or a matrix with 2 columns}
  \item{nugget}{value for the nugget parameter when \code{fix.nugget =
      TRUE}. Defaults to zero.}
  \item{cov.model}{correlation function model. See
    documentation for the function \code{\link{cov.spatial}}. 
    Defaults are equivalent to the \emph{exponential} model.  }
  \item{kappa}{extra smoothness parameter required if one of the following correlation
    function is used: \code{"matern"}, \code{"powered.exponential"},
      \code{"cauchy"}
    or \code{"gneiting.matern"}). For more details see documentation for the
    function \code{\link{cov.spatial}}. In this implementation this
    parameter
    is always regarded as fixed during the ML/REML parameter estimation.}
  \item{fix.lambda}{If \code{TRUE}, the default, the transformation
    parameter is regarded as fixed (known) during the estimation
    process, otherwise the ML/REML for this parameter is also computed.
  }
  \item{lambda}{the value of the Box-Cox transformation parameter.
    Two particular cases are \eqn{\lambda = 1}{lambda = 1} which corresponds to no
    transformation and \eqn{\lambda = 0}{lambda = 0} corresponding to
    the log-transformation.
    The value of the transformation parameter \eqn{\lambda}{lambda} is
    regarded as a fixed if \code{fix.lambda = TRUE} otherwise, it is
    used as a initial value.  }
  \item{method}{\code{"ML"} defines maximum likelihood and \code{"REML"}
    defines restricted maximum likelihood. Defaults to \code{"ML"}. }
  \item{predicted}{a matrix with the values for the  components of the predicted
    values  at the data locations as defined by the model adopted.  }
  \item{residuals}{a matrix with values for the components of the random part
    of the model  at the data locations.  }
  \item{automatic.refit}{Defaults to \code{FALSE}. If set to \code{TRUE}
    the model is automatically changed and re-fitted in the following cases:
    \itemize{
      \item if the model includes the nugget and if the estimated relative nugget
      is less than 0.01, the model without nugget is fitted.
      \item if the estimated range parameter is less than \eqn{1e-12} or a value
      provided in the argument \code{range.limits}, a model without
      spatial correlation is fitted.
    }
  }
  \item{minimisation.func}{minimization function to be used. Options are
  \code{"optim"}, the default, \code{"nlmP"} and \code{"nlm"}.  }
  \item{range.limits}{minimum and maximum values allowed for the
    parameter \eqn{\phi}{phi}, the range parameter. Defaults to  \code{c(0, +Inf)}.}
  \item{messages.screen}{logical.  Indicates
    whether or not status messages are printed on the screen (or output device)
    while the function is running.}
  \item{\dots}{additional parameters to be passed to the minimization
    function. Typically \code{control} type arguments which controls the
    behavior of the minimization algorithm. For further details see
    documentation for the minimization functions.  }
}

\details{
  The random field model considered here has the following parameters:
  \itemize{
    \item the mean parameter \eqn{\beta}{beta},
    \item the covariance parameters: partial sill \eqn{\sigma^2}{\sigma^2} and
    range parameter \eqn{\phi}{phi},
    \item the smoothness parameter \eqn{\kappa}{kappa},
    \item the nugget parameter \eqn{\tau^2}{tau^2},
    \item the Box-Cox transformation parameter \eqn{\lambda}{lambda},
    \item and the anisotropy parameters: anisotropy angle
    \eqn{\psi_A}{psi_A} and anisotropy ratio \eqn{\psi_R}{psi_R},
  }
  
  However not all of them can be estimated using this function.
  For this current implementation:
  \enumerate{
    \item the parameters \eqn{\beta}{beta}, \eqn{\sigma^2}{\sigma^2},
    \eqn{\phi}{phi} are always estimated;
    \item the parameters \eqn{\tau^2}{tau^2}, \eqn{\lambda}{lambda} can
    be estimated or considered fixed;
    \item the parameters \eqn{\kappa}{kappa}, \eqn{\alpha}{alpha},
    \eqn{\psi}{psi} are always regarded as fixed values.
  }
  
}

\value{
  An object of the class \code{variomodel} which is a list with the following components:
  \item{cov.model}{a string with the name of the correlation
    function.  }
  \item{nugget}{value of the nugget parameter \eqn{\tau^2}{tau^2}. This
    is the estimated value if \code{fix.nugget = FALSE} or the fixed
    (known) value if \code{fix.nugget = TRUE}.  }
  \item{cov.pars}{a 2 elements vector with the estimated values of the
    partial sill and range
    parameter, respectively.  }
  \item{kappa}{the fixed value of the smoothness parameter required by
    some of the correlation functions.  }
  \item{beta}{estimated mean parameters. Can be a scalar or
    vector depending on the specification of the trend part of the model
    (covariates).  }
  \item{beta.var}{estimated variance (or covariance matrix) of the mean parameter. }
  \item{lambda}{Box-Cox transformation parameter. A fixed (know) value if
    \code{fix.lambda = TRUE} or the estimated value if \code{fix.lambda = FALSE}.}
  \item{loglik}{the value of the maximized likelihood.  }
  \item{npars}{number of estimated parameters.  }
  \item{AIC}{value of the Akaike information criteria.  }
  \item{BIC}{value of the Bayesian information criteria.  }
  \item{trend.ols}{trend (mean parameters) estimated by ordinary least
    squares (i.e. ignoring the spatial correlation).  }
  \item{info.lambda}{stores information about the Box-Cox transformation
    parameter. Typically this is used by other functions which process the
    output of \code{likfit.old()}.  }
  \item{method}{estimation method used, \code{"ML"} (maximum likelihood)
    or \code{"REML"} (restricted maximum likelihood).}
  \item{s2}{estimated residual variance. Only returned if
    \code{predicted = TRUE} or \cr \code{residuals = TRUE}.}
  \item{predicted}{predicted values at data locations. Only returned if
    \code{predicted = TRUE}.}
  \item{residuals}{estimated residuals at data locations. Only returned if
    \code{residuals = TRUE}.}
  \item{info}{results returned by the minimization function.}
  \item{max.dist}{maximum distance between two data locations. Typically
    this value is used by other functions which process results from
    \code{likfit.old()}.  }
  \item{trend.matrix}{the trend (covariates) matrix.  }
  \item{call}{the function call.  }
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\note{\code{proflik.phi}, \code{proflik.ftau}, \code{proflik.nug} are
  auxiliary functions which computed the value of the likelihood for
  models without nugget, with a fixed nugget and including nugget,
  respectively.
  These functions are internally called by the minimization function
  during the numerical 
  estimation of the model parameters.
}

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{summary.variomodel}},
  \code{\link{plot.variogram}},
  \code{\link{lines.variogram}},
  and \code{\link{lines.variomodel}} for graphical output.
  Profile likelihoods can be computed by \code{\link{proflik}}.
  \code{\link{olsfit}} and 
  \code{\link{wlsfit}} performs variogram based estimation while 
  \code{\link{krige.bayes}} performs Bayesian inference.
  For details about the minimization functions see the documentation for
  \code{\link{optim}} and \code{\link{nlm}}.  
}

\examples{
if(is.R()) data(s100)
ml <- likfit(s100, ini=c(.5, .5), fix.nug=TRUE)
summary(ml)
reml <- likfit(s100, ini=c(.5, .5), fix.nug=TRUE, met="REML")
summary(reml)
plot(variog(s100))
lines(ml)
lines(reml, lty=2)


}

\keyword{spatial}
