\name{krige.bayes}

\alias{krige.bayes}
\alias{model.control}
\alias{prior.control}

\title{Bayesian Analysis for Gaussian Geostatistical Models}

\description{
  The function \code{krige.bayes} performs Bayesian analysis of
  geostatistical data allowing specifications of
  different levels of uncertainty in the model parameters.\cr
  It returns results on the posterior distributions for the model
  parameters and on the predictive distributions for prediction
  locations (if provided).
}

\usage{
krige.bayes(geodata, coords = geodata$coords, data = geodata$data,
            locations = "no", borders = NULL, model = model.control(),
            prior = prior.control(), output = output.control())

model.control(trend.d = "cte", trend.l = "cte", cov.model = "matern",
              kappa = 0.5, aniso.pars = NULL, lambda = 1)

prior.control(beta.prior = c("flat", "normal", "fixed"),
              beta = NULL, beta.var.std = NULL,
              sigmasq.prior = c("reciprocal", "uniform",
                                "sc.inv.chisq", "fixed"),
              sigmasq = NULL, df.sigmasq = NULL,
              phi.prior = c("uniform", "exponential","fixed",
                            "squared.reciprocal", "reciprocal"),
              phi = NULL, phi.discrete = NULL,
              tausq.rel.prior = c("fixed", "uniform"),
              tausq.rel = 0, tausq.rel.discrete = NULL)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix where each row has the 2-D
    coordinates of the \eqn{n} data locations.
    By default it takes the
    component \code{coords} of the argument \code{geodata}, if provided.
  }
  \item{data}{a vector with \emph{n} data values. By default it takes the
    component \code{data} of the argument \code{geodata}, if provided.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data-frame with the 2-D
    coordinates of the \eqn{N} prediction locations.
    Defaults to \code{"no"} in which case
    the function returns only results on the posterior distributions of
    the model parameters.  }
  \item{borders}{optional. If a two column matrix defining a polygon is
    provided the prediction is performed only at locations inside this polygon. }
  \item{model}{defines model components. See section DETAILS below.  }
  \item{prior}{specification of priors for the model parameters. See
    section DETAILS below.  }
  \item{output}{Specify output options. See documentation for
    \code{\link{output.control}}.  }
  
  \item{trend.d}{specifies the trend (covariates) values at the data
    locations. Possible values are:
    \code{"cte"} - model with constant mean,
    \code{"1st"} - trend is defined as a first order polynomial on
    the coordinates,
    \code{"2nd"} - trend is defined as a second order polynomial on the
    coordinates,  a formula of the type \code{~X}, where \code{X}
    is a matrix with covariates (external trend) at data locations. 
    Defaults to \code{"cte"}.  }
  
  \item{trend.l}{specifies the trend (covariates) at the prediction
    locations. Must be of the same type as defined for \code{trend.d}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  
  \item{cov.model}{string indicating the name of the model for the correlation function. Further details in the
    documentation for \code{\link{cov.spatial}}. }
  
  \item{kappa}{additional smoothness parameter. Only used if the
    correlation function is one of: \code{"matern"},
    \code{"powered.exponential"}, \code{"cauchy"} or
    \code{"gneiting.matern"}. In the current implementation this
    parameter is always regarded as fixed during the Bayesian analysis.  }
  
  \item{aniso.pars}{fixed parameters for geometric anisotropy
    correction. If \code{aniso.pars = FALSE} no correction is made, otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link{coords.aniso}}.  }
  
  \item{lambda}{numerical value of the Box-Cox transformation parameter.
    The value \eqn{\lambda = 1}{lambda = 1} corresponds to
    no transformation. The Box-Cox parameter \eqn{\lambda}{lambda}
    is always regarded as fixed and
    data transformation is 
    performed before the analysis.
    Prediction results are back-transformed and
    returned is the same scale as for the original data.
    For \eqn{\lambda = 0}{lambda = 0} the log-transformation is performed.
    If \eqn{\lambda < 0}{lambda < 0} the mean predictor doesn't make sense
    (the resulting distribution has no expectation).  }


  
  \item{beta.prior}{prior distribution for the mean (vector)
    parameter \eqn{\beta}{beta}. The options are "flat" (default),
    "normal" or "fixed" (known mean).  }
  
  \item{beta}{mean hyperparameter for the distribution of the mean (vector) parameter        \eqn{\beta}{beta}. Only used if \code{beta.prior = "normal"} or
    \code{beta.prior = "fixed"}. For the later \code{beta} defines the value of
    the known mean.  }
  
  \item{beta.var.std}{standardised (co)variance hyperparameter(s)
    for the prior for the mean
    (vector) parameter \eqn{\beta}{beta}.
    The (co)variance matrix for\eqn{\beta}{beta} is given by the
    multiplication of this matrix by \eqn{\sigma^2}{sigma^2}.
    Only used if \code{beta.prior = "normal"}.  }
  
  \item{sigmasq.prior}{specifies the prior for the parameter
    \eqn{\sigma^2}{sigma^2}. If \code{"reciprocal"} (the default), the prior 
    \eqn{\frac{1}{\sigma^2}}{1/sigma^2} is used. Otherwise the
    parameter is regarded as fixed.  }
  
  \item{sigmasq}{fixed value of the sill parameter
    \eqn{\sigma^2}{sigma^2}. Only used if 
    \code{sigmasq.prior = FALSE}.  }
  \item{df.sigmasq}{numerical. Number of degrees of freedom for the
    prior for the parameter \eqn{\sigma^2}{sigma^2}. Only used if
    \code{sigmasq.prior = "sc.inv.chisq"}.  }
  \item{phi.prior}{prior distribution for the range parameter
    \eqn{\phi}{phi}. Options are:
    \code{\"uniform\"}, (\eqn{p(\phi) \propto 1}{p(phi) propto 1}),
    \code{\"exponential\"} ,(\eqn{p(\phi) = \frac{1}{\nu} \exp(- \frac{1}{\nu} *
      \phi)}{(-1/nu) exp(- (1/nu) * phi)})
    \code{\"reciprocal\"} (\eqn{p(\phi) \propto
      \frac{1}{\phi}}{p(\phi) propto 1/phi}),
    \code{\"squared.reciprocal\"} (\eqn{p(\phi) \propto \frac{1}{\phi^2}}{p(\phi) propto 1/phi^2}),
    \code{\"fixed\"} (known/estimated value of  \eqn{\phi}{phi}).
    If the latter is defined the parameter \eqn{\phi}{phi} is
    fixed  when performing prediction.
    For the exponential prior the argument \code{phi} provides the
    the hyperparameter \eqn{\nu}{nu} and corresponds to the
    expected value for this distribution.  }
  
  
  \item{phi}{fixed value of the range parameter \eqn{\phi}{phi}. Only needed if
    \code{phi.prior = "fixed"}.  }
  
  \item{phi.discrete}{support points of the discrete prior
    for the range parameter \eqn{\phi}{phi}.  }
  
  \item{tausq.rel.prior}{specifies a prior distribution for the
    relative nugget parameter
    \eqn{\frac{\tau^2}{\sigma^2}}{tau^2/sigma^2}.
    If \code{tausq.rel.prior = "fixed"} the relative nugget is
    considered known (fixed) with value given by the argument
    \code{tausq.rel}.
    If \code{tausq.rel.prior = "uniform"} a discrete uniform prior is used
    with support points given by the argument
    \code{tausq.rel.discrete}.  }
  
  \item{tausq.rel}{fixed value for the relative nugget parameter.
    Only used if
    \code{tausq.rel.prior = "fixed"}.  }
  
  \item{tausq.rel.discrete}{support points of the discrete prior
    for the relative nugget parameter
    \eqn{\frac{\tau^2}{\sigma^2}}{tau^2/sigma^2}.  }
  
  
  
}


\details{
  \code{krige.bayes} is a generic function for Bayesian geostatistical 
  analysis of (transformed) Gaussian where predictions take into account the  parameter
  uncertainty.
  
  It can be set to run conventional kriging methods which
  use known parameters or \emph{plug-in} estimates. However, the
  functions \code{krige.conv} and \code{ksline} are preferable for
  prediction with fixed parameters.
  
  The basis for the Bayesian algorithm is to discretize the prior
  distribution for the parameters \eqn{\phi}{phi} and \eqn{\tau^2_{rel}
    =\frac{\tau^2}{\sigma^2}}{tau_rel = tau/sigma}.
  The Tech. Report referenced below provides details on the results used 
  in the current implementation.\cr
  
  \bold{CONTROL FUNCTIONS}
  
  The function call includes auxiliary control functions which allows
  the user to specify and/or change the specification of model
  components
  (using \code{model.control}), prior
  distributions (using \code{prior.control}) and 
  output options (using \code{output.control}).
  Default options are available in most of the cases.
  
  
}


\value{An object with \code{\link{class}} \code{"krige.bayes"} and
  \code{"kriging"}.
  The attribute \code{prediction.locations} containing the name of the
  object with the coordinates of the prediction locations (argument
  \code{locations}) is assigned to the object.
  Returns a list with the following components:
  
  \item{posterior}{results on on the posterior distribution of the
    model parameters. A list with the following possible components:\cr  }
  
  \itemize{   
    \item{beta}{summary information on the posterior distribution
      of the mean parameter \eqn{\beta}{beta}.  }
    \item{sigmasq}{summary information on the posterior distribution
      of the variance parameter \eqn{\sigma^2}{sigma^2} (partial sill).  }
    \item{phi}{summary information on the posterior distribution
      of the correlation parameter \eqn{\phi}{phi} (range parameter) .  }
    \item{tausq.rel}{summary information on the posterior distribution
      of the relative nugget variance parameter
      \eqn{\tau^2_{rel}}{tau^2_rel}.  }
    \item{joint.phi.tausq.rel}{information on discrete the joint
      distribution of these parameters.  }
    \item{sample}{a data.frame with a sample from the posterior
      distribution.
      Each column correspondes to one of the basic model parameters.  }
    \cr
    \cr
  }
  
  \cr
  \cr
  
  \item{predictive}{results on the predictive distribution at the
    prediction locations, if provided.
    A list with the following possible components:\cr  }
  
  \itemize{    
    \item{mean}{expected values.  }
    \item{variance}{expected variance.  }
    \item{distribution}{type of posterior distribution.  }
    \item{mean.simulations}{mean of the simulations at each locations.  }
    \item{variance.simulations}{variance of the simulations at each locations.  }
    \item{quantiles.simulations}{quantiles computed from the the simulations.  }
    \item{probabilities.simulations}{probabilities computed from the simulations.  }
    \item{simulations}{simulations from the predictive distribution.  }
  }
  
  \cr
  \cr
  
  \item{prior}{a list with information on the prior distribution
    and hyper-parameters of the model parameters (\eqn{\beta, \sigma^2,
      \phi, \tau^2_{rel}}{beta, sigma^2, phi, tau^2_rel}).
  }
  
  \item{model}{model specification as defined by \code{model.control}.  }
  
  \item{.Random.seed}{system random seed before running the function.
    Allows reproduction of results. If
    the \code{\link{.Random.seed}} is set to this value and the function is run
    again, it will produce exactly the same results.  }
  
  \item{max.dist}{maximum distance found between two data locations.  }
  \item{call}{the function call.  }
  
}


\references{
  The technical details about the  implementation of \code{krige.bayes} can be
  found at:
  \cr
  
  Ribeiro, P.J. Jr. and Diggle, P.J. (1999) \emph{Bayesian inference in
    Gaussian model-based geostatistics}. Tech. Report ST-99-08, Dept
  Maths and Stats, Lancaster University.\cr
  Available at: \cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/publications.html}
  
  Further information about \bold{geoR} can be found at:
  \cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
  
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}


\seealso{\code{\link{lines.krige.bayes}},
  \code{\link{image.krige.bayes}} and
  \code{\link{persp.krige.bayes}} for graphical output of the results.  
  \code{\link{krige.conv}} and
  \code{\link{ksline}} for conventional kriging methods.
}


\examples{
\dontrun{
# generating a simulated data-set
ex.data <- grf(50, cov.pars=c(10, .25))
#
# defining the prediction grid:
ex.grid <- as.matrix(expand.grid(seq(0,1,l=21), seq(0,1,l=21)))
#
# computing Bayesian posterior and predictive distributions
# (this can take some time to run)
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, prior =
                 prior.control(phi.discrete=seq(0, 2, l=21)))
#
# Prior and posterior for the parameter phi
plot(ex.bayes, type="h", tausq.rel = FALSE, col=c("red", "blue"))
#
# Plot histograms of samples from the posterior
par(mfrow=c(2,1))
hist(ex.bayes)
par(mfrow=c(1,1))
#
# Ploting theoretical and empirical variograms
plot(ex.data)
# adding lines with fitted variograms
lines(ex.bayes)
lines(ex.bayes, summ="median", lty=2)
lines(ex.bayes, summ="mean", lwd=2, lty=2)
#
# Ploting prediction some results
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(mgp = c(2,1,0))
image(ex.bayes, main="predicted values")
image(ex.bayes, val="variance", main="prediction variance")
image(ex.bayes, val= "simulation", number.col=1,
      main="a simulation from the \npredictive distribution")
image(ex.bayes, val= "simulation", number.col=2,
      main="another simulation from \nthe predictive distribution")
par(op)
}
\testonly{
ex.data <- grf(50, cov.pars=c(10, .25))
ex.post <- krige.bayes(ex.data)
ex.grid <- as.matrix(expand.grid(seq(0,1,l=6), seq(0,1,l=6)))
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, prior =
                 prior.control(phi.discrete=seq(0, 2, l=3),
                 tausq.rel.discrete=seq(0, 2, l=3)),
                 output=output.control(n.post=100))
plot(ex.data)
lines(ex.bayes)
plot(ex.bayes)
lines(ex.bayes, summ="median", lty=2)
lines(ex.bayes, summ="mean", lwd=2, lty=2)
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(mgp = c(2,1,0))
image(ex.bayes, main="predicted values")
image(ex.bayes, val="variance", main="prediction variance")
image(ex.bayes, val= "simulation", number.col=1,
      main="a simulation from the \npredictive distribution")
image(ex.bayes, val= "simulation", number.col=2,
      main="another simulation from \nthe predictive distribution")
par(op)
}

}

\keyword{spatial}
\keyword{models}

