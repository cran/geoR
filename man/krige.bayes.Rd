\name{krige.bayes}

\alias{krige.bayes}
\alias{krige.bayes.aux1}
\alias{krige.bayes.aux2}
\alias{krige.bayes.aux3}
\alias{krige.bayes.aux10}
\alias{krige.bayes.aux20}
\alias{krige.bayes.messages}
\alias{model.control}
\alias{prior.control}
\alias{output.control}

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
            locations = "no",
     model = model.control(trend.d = "cte", trend.l = "cte",
                   cov.model = "exponential", kappa = 0.5,
                   aniso.pars = NULL, lambda = 1),
     prior = prior.control(beta.prior = c("flat", "normal",
                                          "fixed"),
                   beta = NULL, beta.var = NULL,
                   sill.prior = c("reciprocal", "fixed"),
                   sill = NULL,
	           range.prior = c("uniform", "exponential",
                                   "fixed", "squared.reciprocal",
                                   "reciprocal"),
                   exponential.prior.par = 1, range = NULL,
                   range.rel.discrete = NULL,
                   nugget.rel.prior = c("fixed", "uniform"),
                   nugget.fixed = 0, nugget.discrete = NULL),
     output = output.control(n.posterior = 1000,
                   n.predictive = NULL, moments = TRUE,
                   simulations.predictive = TRUE,
	           keep.simulations = TRUE, mean.estimator = TRUE,
                   quantile.estimator =  NULL,
                   probability.estimator = NULL,
                   signal = FALSE, messages.screen = TRUE))
}

%- maybe also `usage' for other objects documented here.

\section{Auxiliary functions}{
  The functions of type \code{krige.bayes.aux} and \code{control} are auxiliary functions called by \code{krige.bayes}. 
}

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
  \item{model}{defines model components. See section DETAILS below.  }
  \item{prior}{specification of priors for the model parameters. See
    section DETAILS below.  }
  \item{output}{Defines output options. See section DETAILS below.  }
}

\details{
  \code{krige.bayes} is a generic function for Bayesian geostatistical 
  analysis where predictions can take into account the  parameter
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
  The arguments for the control functions are as follows.\cr
  
  \bold{ARGUMENTS FOR CONTROL FUNCTIONS}
  
  \emph{Arguments for \code{model = model.control(\dots)} :}

  \describe{
    
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
      
    }


    \emph{Arguments for  \code{prior = prior.control(\dots)} :}
    
    \describe{

      \item{beta.prior}{prior distribution for the mean (vector)
        parameter \eqn{\beta}{beta}. The options are "flat" (default),
	"normal" or "fixed" (known mean).  }
      
      \item{beta}{hyperparameter for the distribution of the mean (vector) parameter        \eqn{\beta}{beta}. Only used if \code{beta.prior = "normal"} or
        \code{beta.prior = "fixed"}. For the later \code{beta} defines the value of
        the known mean.  }
      
      \item{beta.var}{(co)variance hyperparameters for the prior of the mean
        (vector) parameter \eqn{\beta}{beta}. Only used if \code{beta.prior =
          "normal"}.  }
      
      \item{sill.prior}{specifies the prior for the parameter
	\eqn{\sigma^2}{sigma^2}. If \code{"reciprocal"} (the default), the prior 
        \eqn{\frac{1}{\sigma^2}}{1/sigma^2} is used. Otherwise the
	parameter is regarded as fixed.  }
      
      \item{sill}{fixed value of the sill parameter \eqn{\sigma^2}{sigma^2}. Only used if
        \code{sill.prior = FALSE}.  }
      
      \item{range.prior}{prior distribution for the range parameter
        \eqn{\phi}{phi}. Options are: \code{\"uniform\"}, (\eqn{\propto
          1}),
        \code{\"exponential\"} ,(\eqn{\exp(-\nu * \phi)}{exp(- nu * phi)})
        \code{\"reciprocal\"} (\eqn{\frac{1}{\phi}}{1/phi}),
        \code{\"squared.reciprocal\"} (\eqn{\frac{1}{\phi^2}}{1/phi^2}),               \code{\"fixed\"} (known/estimated value of  \eqn{\phi}{phi}).
        If the latter is defined the parameter \eqn{\phi}{phi} is
        fixed  when performing prediction.  }
      
      \item{exponential.prior.par}{parameter \eqn{\nu}{nu} for the exponential
        prior for the parameter \eqn{\phi}{phi}. Only used if
	\code{prior.phi = "exponential"}.  }
      
      \item{range}{fixed value of the range parameter \eqn{\phi}{phi}. Only needed if
        \code{range.prior = "fixed"}.  }
      
      \item{range.discrete}{support points for the discretisation of the prior
        for the range parameter \eqn{\phi}{phi}.  }
      
      \item{nugget.rel.prior}{specifies a prior distribution for the
	relative nugget
        parameter \eqn{\frac{\tau^2}{\sigma^2}}{tau^2/sigma^2}. If \code{nugget.rel.prior = "fixed"} the nugget is
        considered known (fixed) with value given by the argument
        \code{nugget.fixed}.
        If \code{nugget.rel.prior = "uniform"} a discrete uniform prior is used
        with support points given by the argument
	\code{nugget.rel.discrete}.  }
      
      \item{nugget.fixed}{fixed value for the nugget parameter. Only used if
        \code{prior.nugget = "fixed"}.  }
      
      \item{nugget.rel.discrete}{support points for the discretisation of the prior
        for the relative nugget parameter
        \eqn{\frac{\tau^2}{\sigma^2}}{tau^2/sigma^2}.  }
      
    }
    
    
    \emph{Arguments for  \code{output = output.control(\dots)} :}
    
    \describe{
      
      \item{n.posterior}{number of samples to be taken from the posterior
	distribution.  }
      
      \item{n.predictive}{number of samples to be taken from the
	predictive distribution. By default equals to
	\code{n.posterior}.  }
      
      \item{moments}{logical. If \code{TRUE} moments of the predictive distribution
	are computed analytically (without sampling).
	Valid only if \code{lambda = 1}
	or \code{lambda = 0}.  }
      
      \item{simulations.predictive}{logical. Defines whether simulations are
	drawn from the predictive distribution. Only valid if prediction
	locations are provided on the argument \code{locations}.  } 
      
      \item{keep.simulations}{logical. Indicates
	whether or not the samples of the predictive distributions are
	returned.
	Only valid if prediction locations are provided on the argument
	\code{locations}.  }
      
      \item{mean.estimator}{logical. Indicates
	whether or not the mean and variances of the predictive distributions are
	computed and returned. If \code{TRUE} the objects \code{predict.mean} and
	\code{krige.var} are included in the output.
	Only valid if prediction locations are provided on the argument
	\code{locations}.  } 
      
      \item{quantile.estimator}{indicates
	whether or not quantiles of the predictive distributions are
	computed and returned.
	If a vector with numbers in the interval
	\eqn{[0,1]}{[0,1]} is provided the output includes the object
	\code{quantiles}, which contains values of corresponding
	estimated quantiles. For
	example, if \code{estimator = c(0.25,0.50,0.75)} the function
	returns the quartiles of the distributions at each of the
	prediction locations.
	If \code{quantile.estimator = TRUE}, the default, the vector \code{c(0.025, 0.5,
	  0.975)}, is assumed. 
	A measure of uncertainty for the predictions,  which is analogous to the kriging
	standard error, can be computed by \eqn{(quantile0.975 - quantile0.025)/4}.
	Only used if prediction locations are provided in the argument
	\code{locations}.  }
      
      \item{probability.estimator}{the default is \code{FALSE} for which
	case
	nothing is computed. If some
	cutoff values are provided instead, an object called
	\code{probability} is included in the output.  This object contains, for each
	prediction location, the probability that the variable is less
	than or equal to the
	cutoff value given in the argument.  }
      
      \item{signal}{logical. If \code{TRUE} the signal is predicted, otherwise the
	variable is predicted. If no transformation is performed the
	expectations are the same in both cases and
	the kriging variances are different, if the nugget is different of zero.  }
      
      \item{messages.screen}{logical. Indicates
	whether or not status messages are printed on the screen (or output other device)
	while the function is running.  }
      
    }

}

    

  
\value{An object of the \code{\link{class}} \code{"krige.bayes"} which
  is a list with the following components:
  
  \item{posterior}{A list with results for the posterior distribution of the
    model parameters. The components are:\cr  }
  
  \itemize{   
    \item{beta.summary}{summary for the posterior distribution
      of the mean parameter \eqn{\beta}{beta}.  }
    \item{sigmasq.summary}{summary for the posterior distribution
      of the variance parameter \eqn{\sigma^2}{sigma^2} (partial sill).  }
    \item{phi.summary}{summary for the posterior distribution
      of the correlation parameter \eqn{\phi}{phi} (range parameter) .  }
    \item{tausq.summary}{summary for the posterior distribution
      of the nugget variance parameter \eqn{\tau^2}{tau^2}.  }
    \item{beta.samples}{samples from the posterior distribution of the
      mean parameter \eqn{\beta}{beta}.  }
    \item{sigmasq.samples}{samples from the posterior distribution
      of the variance parameter \eqn{\sigma^2}{sigma^2}.  }
    \item{phi.samples}{samples from the posterior distribution
      of the correlation parameter \eqn{\phi}{phi}.  }
    \item{tausq.samples}{samples from the posterior distribution
      of the nugget variance parameter \eqn{\tau^2}{tau^2}.  }
    \item{phi.marginal}{samples from the marginal posterior distribution
      of the correlation parameter \eqn{\phi}{phi}, resulting from
      averaging the posterior  over the
      distribution of \eqn{(\beta, \sigma^2)}{(beta, sigma^2)}.  }
    \item{nugget.marginal}{samples from the marginal posterior
      distribution of the nugget variance parameter \eqn{\tau^2}{tau^2},
      resulting from averaging the posterior over the
      distribution of \eqn{(\beta, \sigma^2)}{(beta, sigma^2)}.  }
    \cr
    \cr
  }

  \cr
  \cr
  
  \item{predictive}{A list with results for the predictive distribution of the
    prediction locations (if provided). The
    components are:
    \cr
    \cr
  }
  
  \itemize{    
    \item{moments}{a numerical matrix. The columns contains the
      estimated
      moments of the
      predictive distribution, at each prediction location}
    
    \item{simulations}{a numerical matrix. Each column has a simulation
      from the predictive distribution. Returned
      only if \code{keep.simulations = TRUE}.  }
    \item{mean.simulations}{a vector with the estimated mean at the prediction
      locations computed by averaging over the simulations. Returned
      only if \code{mean.estimator = TRUE}.  }
    \item{variance.simulations}{a vector with the estimated variance at the
      prediction locations, computed using the simulations. Returned
      only if \code{mean.estimator = TRUE}. }
    \item{quantile}{A matrix or vector with quantile estimators. Returned
      only if the argument \code{quantile.estimator} is used.  }
    \item{probability}{A matrix or vector with probability estimators. Returned
      only if the argument \code{probability.estimator} is used.  }
    \item{type.prediction}{information on the type of prediction
      performed.  }
    \cr
    \cr
  }

  \cr
  \cr
  
  \item{message.prediction}{information about the parameter uncertainty
    taken into account. Indicates which parameters has been regarded as
    random during the analysis. }
  
  \item{.Random.seed}{system random seed before running the function.
    Allows reproduction of results. If
    the \code{\link{.Random.seed}} is set to this value and the function is run
    again, it will produce exactly the same results.  }
  
  
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
ex.grid <- as.matrix(expand.grid(seq(0,1,l=11), seq(0,1,l=11)))
#
# computing Bayesian posterior and predictive distributions
# (this can take some time to run)
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, prior =
                 prior.control(range.discrete=seq(0, 2, l=51)))
#
# Ploting theoretical amd empirical variograms
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
image.krige.bayes(ex.bayes, loc=ex.grid, main="predicted values")
image.krige.bayes(ex.bayes, val="moments.variance",
                  loc=ex.grid, main="prediction variance")
image.krige.bayes(ex.bayes, val= "simulation", number.col=1,
                  loc=ex.grid,
      main="a simulation from the \npredictive distribution")
image.krige.bayes(ex.bayes, val= "simulation", number.col=2,
                  loc=ex.grid,
     main="another simulation from \nthe predictive distribution")
par(op)
}
\testonly{
ex.data <- grf(50, cov.pars=c(10, .25))
ex.grid <- as.matrix(expand.grid(seq(0,1,l=11), seq(0,1,l=11)))
ex.bayes <- krige.bayes(ex.data, loc=ex.grid, prior =
                 prior.control(range.discrete=seq(0, 2, l=3),
                 nugget.rel.prior = "uniform",
                 nugget.rel.discrete=seq(0,.5, l=2)),
                 output=output.control(n.post=100))
plot(ex.data)
lines(ex.bayes)
lines(ex.bayes, summ="median", lty=2)
lines(ex.bayes, summ="mean", lwd=2, lty=2)
op <- par(no.readonly = TRUE)
par(mfrow=c(2,2))
par(mar=c(3,3,1,1))
par(mgp = c(2,1,0))
image.krige.bayes(ex.bayes, loc=ex.grid, main="predicted values")
image.krige.bayes(ex.bayes, val="moments.variance",
                  loc=ex.grid, main="prediction variance")
image.krige.bayes(ex.bayes, val= "simulation", number.col=1,
                  loc=ex.grid,
      main="a simulation from the \npredictive distribution")
image.krige.bayes(ex.bayes, val= "simulation", number.col=2,
                  loc=ex.grid,
     main="another simulation from \nthe predictive distribution")
par(op)
}

}

\keyword{spatial}
\keyword{models}

