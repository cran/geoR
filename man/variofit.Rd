\name{variofit}

\alias{variofit}
\alias{.loss.vario}

\concept{variogram parameter estimation}

\title{ Variogram Based Parameter Estimation }
\description{
  Estimate covariance parameters by fitting a parametric model to a
  empirical variogram. Variograms models can be fitted by using weighted or ordinary least squares.
}
\usage{
variofit(vario, ini.cov.pars, cov.model,
         fix.nugget = FALSE, nugget = 0,
         fix.kappa = TRUE, kappa = 0.5,
         simul.number = NULL, max.dist = vario$max.dist,
         weights, minimisation.function,
         limits = pars.limits(), messages, \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{vario}{an object of the class \code{"variogram"}, typically an output of the function
    \code{\link{variog}}. The object is a list with information about the
    empirical variogram.  }
  \item{ini.cov.pars}{initial values for the covariance parameters:
    \eqn{\sigma^2}{sigma^2} (partial sill) and \eqn{\phi}{phi} (range
    parameter). See \code{DETAILS} below.  }
  \item{cov.model}{a string with the name of the correlation
    function. For further details see documentation for
    \code{\link{cov.spatial}}. For the linear model use \code{cov.model
      = "linear"}. Read values from \code{variomodel} object passed \code{ini.cov.pars},
     otherwise default is the \emph{exponential} model.  }
  \item{fix.nugget}{logical, indicating whether the parameter
    \eqn{\tau^2}{tau^2} (nugget variance) should be regarded as fixed
    (\code{fix.nugget = TRUE}) or should
    be estimated (\code{fix.nugget = FALSE}). Defaults to
    \code{FALSE}.  }    
  \item{nugget}{value for the nugget parameter.  Regarded as a
    fixed values if \code{fix.nugget = TRUE} or as a initial value for the
    minimization algorithm if \code{fix.nugget = FALSE}. 
    Defaults to zero.  }
  \item{fix.kappa}{logical, indicating whether the parameter
    \eqn{\kappa}{kappa} should be regarded as fixed or 
    be estimated. Defaults to \code{TRUE}.  }    
  \item{kappa}{value of the smoothness parameter.  Regarded as a
    fixed values if \code{fix.kappa = TRUE} or as a initial value for the
    minimization algorithm if \code{fix.kappa = FALSE}. Only required if
    one of the following correlation functions is used: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"}
    and \code{"gneiting.matern"}. Defaults to \eqn{0.5}.  }
  \item{simul.number}{number of simulation. To be used when the object passed to the
    argument \code{vario} has empirical variograms for more than one
    data-set (or simulation). Indicates to which one the model will be
    fitted.  }
  \item{max.dist}{maximum distance considered when fitting the
    variogram. Defaults to \code{vario$max.dist}.  }
  \item{weights}{type weights used in the loss function. See
    \code{DETAILS} below.  } 
  \item{limits}{values defining lower and upper limits for the model
    parameters used in the numerical minimisation.
    Only valid if \code{minimisation.function = "optim"}.
    The auxiliary function \code{\link{pars.limits}} is called to set the limits.  }
  \item{minimisation.function}{minimization function used to estimate
    the parameters. Options are \code{"optim"}, \code{"nlm"}.
    If \code{weights = "equal"} the option 
    \code{"nls"} is also valid and det as default.
    Otherwise defaults to \code{"optim"}.  }
  \item{messages}{logical. Indicates
    whether or not status messages are printed on the screen (or other output device)
    while the function is running.  }
  \item{\dots}{further parameters to be passed to the minimization
    function. Typically arguments of the type \code{control()} which controls the
    behavior of the minimization algorithm. See documentation for the 
    selected minimization function for further details.  }
}

\details{

  \bold{Numerical minimization}

  The parameter values are found by numerical optimization using one of
  the functions: \code{\link{optim}}, \code{\link{nlm}} and \code{\link{nls}}.
  In given circunstances the algorithm may not converge to correct
  parameter values when called with default options and the user may
  need to pass extra options for the optimizers. For instance the
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
  
  \bold{Initial values}
  
  The algorithms for minimization functions require initial values of
  the parameters.
  
  A unique initial value is used if a vector is provided in the argument
  \code{ini.cov.pars}. The elements are initial values for
  \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.
  This vector is concatenated with the value of the
  argument \code{nugget} if \code{fix.nugget = FALSE} and \code{kappa}
  if \code{fix.kappa = TRUE}.
  
  Specification of multiple initial values is also possible.
  If this is the case, the function
  searches for the one which minimizes the loss function and uses this as
  the initial value for the minimization algorithm.
  Multiple initial values are specified  by providing a matrix in the
  argument
  \code{ini.cov.pars} and/or, vectors in the arguments 
  \code{nugget} and \code{kappa} (if included in the estimation).
  If \code{ini.cov.pars} is a matrix, the first column has values of
  \eqn{\sigma^2}{sigma^2} and the second has values of  \eqn{\phi}{phi}. 

  Alternatively the argument \code{ini.cov.pars} can take an object of
  the class \code{eyefit} or \code{variomodel}. This allows the usage 
  of an output of the functions \code{\link{eyefit}}, \code{\link{variofit}} or
  \code{\link{likfit}} be used as initial value.
  
  If \code{minimisation.function = "nls"} only the values of
  \eqn{\phi}{phi} and \eqn{\kappa}{kappa} (if this is included in the
  estimation) are used. Values for the remaning are not need by the algorithm.
  
  If \code{cov.model = "linear"} only the value of
  \eqn{\sigma^2}{sigma^2} is used. Values for the
  remaning are not need by this algorithm.
  
  If \code{cov.model = "pure.nugget"} no initial values are needed since
  no minimisation function is used.

  
  \bold{Weights}
  
  The different options for the argument \code{weights}
  are used to define the loss function to be minimised.
  The available options are as follows.
  
  \describe{
    \item{\code{"npairs"}}{indicates that the weights are given by the
      number of pairs in each bin.
      This is the default option unless \code{variog$output.type ==
      "cloud"}.
      The loss function is:
      \deqn{LOSS(\theta) = \sum_k n_k [(\hat{\gamma}_k) -
	\gamma_k(\theta)]^2}{LOSS(theta) = sum_k n_k (hat(gamma) -
	gamma(theta))^2}
    }
    
    \item{\code{"cressie"}}{weights as suggested by Cressie (1985).
      \deqn{LOSS(\theta) = \sum_k n_k [\frac{\hat{\gamma}_k -
	  \gamma_k(\theta)}{\gamma_k(\theta)}]^2}{LOSS(theta) = \sum_k n_k [(hat(gamma_k) -
	  gamma_k(theta))/{gamma_k(theta)}]^2}
    }
    
    \item{\code{"equal"}}{equal values for the weights. For this case
      the estimation corresponds to the ordinary least squares variogram
      fitting. This is the default option if \code{variog$output.type ==
      "cloud"}.
      \deqn{LOSS(\theta) = \sum_k [(\hat{\gamma}_k) -
	\gamma_k(\theta)]^2}{LOSS(\theta) = \sum_k (\hat(\gamma) -
	\gamma(\theta))^2}
    }
    

    Where \eqn{\theta}{theta} is the vector with the variogram parameters
    and
    for each \eqn{k^{th}}{kth}-bin
    \eqn{n_k}{n_k} is the number of
    pairs, \eqn{(\hat{\gamma}_k)}{hat(gamma_k)} is the
    value of the empirical variogram and
    \eqn{\gamma_k(\theta)}{gamma_k(theta)}
    is the value of the theoretical variogram.
    
    See also Cressie (1993) and Barry, Crowder and Diggle (1997) for further discussions on 
    methods to estimate the variogram parameters.
  }
}

\value{
  An object of the \code{\link{class}} \code{"variomodel"} and \code{"variofit"} which is list with the following components:
  \item{nugget}{value of the nugget parameter. An estimated value if
    \code{fix.nugget = FALSE} or a fixed value if \code{fix.nugget = TRUE}.  }
  \item{cov.pars}{a two elements vector with estimated values of the covariance
    parameters \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi},  respectively.  }
  \item{cov.model}{a string with the name of the correlation function.  }
  \item{kappa}{fixed value of the smoothness parameter.  }
  \item{value}{minimized value of the loss function.  }
  \item{max.dist}{maximum distance considered in the variogram fitting.
  }
  \item{minimisation.function}{minimization function used.  }
  \item{weights}{a string indicating the
    type of weights used for the variogram fitting.  }
  \item{method}{a string indicating the
    type of variogram fitting method (OLS or WLS).  }
  \item{fix.kappa}{logical indicating whether the parameter \eqn{\kappa}{kappa} was
    fixed.  }
  \item{fix.nugget}{logical indicating whether the nugget parameter was
    fixed.  }
  \item{lambda}{transformation parameters inherith from the object
    provided in the argument \code{vario}.  }
  \item{message}{status messages returned by the function.  }  
  \item{call}{the function call.  }
}

\references{
  Barry, J.T., Crowder, M.J. and Diggle, P.J. (1997) Parametric
  estimation of the variogram. \emph{Tech. Report, Dept Maths & Stats,
    Lancaster University}.
  
  Cressie, N.A.C (1985) \emph{Mathematical Geology}. \bold{17}, 563-586. 

  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.
  
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{cov.spatial}} for a detailed description of the
  available correlation (variogram) functions,
  \code{\link{likfit}} for maximum
  and restricted maximum likelihood estimation,
  \code{\link{lines.variomodel}} for graphical output of the fitted
  model. For details on the minimization functions see \code{\link{optim}},
  \code{\link{nlm}} and \code{\link{nls}}.
}

\examples{
vario100 <- variog(s100, max.dist=1)
ini.vals <- expand.grid(seq(0,1,l=5), seq(0,1,l=5))
ols <- variofit(vario100, ini=ini.vals, fix.nug=TRUE, wei="equal")
summary(ols)
wls <- variofit(vario100, ini=ini.vals, fix.nug=TRUE)
summary(wls)
plot(vario100)
lines(wls)
lines(ols, lty=2)

\testonly{
vr <- variog(s100, max.dist=1)
## OLS#
o1 <- variofit(vr, ini = c(.5, .5), fix.nug=TRUE, wei = "equal")
o2 <- variofit(vr, ini = c(.5, .5), wei = "equal")
o3 <- variofit(vr, ini = c(.5, .5), fix.nug=TRUE,
      fix.kappa = FALSE, wei = "equal")
o4 <- variofit(vr, ini = c(.5, .5), fix.kappa = FALSE, wei = "equal")
## WLS
w1 <- variofit(vr, ini = c(.5, .5), fix.nug=TRUE)
w2 <- variofit(vr, ini = c(.5, .5))
w3 <- variofit(vr, ini = c(.5, .5), fix.nug=TRUE, fix.kappa = FALSE)
w4 <- variofit(vr, ini = c(.5, .5), fix.kappa = FALSE)
}

}

\keyword{spatial}

