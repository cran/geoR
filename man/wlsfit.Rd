\name{wlsfit}

\alias{wlsfit}
\alias{loss.wlsvario}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Variogram Fitting by Weighted Least Squares}

\description{
  Fits a variogram model to an empirical variogram.
  Variogram parameters are estimated by weighted least squares.
}

\usage{
wlsfit(vario, ini.cov.pars, cov.model = "matern",
       fix.nugget = FALSE, nugget = 0, kappa = 0.5,
       simul.number = NULL, max.dist = "all",
       minimisation.function = c("optim", "nlm"),
       lower = 0, weight = c("npairs", "cressie"),
       messages.screen = TRUE, ...)

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
    \code{\link{cov.spatial}}.
    Defaults are equivalent to the \emph{exponential} model.  }
  \item{fix.nugget}{logical, indicating whether the parameter
    \eqn{\tau^2}{tau^2} (nugget variance) should be regarded as fixed
    (\code{fix.nugget = TRUE}) or should
    be estimated (\code{fix.nugget = FALSE}). Defaults to
    \code{FALSE}.  }    
  \item{nugget}{value for the nugget parameter.  Regarded as a
    fixed values if \code{fix.nugget = TRUE} or as a initial value for the
    minimization algorithm if \code{fix.nugget = FALSE}. 
    Defaults to zero.  }
  \item{kappa}{fixed value of the smoothness parameter, only required by
    the following correlation functions: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"}
    and \code{"gneiting.matern"}.  }
  \item{simul.number}{number of simulation. To be used when the object passed to the
    argument \code{vario} has empirical variograms for more than one
    data-set (or simulation). Indicates to which one the model will be
    fitted.  }
  \item{max.dist}{maximum distance considered when fitting the
    variogram. Defaults to \code{vario$max.dist}.  }
  \item{minimisation.function}{minimization function used to estimate
    the parameters. Options are \code{"optim"} and \code{"nlm"}.
    Defaults to \code{"optim"}.  }
  \item{lower}{lower limits for the parameter values. Defaults to zero.}
  \item{weight}{type weights used in the loss function. See
    \code{DETAILS} below.} 
  \item{messages.screen}{logical. Indicates
    whether or not status messages are printed on the screen (or other output device)
    while the function is running.  }
  \item{\dots}{further parameters to be passed to the minimization
    function. Typically arguments of the type \code{control()} which controls the
    behavior of the minimization algorithm. See documentation for the 
    selected minimization function for further details.  }
}

\details{
  \bold{Initial values}
  
  The algorithms for minimization functions require initial values of
  the parameters.
  
  A unique initial value is used if a vector is provided in the argument
  \code{ini.cov.pars}. The elements are initial values for
  \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.
  This vector is concatenated with the value of the
  argument \code{nugget} if \code{fix.nugget = FALSE}.
  
  Specification of multiple initial values is also possible.
  If this is the case, the function
  searches for the one which minimizes the loss function and uses this as
  the initial value for the minimization algorithm.
  Multiple initial values are specified  by providing a matrix in the
  argument
  \code{ini.cov.pars} and/or, if \code{fix.nugget = FALSE},
  providing a vector with length greater than one for the argument
  \code{nugget}.
  If \code{ini.cov.pars} is a matrix, the first column has values of
  \eqn{\sigma^2}{sigma^2} and the second has values of  \eqn{\phi}{phi}.\cr


  
  \bold{Weights}
  
  Two different types of weights can be used within the loss function:
  \describe{
    \item{\code{"npairs"}}{indicating that the weights are given by the
      number of pairs in each bin.  }
    
    \item{\code{"cressie"}}{weights
      as suggested by Cressie (1993, p.95).  } 
    
    See also Barry, Crowder and Diggle (1997) for a discussion on the
    methods to estimate variogram parameters.
  }
}


\value{
   An object of the \code{\link{class}} \code{"variomodel"} which is list with the following components:
  \item{nugget}{value of the nugget parameter. An estimated value if
    \code{fix.nugget = FALSE} or a fixed value if \code{fix.nugget = TRUE}.  }
  \item{cov.pars}{a 2 elements vector with estimated values of the covariance
    parameters \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.  }
  \item{cov.model}{a string with the name of the correlation function.  }
  \item{kappa}{fixed value of the smoothness parameter.  }
  \item{value}{minimized value of the loss function.  }
  \item{max.dist}{maximum distance considered in the variogram fitting.
  }
  \item{minimisation.function}{minimization function used.  }
  \item{message}{status messages returned by the function.  }
  \item{method}{a string \code{"WLS"} indicating the
    estimation method.  }
  \item{fix.nugget}{logical indicating whether the nugget parameter was
    fixed.  }
  \item{lambda}{transformation parameters inherith from the object
    provided in the argument \code{vario}.  }
  \item{call}{the function call.  }
 
}

\references{

  Barry, J.T., Crowder, M.J. and Diggle, P.J. (1997) Parametric
  estimation of the variogram. \emph{Tech. Report, Dept Maths & Stats,
    Lancaster University}.

  Cressie, N.A.C (1993) \emph{Statistics for Spatial Data}. New York: Wiley.

  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{cov.spatial}} for a detailed description of the
  available correlation (variogram) functions, \code{\link{olsfit}} for
  ordinary least squares variogram fit, \code{\link{likfit}} for maximum
  and restricted maximum likelihood estimation,
  \code{\link{lines.variomodel}} for graphical output of the fitted
  model. For details on the minimization functions see
  documentation for \code{\link{optim}} and
  \code{\link{nlm}}.
}


\examples{
if(is.R()) data(s100)
vario100 <- variog(s100, max.dist=1)
ini.vals <- expand.grid(seq(0,1,l=5), seq(0,1,l=5))
wls <- wlsfit(vario100, ini=ini.vals, fix.nug=TRUE)
summary(wls)
plot(vario100)
lines(wls)

}

\keyword{spatial}
