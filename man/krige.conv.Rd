\name{krige.conv}

\alias{krige.conv}
\alias{krige.control}
%- Also NEED an `\alias' for EACH other topic documented here.

\concept{kriging}
\concept{spatial interpolation}

\title{Spatial Prediction -- Conventional Kriging}


\description{
  This function performs spatial prediction for fixed covariance
  parameters using global neighbourhood.
  
  Options available implement the following types of kriging:
  \emph{SK} (simple kriging), \emph{OK} (ordinary kriging),
  \emph{KTE} (external trend kriging) and \emph{UK} (universal kriging). 
  
}

\usage{
krige.conv(geodata, coords=geodata$coords, data=geodata$data,
           locations, borders, krige, output)

krige.control(type.krige = "ok", trend.d = "cte", trend.l = "cte",
            obj.model = NULL, beta, cov.model, cov.pars, kappa,
            nugget, micro.scale = 0, dist.epsilon = 1e-10, 
            aniso.pars, lambda)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. Typically an object of the \code{\link{class}}
    \code{"geodata"} - a \pkg{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix or data-frame with the 2-D
    coordinates of the \eqn{n} data locations.
    By default it takes the
    component \code{coords} of the argument \code{geodata}, if provided.  }
  \item{data}{a vector with \emph{n} data values. By default it takes the
    component \code{data} of the argument \code{geodata}, if provided.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data-frame with the 2-D
    coordinates of the \eqn{N} prediction locations, or a list for which
    the first two components are used. Input is internally checked by the
    function \code{check.locations}.  }
  \item{borders}{optional. By default reads the element \code{borders}
    from the \code{geodata} object, if present.
    Setting to \code{NULL} prevents this behavior.
    If a two column matrix defining a polygon is
    provided the prediction is performed only at locations inside this polygon. }
  \item{krige}{a list defining the model components and the type of
    kriging. It can take an output to a call to \code{krige.control} or
    a list with elements as for the arguments in \code{krige.control}.
    Default values are assumed for arguments  or list elements not provided.
    See arguments for `krige.control'.  }
  \item{output}{a list specifying output options.
    It can take an output to a call to \code{output.control} or
    a list with elements as for the arguments in \code{output.control}.
    Default values are assumed for arguments not provided.
    See documentation for
    \code{\link{output.control}} for further details.  }
  \item{type.krige}{type of kriging to be performed. Options are
    \code{"SK", "OK"} corresponding to simple or ordinary
    kriging. Kriging with external trend and universal kriging can be
    defined setting \code{type.krige = "OK"} and specifying the
    trend model using the arguments \code{trend.d} and
    \code{trend.l}.  } 
  \item{trend.d}{specifies the trend (covariate) values at the data
    locations.
    See documentation of \code{\link{trend.spatial}} for
    further details.
    Defaults to \code{"cte"}.   }
  \item{trend.l}{specifies the trend (covariate) values at prediction
    locations. It must be of the same type as for \code{trend.d}.
    Only used if prediction locations are provided in the argument
    \code{locations}.  }
  \item{obj.model}{a list with the model parameters. Typically an
    output of \code{\link{likfit}} or
    \code{\link{variofit}}.  }
  \item{beta}{numerical value of the mean (vector) parameter.
    Only used if \code{type.krige="SK"}.}  
  \item{cov.model}{string indicating the name of the model for the
    correlation function. Further details  can be found in the
    documentation of the function
    \code{\link{cov.spatial}}.  }
  \item{cov.pars}{a 2 elements vector with values of the covariance parameters \eqn{\sigma^2}{sigma^2}
    (partial sill)  and \eqn{\phi}{phi} (range parameter),
    respectively.  }
 %  \item{cov.pars}{a vector with 2 elements or an \eqn{n \times 2}{n x
%      2} matrix with the covariance parameters \eqn{\sigma^2}{sigma^2}
%    (partial sill)  and \eqn{\phi}{phi} (range parameter). If a vector, the elements are the values of
%    \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.
%    If a matrix, corresponding to a model with several structures, the
%    values of \eqn{\sigma^2}{sigma^2} 
%    are in the first column and the values of \eqn{\phi}{phi} are in the second.  }
  \item{kappa}{additional smoothness parameter required by the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} and
    \code{"gneiting.matern"}.   }
  \item{nugget}{the value of the nugget variance parameter \eqn{\tau^2}{tau^2}. Defaults to zero.  }
  \item{micro.scale}{micro-scale variance. If different from zero, the
    nugget variance is divided into 2 terms: \emph{micro-scale variance}
    and \emph{measurement error}. This affect the precision of the predictions. 
    Often in practice, these two variance components are indistinguishable but the
    distinction can be made here if justifiable. See the section
    \code{DETAILS} in the documentation of \code{\link{output.control}}. } 
  \item{dist.epsilon}{a numeric value. Locations which are separated
    by a distance less than this value
    are considered co-located.  }
  \item{aniso.pars}{parameters for geometric anisotropy
    correction. If \code{aniso.pars = FALSE} no correction is made, otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link{coords.aniso}}.  }
  \item{lambda}{numeric value of the Box-Cox transformation parameter.
    The value \eqn{\lambda = 1}{lambda = 1} corresponds to
    no transformation and \eqn{\lambda = 0}{lambda = 0} corresponds to
    the log-transformation. 
    Prediction results are back-transformed and
    returned is the same scale as for the original data.  }
}

\details{
  According to the arguments provided,  one of the following different
  types of kriging: \emph{SK}, \emph{OK}, \emph{UK} or  \emph{KTE}
  is performed.
  Defaults correspond to ordinary kriging.

}


\value{
  An object of the \code{\link{class}} \code{kriging}.
  The attribute \code{prediction.locations} containing the name of the
  object with the coordinates of the prediction locations (argument
  \code{locations}) is assigned to the object.
  Returns  a  list with the following components:
  
  \item{predict}{a vector with predicted values.  }
  \item{krige.var}{a vector with predicted variances.  }
  \item{beta.est}{estimates of the \eqn{\beta}{beta}, parameter
    implicit in kriging procedure. Not included if \code{type.krige = "SK"}.  }
  \item{simulations}{an \eqn{ni \times n.sim}{ni x n.sim} matrix where \eqn{ni} is the
    number of prediction locations. Each column corresponds to a
    conditional simulation of the predictive distribution.
    Only returned if \code{n.sim > 0}.  }
  \item{message}{messages about the type of prediction performed.  }
  \item{call}{the function call.  }
  Other results can be included depending on the options passed to
  \code{\link{output.control}}. 
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR/}.
}

%\note{}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{output.control}} sets output options,
  \code{\link{image.kriging}} and \code{\link{persp.kriging}}
  for graphical output of the results,
  \code{\link{krige.bayes}} for Bayesian prediction and \code{\link{ksline}}
  for a different implementation of kriging allowing for moving
  neighborhood. For model fitting see \code{\link{likfit}} or \code{\link{variofit}}.
}


\examples{
\dontrun{
# Defining a prediction grid
loci <- expand.grid(seq(0,1,l=21), seq(0,1,l=21))
# predicting by ordinary kriging
kc <- krige.conv(s100, loc=loci,
                 krige=krige.control(cov.pars=c(1, .25)))
# mapping point estimates and variances
par.ori <- par(no.readonly = TRUE)
par(mfrow=c(1,2), mar=c(3.5,3.5,1,0), mgp=c(1.5,.5,0))
image(kc, main="kriging estimates")
image(kc, val=sqrt(kc$krige.var), main="kriging std. errors")
# Now setting the output to simulate from the predictive
# (obtaining conditional simulations),
# and to compute quantile and probability estimators
s.out <- output.control(n.predictive = 1000, quant=0.9, thres=2)
set.seed(123)
kc <- krige.conv(s100, loc = loci,
         krige = krige.control(cov.pars = c(1,.25)),
         output = s.out)
par(mfrow=c(2,2))
image(kc, val=kc$simul[,1], main="a cond. simul.")
image(kc, val=kc$simul[,1], main="another cond. simul.")
image(kc, val=(1 - kc$prob), main="Map of P(Y > 2)")
image(kc, val=kc$quant, main="Map of y s.t. P(Y < y) = 0.9")
par(par.ori)
}
}
\keyword{spatial}

