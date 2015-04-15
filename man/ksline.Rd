\name{ksline}

\alias{ksline}
\alias{.ksline.aux.1}

\title{Spatial Prediction -- Conventional Kriging}

\description{
  This function performs spatial prediction for given covariance
  parameters. Options implement
  the following kriging types:
  \emph{SK} (simple kriging), \emph{OK} (ordinary kriging),
  \emph{KTE} (external trend kriging) and \emph{UK} (universal kriging).
  
  The function \code{\link{krige.conv}} should be preferred, unless
  moving neighborhood is to be used.
}

\usage{
ksline(geodata, coords = geodata$coords, data = geodata$data,
       locations, borders = NULL, 
       cov.model = "matern",
       cov.pars=stop("covariance parameters (sigmasq and phi) needed"),
       kappa = 0.5, nugget = 0, micro.scale = 0,
       lambda = 1, m0 = "ok", nwin = "full",
       n.samples.backtransform = 500, trend = 1, d = 2,
       ktedata = NULL, ktelocations = NULL, aniso.pars = NULL,
       signal = FALSE, dist.epsilon = 1e-10, messages)
}


\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. Typically an object of the \code{\link{class}}
    \code{"geodata"} - a \pkg{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix where each row has the 2-D
    coordinates of the \eqn{n} data locations.
    By default it takes the
    component \code{coords} of the argument \code{geodata}, if
    provided.  }
  \item{data}{a vector with \emph{n} data values. By default it takes the
    component \code{data} of the argument \code{geodata}, if provided.  }
  \item{locations}{an \eqn{N \times 2}{N x 2} matrix or data-frame with the 2-D
    coordinates of the \eqn{N} prediction locations, or a list for which
    the first two components are used. Input is internally checked by the
    function \code{check.locations}.  } 
  \item{borders}{optional. If a two column matrix defining a polygon is
    provided the prediction is performed only at locations inside this polygon. }
  \item{cov.pars}{a vector with 2 elements or an \eqn{n \times 2}{n x
      2} matrix with the covariance parameters \eqn{\sigma^2}{sigma^2}
    (partial sill)  and \eqn{\phi}{phi} (range parameter). If a vector, the elements are the values of
    \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi}, respectively.
    If a matrix, corresponding to a model with several structures, the
    values of \eqn{\sigma^2}{sigma^2} 
    are in the first column and the values of \eqn{\phi}{phi} are in the second.  }
  \item{nugget}{the value of the nugget variance parameter \eqn{\tau^2}{tau^2}. Defaults to zero.  }
  \item{micro.scale}{micro-scale variance. If different from zero, the
    nugget variance is divided into 2 terms: \emph{micro-scale variance}
    and \emph{measurement error}. This might affect the precision of
    the predictions. 
    In practice, these two variance components are usually indistinguishable but the
    distinction can be made here if justifiable.  } 
  \item{cov.model}{string indicating the name of the model for the
    correlation function. Further details  in the
    documentation for
    \code{\link{cov.spatial}}.  
    Defaults are equivalent to the \emph{exponential} model.  }
  \item{kappa}{additional smoothness parameter required by the following correlation
    functions: \code{"matern"}, 
    \code{"powered.exponential"},  
    \code{"cauchy"}  and  
     \code{ "gneiting.matern"}.   }
  \item{lambda}{numeric value of the Box-Cox transformation parameter.
    The value \eqn{\lambda = 1}{lambda = 1} corresponds to
    no transformation and \eqn{\lambda = 0}{lambda = 0} corresponds to
    the log-transformation. 
    Prediction results are back-transformed and
    returned is the same scale as for the original data.  }
  \item{m0}{The default value \code{"ok"} indicates that ordinary
    kriging will be performed. Other options are \code{"kt"} for kriging
    with a trend model (universal kriging) and \code{"kte"} for kriging
    with external trend (covariates). If a numeric value is provided it
    is assumed to be the value of a know mean and simple kriging is then
    performed.
    If \code{"av"} the arithmetic mean of the data is assumed to be the
    know mean for simple kriging algorithm.  }
  \item{nwin}{If \code{"full"} \emph{global neighborhood} is used
    i.e., all data values are used in the
    prediction of every prediction location.
    An integer number defines the \emph{moving neighborhood} algorithm.
    The number provided is used as  
    the number closest neighbors to be used for the prediction at each
    location.
    Defaults to \code{"full"}.  }
  \item{n.samples.backtransform}{number of samples used in the
    back-transformation. When transformations are used
    (specified by an argument \code{lambda}), back-transformations
    are usually performed by sampling from the predictive distribution and
    then back-transforming the sampled values. The exceptions are for
    \eqn{\lambda = 0}{lambda = 0} (log-transformation)  and
    \eqn{\lambda = 1}{lambda = 1} (no transformation).  }
  \item{trend}{only required if \code{m0 = "kt"} (universal kriging).
    Possible values are \eqn{1} or \eqn{2}, corresponding to a first or second degree
    polynomial trend on the coordinates, respectively.  }
  \item{d}{spatial dimension, \eqn{1} defines a prediction on a line, \eqn{2} on
    a plane (the default).}
  \item{ktedata}{only required if \code{m0 = "kte"}. A vector or matrix
    with the values of the
    external trend (covariates) at the data locations.  }
  \item{ktelocations}{only required if \code{m0 = "kte"}. A vector or
    matrix with the values of the
    external trend (covariates) at the prediction locations.  }
  \item{aniso.pars}{parameters for geometric anisotropy
    correction. If \code{aniso.pars = FALSE} no correction is made, otherwise
    a two elements vector with values for the anisotropy parameters
    must be provided. Anisotropy correction consists of a
    transformation of the data and prediction coordinates performed
    by the function \code{\link{coords.aniso}}.  }
  \item{signal}{logical. If \code{TRUE} the signal is predicted, otherwise the
    variable is predicted. If no transformation is performed the
    expectations are the same in both cases and the difference is only for
    values of the kriging variance, if the value of the nugget is
    different from zero.  }
  \item{dist.epsilon}{a numeric value. Points which are separated
    by a distance less than this value
    are considered co-located.  }
  \item{messages}{logical. Indicates
    whether or not status messages are printed on the screen (or other output device)
    while the function is running.}
}

%\details{
%}

\value{
  An object of the \code{\link{class}} \code{kriging} which is a  list
  with the following components:
  \item{predict}{the predicted values.  }
  \item{krige.var}{the kriging variances.  }
  \item{dif}{the difference between the predicted value and the global
    mean. Represents the contribution to the neighboring data to
    the prediction at each point.  }
  \item{summary}{values of the arithmetic and weighted mean of the data and standard
    deviations. The weighted mean corresponds to the estimated value of
    the global mean. }
  \item{ktrend}{the matrix with trend if \code{m0 = "kt"}
    (universal kriging).  }
  \item{ktetrend}{the matrix with trend if \code{m0 = "kte"}
    (external trend kriging).}
  \item{beta}{the value of the mean which is implicitly estimated for
    \code{m0 = "ok", "kte"} or \code{"kt"}.}
  \item{wofmean}{weight of mean. The predicted value is an
    weighted average between the global mean and the values at the
    neighboring locations.
    The value returned is the weight of the mean.  }
  \item{locations}{the coordinates of the prediction locations.  }
  \item{message}{status messages returned by the algorithm.  }
  \item{call}{the function call.  }
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\note{
  This is a preliminary and inefficient function implementing kriging methods.
  For predictions using global neighborhood the function
  \code{\link{krige.conv}} should be used instead.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{krige.conv}} for a more efficient implementation of
  conventional kriging methods, \cr \code{\link{krige.bayes}} for Bayesian prediction.}


\examples{
loci <- expand.grid(seq(0,1,l=31), seq(0,1,l=31))
kc <- ksline(s100, loc=loci, cov.pars=c(1, .25))
par(mfrow=c(1,2))
image(kc, main="kriging estimates")
image(kc, val=sqrt(kc$krige.var), main="kriging std. errors")
}

\keyword{spatial}%-- one or more ...


