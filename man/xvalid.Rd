\name{xvalid}

\alias{xvalid}
%- Also NEED an `\alias' for EACH other topic documented here.

\concept{cross-validation}


\title{Cross-validation using kriging}

\description{
  This is a function to perform model validation.
  Options include \emph{leaving-one-out} cross-validation where.
  each data location is removed from the data set and the variable at
  this location is predicted using the remaining locations, for as given
  model. This can be
  done for all or some of the locations.
  Alternativelly, other validation locations which are not the same as the
  original data locations can be used.   
}

\usage{
xvalid(geodata, coords = geodata$coords, data = geodata$data,
       model, reestimate = FALSE, variog.obj = NULL,
       output.reestimate = FALSE, locations.xvalid = "all",
       data.xvalid = NULL, messages, \dots)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing element \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \pkg{geoR} data-set.
    If not provided the arguments
    \code{coords} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing
    coordinates of the \eqn{n} data locations in each row.
    Defaults to \code{geodata$coords}, if provided.}
  \item{data}{a vector or matrix with data values.
    If a matrix is provided, each column is regarded as one variable or realization.
    Defaults to \code{geodata$data}, if provided.}
  \item{model}{an object containing information on a fitted
    model. Typically an output of \code{\link{likfit}},
    \code{\link{variofit}}.  }
  \item{reestimate}{logical. Indicates whether or not the model
    parameters should be re-estimated for each point removed from the data-set. }
  \item{variog.obj}{on object with the empirical variogram, typically an
    output  of the function \code{\link{variog}}. Only used if
    \code{reestimate = TRUE} and the object passed to the argument
    \code{model} is the result of a variogram based estimation, i.e. if
    the model was fitted by \code{\link{variofit}}.  }
  \item{output.reestimate}{logical. Only valid if \code{reestimate =
      TRUE}. Specifies whether the re-estimated parameters are returned.  }
  \item{locations.xvalid}{there are three possible specifications for
    this argument: \code{"all"} indicates the \emph{leaving-on-out}
    method is used. The second possibility is to use only a sub-set of
    the data for cross-validation. For this case  a vector should be
    provided with numbers indicating at which locations the
    cross-validation should be performed. The third option is to perform
    validation on a different set of data. For this a matrix with the
    coordinates of the validation points should be provided and the the
    argument \code{locations.data} should also be provided.  }
  \item{data.xvalid}{data values at the validation locations.
    Only used if the validation locations are not the same or a
    subset of the original data coordinates.  }
  \item{messages}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{further arguments to the minimization functions used by
    \code{\link{likfit}}, \code{\link{variofit}}.  }
}

\details{
  The cross-validation uses the function \code{krige.conv} to predict
  at each location.
  
  For models fitted by \code{\link{variofit}} the
  parameters \eqn{\kappa}{kappa}, \eqn{\psi_A}{psiA}, \eqn{\psi_R}{psiR}
  and \eqn{\lambda}{lambda} are always regarded as fixed.
  
  See documentation of the function \code{\link{likfit}} for more
  details on the model and its parameters.
}

\value{
  An object of the \code{\link{class}} \code{"xvalid"} which is a list with the following components:
  \item{data}{the original data.  }
  \item{predicted}{the values predicted by cross-validation.  }
  \item{krige.var}{the cross-validation prediction variance.  }
  \item{error}{difference \code{data - predicted}.   }
  \item{std.error}{the errors divided by the square root of the
    prediction variances.  }
  \item{prob}{the cumulative probability at original value under a
    normal distribution with parameters given by the cross-validation results.  }
  
  If \code{reestimate = TRUE} and \code{output = TRUE} additional
  columns are added to the data-frame. Each column will contain the
  values of the re-estimated parameters.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.est.ufpr.br/geoR}.
}

\author{Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{plot.xvalid}} for plotting of the results, \code{\link{likfit}},
  \code{\link{variofit}} for parameter estimation and
  \code{\link{krige.conv}} for the kriging method used for predictions.
}


\examples{
data(s100)
#
# Maximum likelihood estimation
#
s100.ml <- likfit(s100, ini = c(.5, .5), fix.nug = TRUE)
#
# Weighted least squares estimation
#
s100.var <- variog(s100, max.dist = 1)
s100.wls <- variofit(s100.var, ini = c(.5, .5), fix.nug = TRUE)
#
# Now, performing cross-validation without reestimating the model
#
s100.xv.ml <- xvalid(s100, model = s100.ml)
s100.xv.wls <- xvalid(s100, model = s100.wls)
##
## Plotting results
##
par.ori <- par(no.readonly = TRUE)
##
par(mfcol=c(5,2), mar=c(2.3,2.3,.5,.5), mgp=c(1.3, .6, 0))
plot(s100.xv.ml)
par(mfcol=c(5,2))
plot(s100.xv.wls)
##
par(par.ori)
#
\testonly{
data(s100)
  vr <- variog(s100, max.dist=1)
  ## OLS#
  o1 <- variofit(vr, ini=c(.5, .5), fix.nug=TRUE, wei="equal")
  xvo1 <- xvalid(s100, model=o1, variog.obj=vr, loc=sample(1:100,5))
  o2 <- variofit(vr, ini=c(.5, .5), wei = "equal")
  xvo2 <- xvalid(s100, model=o2, variog.obj=vr, loc=sample(1:100,5))
  o3 <- variofit(vr, ini=c(.5, .5), fix.nug=TRUE, fix.kappa = FALSE, wei = "equal")
  xvo3 <- xvalid(s100, model=o3, variog.obj=vr, loc=sample(1:100,5))
  #o4 <- variofit(vr, ini=c(.5, .5), fix.kappa = FALSE, wei = "equal")
  #xvo4 <- xvalid(s100, model=o4, variog.obj=vr, loc=sample(1:100,5))
  ## WLS
  w1 <- variofit(vr, ini=c(.5, .5), fix.nug=TRUE)
  xvw1 <- xvalid(s100, model=w1, variog.obj=vr, loc=sample(1:100,5))
  w2 <- variofit(vr, ini=c(.5, .5))
  xvw2 <- xvalid(s100, model=w2, variog.obj=vr, loc=sample(1:100,5))
  w3 <- variofit(vr, ini=c(.5, .5), fix.nug=TRUE, fix.kappa = FALSE)
  xvw3 <- xvalid(s100, model=w3, variog.obj=vr, loc=sample(1:100,5))
  w4 <- variofit(vr, ini=c(.5, .5), fix.kappa = FALSE)
  xvw4 <- xvalid(s100, model=w4, variog.obj=vr, loc=sample(1:100,5))
  # ML
  m1 <- likfit(s100, ini=c(.5,.5), fix.nug=FALSE, cov.model="mat", kappa=1)
  xvm1 <- xvalid(s100, model=m1, loc=sample(1:100,5))
  ap <- grf(10, cov.pars=c(1, .25))
  xvm2 <- xvalid(s100, model=m1, locations.xvalid=ap$coords, data.xvalid=ap$data)
  bor <- s100$coords[chull(s100$coords),]
  par(mfcol=c(5,2),mar=c(3,3,1,0),mgp=c(2,.5,0))  
  plot(xvm2, borders=bor)
  # RML
  rm1 <- likfit(s100, ini=c(.5,.5), fix.nug=FALSE, met = "REML", cov.model="mat", kappa=1)
  xvrm1 <- xvalid(s100, model=m1, loc=sample(1:100,5))
}
}

\keyword{spatial}%-- one or more ...
