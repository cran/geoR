\name{xvalid}

\alias{xvalid}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Cross-validation using kriging}

\description{
  This functions performs \emph{leaving-one-out} cross-validation for geostatistical data.
  Each data location is removed from the data set and the variable at
  this location is predicted using the remaining locations.
  The prediction is based on a model previously fitted. 
}

\usage{
xvalid(geodata, coords = geodata$coords, data = geodata$data,
       model, reestimate = FALSE, variog.obj = NULL,
       output.reestimate = FALSE, locations.xvalid = "all",
       messages.screen = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing element \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data-set.
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
    \code{\link{wlsfit}} or \code{\link{olsfit}}.  }
  \item{reestimate}{logical. Indicates whether or not the model
    parameters should be re-estimated for each point removed from the data-set. }
  \item{variog.obj}{on object with the empirical variogram, typically an
    output  of the function \code{\link{variog}}. Only used if
    \code{reestimate = TRUE} and the object passed to the argument
    \code{model} is the result of a variogram based estimation, i.e. if
    the model was fitted by \code{\link{wlsfit}} or \code{\link{olsfit}}.  }
  \item{output.reestimate}{logical. Only valid if \code{reestimate =
      TRUE}. Specifies whether the re-estimated parameters are returned.  }
  \item{locations.xvalid}{a vector with numbers indicating at which
    locations the cross-validation should be performed. Default to all
    locations.  }
  \item{messages.screen}{logical. Indicates whether
    status messages should be printed on the screen (or output device)
    while the function is running.  }
  \item{\dots}{further arguments to the minimization functions used by \code{\link{likfit}},
    \code{\link{wlsfit}} or \code{\link{olsfit}}.  }
}

\details{
  The cross-validation uses the function \code{krige.conv} to predict
  at each location.
  
  For models fitted by \code{\link{wlsfit}} or \code{\link{olsfit}} the
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
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{plot.xvalid}} for plotting of the results, \code{\link{likfit}}, \code{\link{wlsfit}},
  \code{\link{olsfit}} for parameter estimation and
  \code{\link{krige.conv}} for the kriging method used for predictions.
}


\examples{
if(is.R()) data(s100)
#
# Maximum likelihood estimation
#
s100.ml <- likfit(s100, ini = c(.5, .5), fix.nug = TRUE)
#
# Weighted least squares estimation
#
s100.var <- variog(s100, max.dist = 1)
s100.wls <- wlsfit(s100.var, ini = c(.5, .5), fix.nug = TRUE)
#
# Now, performing cross-validation
#
s100.xv.ml <- xvalid(s100, model = s100.ml)
s100.xv.wls <- xvalid(s100, model = s100.wls)
}
\keyword{spatial}%-- one or more ...
