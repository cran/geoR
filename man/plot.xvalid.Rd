\name{plot.xvalid}

\alias{plot.xvalid}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Cross-Validation Results}

\description{
  This function produces ten plots with the results produced by the
  cross-validation function \code{\link{xvalid}}.
}

\usage{
plot.xvalid(obj, valid.obj, coords = valid.obj$coords,
            error.plots = TRUE, std.error.plots = TRUE,
            ask = TRUE)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{an object of the class \code{"xvalid"}, typically an output
    from the function \code{\link{xvalid}}.  }
  \item{valid.obj}{a list containing at least the element \code{coords}
    as described next. Typically an object of the class
    \code{"geodata"} - a \bold{geoR} data-set.
    The argument \code{coords} can be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix containing
    coordinates of the (cross-)validation locations in each row.
    Defaults to \code{valid.obj$coords}, if provided.}
  \item{error.plots}{logical. Defines whether the plots for the errors
    (\eqn{error = data - predicted}) will
    be produced.  }
  \item{error.plots}{logical. Defines whether the plots for the
    standardised errors will
    be produced.  }
  \item{ask}{logical. Defines whether or not the user is prompted before
    each plot is produced.}
}
%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}
\value{
  No value returned. Plots are produced on the current graphics device.
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{xvalid}} for the cross-validation computations.  }

\examples{
if(is.R()) data(s100)
wls <- variofit(variog(s100, max.dist = 1), ini = c(.5, .5), fix.n = TRUE)
xvl <- xvalid(s100, model = wls)
#
op <- par(no.readonly = TRUE)
par(mfcol = c(3,2))
par(mar = c(3,3,0,1))
par(mgp = c(2,1,0))
plot(xvl, s100, error = FALSE, ask = FALSE)
plot(xvl, s100, std.err = FALSE, ask = FALSE)
par(op)


}

\keyword{spatial}%-- one or more ...
\keyword{dplot}
