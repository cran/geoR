\name{plot.xvalid}

\alias{plot.xvalid}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Plot Cross-Validation Results}

\description{
  This function produces ten plots with the results produced by the
  cross-validation function \code{\link{xvalid}}.
}

\usage{
\method{plot}{xvalid}(x, coords, borders = NULL, ask = TRUE,
            error = TRUE, std.error = TRUE, data.predicted = TRUE,
            pp = TRUE, map = TRUE, histogram = TRUE,
            error.predicted = TRUE, error.data = TRUE, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{an object of the class \code{"xvalid"}, typically an output
    from the function \code{\link{xvalid}}.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} object containing
    coordinates of the (cross-)validation locations.  }
  \item{borders}{optional. Takes a two column matrix or data-frame with
    coordinates of the borders. If provided the borders are included in
    the errors maps.  }
  \item{ask}{logical. Defines whether or not the user is prompted before
    each plot is produced.}
  \item{error}{logical. Defines whether the plots for the errors
    (\eqn{error = data - predicted}) will
    be produced.  }
  \item{std.error}{logical. Defines whether the plots for the
    standardised errors will be produced.  }
  \item{data.predicted}{logical defining whether a plot of data versus
    predicted should be displayed. Defaults to \code{TRUE}.  }
  \item{pp}{logical defining whether a \emph{pp} plot
    should be displayed. Defaults to \code{TRUE}.  }
  \item{map}{logical defining whether a map of the errors
    should be displayed. Defaults to \code{TRUE}.  }
  \item{histogram}{logical defining whether a histogram of the errors
    should be displayed. Defaults to \code{TRUE}.  }
  \item{error.predicted}{logical defining whether a plot of errors versus
    predicted should be displayed. Defaults to \code{TRUE}.  }
  \item{error.data}{logical defining whether a plot of errors versus
    data should be displayed. Defaults to \code{TRUE}.  }
  \item{\dots}{other arguments to be passed to the function
    \code{\link{plot}}.  }

}

\details{
  The number of plots to be produced will depend on the input options.
  If the graphics device is set to just one plot (something equivalent
  to \option{par(mfcol=c(1,1))}) after each graphic being displayed the user will be prompt
  to press \kbd{<return>} to see the next graphic.

  Alternativaly the user can set the graphical parameter to have several
  plots in one page. With default options for the arguments the maximum
  number of plots (10) is produced and setting \option{par(mfcol=c(5,2))})
  will display them in the same page.

  The \dQuote{errors} for the plots are defined as
  \deqn{error = data - predicted}
  and the plots uses the color blue to indicate positive errors and
  red to indicate negative erros. 
}

\value{
  No value returned. Plots are produced on the current graphics device.  
}

%\note{ ~~further notes~~ }

% ~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{xvalid}} for the cross-validation computations.  }

\examples{
wls <- variofit(variog(s100, max.dist = 1), ini = c(.5, .5), fix.n = TRUE)
xvl <- xvalid(s100, model = wls)
#
op <- par(no.readonly = TRUE)
par(mfcol = c(3,2))
par(mar = c(3,3,0,1))
par(mgp = c(2,1,0))
plot(xvl, error = FALSE, ask = FALSE)
plot(xvl, std.err = FALSE, ask = FALSE)
par(op)
}

\keyword{spatial}%-- one or more ...
\keyword{dplot}
