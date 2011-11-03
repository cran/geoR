\name{lines.variomodel}

\alias{lines.variomodel}
\alias{lines.variomodel.default}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Adds a Line with a Variogram Model to a Variogram Plot}

\description{
  This function adds a line with a variogram model specifyed by the user
  to a current variogram plot.
  The variogram is specifyed either by passing a list with values for the variogram
  elements or using each argument in the function.
}


\usage{
\method{lines}{variomodel}(x, \dots)
\method{lines.variomodel}{default}(x, cov.model, cov.pars, nugget, kappa,
                          max.dist, scaled = FALSE, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{a list with the values for the following components: \code{cov.model}, \code{cov.pars},
    \code{nugget}, \code{kappa} , \code{max.dist}; or a numeric vector
    with values for x-axis values for the variogram (distances).
    This argument is not required if the other arguments in the function
    are provided, otherwise is compulsory.
    If a list is provided the arguments which match the list elements
    are ignored.
  }
  \item{cov.model}{a string with the type of the variogram function. See
    documentation of \code{\link{cov.spatial}} for further details.
  }
  \item{cov.pars}{a vector or matrix with the values for the partial sill
    (\eqn{\sigma^2}{sigmasq}) and range (\eqn{\phi}{phi}) parameters. }
  \item{nugget}{a scalar with the value of the nugget
    (\eqn{\tau^2}{tausq}) parameter.
  }
  \item{kappa}{a scalar with the value of the smoothness
    (\eqn{\kappa}{kappa}) parameters. Only required if \code{cov.model}
    is one of the following: \code{"matern"}, \code{"powered.exponential"},
    \code{"cauchy"} and \code{"gneiting.matern"}
  }
  \item{max.dist}{maximum distance (x-axis) to compute and draw the line
    representing the  variogram model.
    If a list is provided in \code{x} the default is the distance given by
    \code{x$max.dist}. If a vector is provided in \code{x} it takes \code{max(x)}.
  }
  \item{scaled}{logical. If \code{TRUE} the total sill in the plot is
    equals to \eqn{1}.  }
  \item{\dots}{arguments to be passed to the function
    \code{\link{curve}}.  }
}

\details{
  Adds a line with a variogram model  to a plot.
  In conjuction with \code{\link{plot.variogram}} can be
  used for instance to compare sample variograms against fitted models returned by
  \code{\link{variofit}}  and/or \code{\link{likfit}}.
}

\value{
  A line with a variogram model is added to a plot on the current
  graphics device. No values are returned.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.}

%\note{ ~~further notes~~ }

%~Make other sections like WARNING with \section{WARNING }{....} ~

\seealso{\code{\link{lines.variomodel.krige.bayes}},
  \code{\link{lines.variomodel.grf}},
  \code{\link{lines.variomodel.variofit}},
  \code{\link{lines.variomodel.likGRF}},
  \code{\link{plot.variogram}}, \code{\link{lines.variogram}},
  \code{\link{variofit}}, \code{\link{likfit}}, \code{\link{curve}}. 
}

\examples{
# computing and ploting empirical variogram
vario <- variog(s100, max.dist = 1)
plot(vario)
# estimating parameters by weighted least squares
vario.wls <- variofit(vario, ini = c(1, .3), fix.nugget = TRUE)
# adding fitted model to the plot  
lines(vario.wls)
#
# Ploting different variogram models
plot(0:1, 0:1, type="n")
lines.variomodel(cov.model = "exp", cov.pars = c(.7, .25), nug = 0.3, max.dist = 1) 
# an alternative way to do this is:
my.model <- list(cov.model = "exp", cov.pars = c(.7, .25), nugget = 0.3,
max.dist = 1) 
lines.variomodel(my.model, lwd = 2)
# now adding another model
lines.variomodel(cov.m = "mat", cov.p = c(.7, .25), nug = 0.3,
                 max.dist = 1, kappa = 1, lty = 2)
# adding the so-called "nested" models
# two exponential structures
lines.variomodel(seq(0,1,l=101), cov.model="exp",
                 cov.pars=rbind(c(0.6,0.15),c(0.4,0.25)), nug=0, col=2)
## exponential and spherical structures
lines.variomodel(seq(0,1,l=101), cov.model=c("exp", "sph"),
                 cov.pars=rbind(c(0.6,0.15), c(0.4,0.75)), nug=0, col=3)
}
\keyword{spatial}
\keyword{aplot}

