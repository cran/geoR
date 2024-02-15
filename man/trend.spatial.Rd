\name{trend.spatial}

\alias{trend.spatial}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Builds the Trend Matrix}

\description{
  Builds the \emph{trend} matrix in accordance to a specification
  of the mean provided by the user. 
}

\usage{
trend.spatial(trend, geodata, add.to.trend)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{trend}{specifies the mean part of the model.
    See \code{DETAILS} below.  }
  \item{geodata}{optional. An object of the class \code{geodata} as described in
    \code{\link{as.geodata}}. }
  \item{add.to.trend}{optional. Specifies aditional terms to the mean
    part of the model. See details below.}
}

\details{
  The implicity model assumes that there is an underlying process
  with mean \eqn{\mu(x)}{mu(x)}, where \eqn{x = (x_1, x_2)}{x = (x1, x2)} denotes the coordinates
  of a spatial location.
  The argument \code{trend} defines the form of the mean and the
  following options are allowed:
  
  \describe{
    \item{\code{"cte"}}{the mean is assumed to be constant over the region,
      in which case \eqn{\mu(x)= \mu}{mu(x) = mu}. This is the default
      option.  }
    \item{\code{"1st"}}{the mean is assumed to be a first order polynomial
      on the coordinates:
      \deqn{\mu(x)= \beta_0 + \beta_1 x_1 + \beta_2 x_2}{%
	mu(x) = beta0 + beta1*x1 + beta2*x2. }}
    \item{\code{"2nd"}}{the mean is assumed to be a second order polynomial
      on the coordinates:
      \deqn{\mu(x)= \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \beta_3 (x_1)^2 +
	\beta_4 (x_2)^2 + \beta_5 x_1 *  x_2}{%
	\mu(x)= beta0 + beta1*x1 + beta2*x2 + beta3*(x1)^2 +
	beta4*(x2)^2 + beta5*x1*x2.}}
    \item{\code{~ model}}{a model specification. See
      \code{\link{formula}} for further details on how to specify
      a model in \R  using formulas. Notice that the model term before
      the \code{~} is not necessary.
      Typically used to include covariates
      (external trend) in the model.}
  }
  
  Denote by \eqn{x_1} and \eqn{x_2} the spatial coordinates.
  The following specifications are equivalent:
  \itemize{
    \item \code{trend = "1st"} and \code{trend = ~ x1 + x2}
    \item \code{trend = "2nd"} and \code{trend = ~ x1 + x2 + I(x1^2) + I(x2^2) + I(x1*x2)}
  }
  
  \bold{Search path for covariates}\cr
  Typically, functions in the package \pkg{geoR} which calls
  \code{trend.spatial} will have the arguments \code{geodata},
  \code{coords} and \code{data}. 
  
  When the trend is specifed as \code{trend = ~ model}
  the terms included in the model will be searched for in the following
  path sequence (modified in version 1.7-6, no longer attach objects):
  \enumerate{
    \item in the object \code{geodata} (coerced to data-frame)
    \item in the users/session Global environment
    \item in the session search path
  }

  The argument \code{add.to.trend} adds terms to what is specified in
  the argument \code{trend}. This seems redundant but allow
  specifications of the type: \code{trend="2nd", add.trend=~other.covariates}.
  
}

\value{
  An object of the class \code{trend.spatial}
  which is an \eqn{n \times p}{n x p} \emph{trend}
  matrix, where \eqn{n}
  is the number of spatial
  locations and \eqn{p} is the number of mean parameters in the model. 
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR/}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\note{This is an auxiliary function typically called by
  other \pkg{geoR} functions.
}

\seealso{The section \code{DETAILS} in the documentation for
  \code{\link{likfit}} for more about the underlying model.
}

\examples{
# a first order polynomial trend
trend.spatial("1st", sic.100)[1:5,]
# a second order polynomial trend
trend.spatial("2nd", sic.100)[1:5,]
# a trend with a covariate
trend.spatial(~altitude, sic.100)[1:5,]
# a first degree trend plus a covariate
trend.spatial(~coords+altitude, sic.100)[1:5,]
# with produces the same as
trend.spatial("1st", sic.100, add=~altitude)[1:5,]
# and yet another exemple
trend.spatial("2nd", sic.100, add=~altitude)[1:5,]
}

\keyword{spatial}%-- one or more ...

