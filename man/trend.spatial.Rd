\name{trend.spatial}

\alias{trend.spatial}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Builds the Trend Matrix}

\description{
  Builds the \emph{trend} matrix in accordance to a specification
  of the mean provided by the user. 
}

\usage{
trend.spatial(trend, geodata)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{trend}{specifies the mean part of the model.
    See \code{DETAILS} below.  }
  \item{geodata}{optional. An object of the class \code{geodata} as described in
    \code{\link{as.geodata}}. }
}

\details{
  The implicity model assumes that there is an underlying process
  with mean \eqn{\mu(x)}{mu(x)}, where \eqn{x = (x_1, x_2)}{x = (x1, x2)} denotes the coordinates
  of a spatial location.
  The argument \code{trend} defines the form of the mean with the
  following options for usage:

  \itemize{
    \item{\code{"cte"}}{the mean is assumed to be constant over the region,
      in which case \eqn{\mu(x)= \mu}{mu(x) = mu}. This is the default
      option.  }
    \item{\code{"1st"}}{the mean is assumed to be a first degree polynomial
      on the coordinates:
      \deqn{\mu(x)= \beta_0 + \beta_1 x_1 + \beta_2 x_2}{%
	mu(x) = beta0 + beta1*x1 + beta2*x2. }}
    \item{\code{"2nd"}}{the mean is assumed to be a second degree polynomial
      on the coordinates:
      \deqn{\mu(x)= \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \beta_3 (x_1)^2 +
	\beta_4 (x_2)^2 + \beta_5 x_1 *  x_2}{%
	\mu(x)= beta0 + beta1*x1 + beta2*x2 + beta3*(x1)^2 +
	beta4*(x2)^2 + beta5*x1*x2.}}
    \item{\code{~ model}}{a model specification. See
      \code{\link{formula}} for further details on how to specify
      a model using formulas. Notice that the model term before
      \code{~} is not necessary. Tipically used to include covariates
      (external trend) in the model.}
  }
  
  Denote by \eqn{x_1} and \eqn{x_2} the spatial coordinates.
  The following specifications are equivalent:
  \itemize{
    \item \code{trend = "1st"} and \code{trend = ~ x1 + x2}
    \item \code{trend = "2nd"} and \code{trend = ~ x1 + x2 + x1^2 + x2^2 + x1*x2}
  }

  \bold{Search path for covariates}\cr
  Typically, functions in the package \code{geoR} which calls
  \code{trend.spatial} will have the arguments \code{goedata},
  \code{coords} and \code{data}. 

  When the trend is specifed as \code{trend = ~ model}
  the terms included in the model will be searched for in the following
  loactions (in this order):
  \enumerate{
    \item as elements of the list \code{geodata}
    \item as columns in the data-frame \code{geodata$covariates}
    \item as columns in the data-frame \code{geodata$data}
    \item in the search path
  }
}

\value{
  An object of the class \code{trend.spatial}
  which is an \eqn{n \times p}{n x p} \emph{trend}
  matrix, where \eqn{n}
  is the number of spatial
  locations and \eqn{p} is the number of mean parameters in the model. 
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\note{This is an auxiliary function typically called by
  other \bold{geoR} functions.  }

\seealso{The section \code{DETAILS} in \code{\link{likfit}} for more
  on the underlying model.}

\examples{
data(SIC)
# a first degree polynomial trend
trend.spatial("1st", sic.100)[1:5,]
# a second degree polynomial trend
trend.spatial("2nd", sic.100)[1:5,]
# a trend with a covariate
trend.spatial(~altitude, sic.100)[1:5,]
# a first degree trend plus a covariate
trend.spatial(~coords+altitude, sic.100)[1:5,]
}

\keyword{spatial}%-- one or more ...

