\name{trend.spatial}

\alias{trend.spatial}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Builds the Trend Matrix}

\description{
  Builds the trend matrix according to the specification
  of the mean part of the model provided by the user. 
}

\usage{
trend.spatial(trend, geodata)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{trend}{specifies the mean part of the model.
    See \code{DETAILS} below.  }
  \item{geodata}{an object of the class \code{geodata} as described in
    \code{\link{as.geodata}}. }
}

\details{
  The implicty model assumes that there is an underlying process
  with mean \eqn{\mu(x)}{mu(x)}, where \eqn{x = (x_1, x_2)}{x = (x1, x2)} denotes the coordinates
  of a spatial location.
  The argument \code{trend} defines the form of the mean with the
  following options:

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
      \deqn{\mu(x)= \beta_0 + \beta_1 x_1 + \beta_2 x_2 + \beta_1 (x_1)^2 +
	\beta_2 (x_2)^2 + \beta_1 x_1 * \beta_2 x_2}{%
	\mu(x)= beta0 + beta1*x1 + beta2*x2 + beta1*(x1)^2 +
	beta2*(x2)^2 + beta1*x1 * beta2*x2.}}
    \item{\code{~ model}}{a model specification. See
      \code{\link{formula}} for further details on how to specify
      a model using formulas. Notice that the model term before
      \code{~} is not necessary. Tipically used to include covariates
      (external trend) in the model.}
  }
  
  The following specifications are equivalent:
  \itemize{
    \item \code{trend = "1st"} and \code{trend = ~ x1 + x2}
    \item \code{trend = "2nd"} and \code{trend = ~ x1 + x2 + x1^2 +
      x2^2 + x1*x2}
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
  }
}

\value{
  An \eqn{n \times p}{n x p} trend matrix where \eqn{n}
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

\note{This is an auxiliary function called by
  other \bold{geoR} functions.  }

%\seealso{}
%
%\examples{
%
%}

\keyword{spatial}%-- one or more ...

