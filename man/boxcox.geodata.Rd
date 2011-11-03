\name{boxcox.geodata}

\alias{boxcox.geodata}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Box-Cox transformation for geodata objects}
\description{
  Method for Box-Cox transformation for objects of the class
  \code{geodata} assuming the data are independent. 
  Computes and optionally plots profile log-likelihoods for the parameter of the Box-Cox simple power transformation \eqn{y^lambda}.
}
\usage{
\method{boxcox}{geodata}(object, trend = "cte", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{an object of the class geodata. See \code{\link{as.geodata}}. }
  \item{trend}{specifies the mean part of the model. See
    \code{\link{trend.spatial}} for further details.
    Defaults to \code{"cte"}.  }
  \item{\dots}{arguments to be passed for the function
    \code{\link[MASS]{boxcox}}.  }
}

\details{
  This is just a wrapper for the function \code{\link[MASS]{boxcox}}
  facilitating its usage with \code{geodata} objects.

  Notice this assume independent observations which is typically
  not the case for \code{geodata} objects.
}

\value{
  A list of the \code{lambda} vector and the computed profile log-likelihood vector, invisibly if the result is plotted.
}

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link[MASS]{boxcox}} for
  parameter estimation results for independent data and
  \code{\link{likfit}} for parameter estimation
  within the geostatistical model.
}

\examples{
if(require(MASS)){
boxcox(wolfcamp)

data(ca20)
boxcox(ca20, trend = ~altitude)
}
}
\keyword{regression}
\keyword{models}
\keyword{hplot}


