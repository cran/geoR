\name{variog4}
\alias{variog4}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computes Directional Variograms}
\description{
  Computes directional variograms for 4 directions provided by the user.
}
\usage{
variog4(geodata, coords = geodata$coords, data = geodata$data,
        uvec = "default", trend = "cte", lambda = 1,
        option = c("bin", "cloud", "smooth"),
        estimator.type = c("classical", "modulus"),
        nugget.tolerance = 0, max.dist = NULL, pairs.min = 2,
        bin.cloud = FALSE, direction = c(0, pi/4, pi/2, 3*pi/4),
        tolerance = pi/8, unit.angle = c("radians", "degrees"),
        messages.screen = TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  Most arguments are the same as for the function \code{\link{variog}}.
  The only exception is the following.
  \item{direction}{a vector with values of 4 angles, indicating the
    directions for which the variograms will be computed. Defaults to
    \code{c(0, 45, 90, 135)} (degrees). }
}

%\details{
%  ~~ If necessary, more details than the __description__  above ~~
%}

\value{
  The output is an object of the class \code{variog4},
  a list with five components.
  The first four elements are estimated variograms for the directions
  provided and the last is the omnidirectional variogram. 
  Each individual component has is an object of the class \code{variogram},
  an output of the function \code{\link{variog}}.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{variog}} for variogram calculations and
  \code{\link{plot.variog4}} for plotting results}

\examples{
if(is.R()) data(s100)
plot(variog4(s100, max.dist=1))
}
\keyword{spatial}%-- one or more ...
