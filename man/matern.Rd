\name{matern}

\alias{matern}
%- Also NEED an `\alias' for EACH other topic documented here.

\title{Computer Values of the Matern Correlation Function}

\description{
  This function computes values of the \eqn{\mbox{Mat\'{e}rn}}{Matern}
  correlation function
  for given distances and parameters.
}

\usage{
matern(u, phi, kappa)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{u}{a vector, matrix or array with values of the distances
    between pairs of data locations.  }
  \item{phi}{value of the range parameter \eqn{\phi}{phi}.  }
  \item{kappa}{value of the smoothness parameter \eqn{\kappa}{kappa}.  }
}

\details{
  The \eqn{\mbox{Mat\'{e}rn}}{Matern} model is defined as:
  
  \deqn{\rho(u;\phi,\kappa) =\{2^{\kappa-1}
    \Gamma(\kappa)\}^{-1} (u/\phi)^\kappa
    K_\kappa(u/\phi)}{rho(u;phi,kappa) =(2^(kappa-1) Gamma(kappa))^(-1)
    (u/phi)^\kappa K_\kappa(u/\phi)}
  
  where \eqn{\phi}{phi}  and \eqn{\kappa}{kappa} are parameters and 
  \eqn{K_\kappa(\cdot)}{K_\kappa(\ldots)} denotes the modified Bessel function of the third
  kind of order \eqn{\kappa}{kappa}.
  The family is valid for \eqn{\phi>0}{phi > 0} and \eqn{\kappa>0}{kappa > 0}.
}


\value{
  A vector matrix or array, according to the argument \code{u}, with the
  values of the \eqn{\mbox{Mat\'{e}rn}}{Matern}
  correlation function for the given distances.
}

\references{
  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{cov.spatial}} for the correlation functions
  implemented in \bold{geoR}, and \code{\link{besselK}} for computation
  of the Bessel functions.}

\examples{
dval <- seq(0,1,l=101)
#
# Models with fixed range and varying smoothness parameter
#
mat.5 <-  matern(dval, phi = .25, kappa = 0.5)
mat1 <-  matern(dval, phi = .25, kappa = 1)
mat2 <-  matern(dval, phi = .25, kappa = 2)
mat3 <-  matern(dval, phi = .25, kappa = 3)
plot(0:1, 0:1, type="n", xlab="distance", ylab=expression(rho(h)),
main=expression(paste("varying  ", kappa, "  and fixed  ", phi)))
lines(dval, mat.5, lty=2)
lines(dval, mat1)
lines(dval, mat2, lwd=2, lty=2)
lines(dval, mat3, lwd=2)
legend(0.6,1, c(expression(kappa == 0.5), expression(kappa == 1),
    expression(kappa == 2), expression(kappa == 3)),
    lty=c(2,1,2,1), lwd=c(1,1,2,2))
#
# Correlations with equivalent "practical range"
# and varying smoothness parameter
#
mat.5 <-  matern(dval, phi = .25, kappa = 0.5)
mat1 <-  matern(dval, phi = .188, kappa = 1)
mat2 <-  matern(dval, phi = .140, kappa = 2)
mat3 <-  matern(dval, phi = .117, kappa = 3)
plot(0:1, 0:1, type="n", xlab="distance", ylab=expression(rho(h)),
    main="models with the equivalent \"practical\" range")
lines(dval, mat.5, lty=2)
lines(dval, mat1)
lines(dval, mat2, lwd=2, lty=2)
lines(dval, mat3, lwd=2)
legend(0.5,1, c(expression(paste(kappa == 0.5, "  and  ",
    phi == 0.250)),
    expression(paste(kappa == 1, "  and  ", phi == 0.188)),
    expression(paste(kappa == 2, "  and  ", phi == 0.140)),
    expression(paste(kappa == 3, "  and  ", phi == 0.117))),
    lty=c(2,1,2,1), lwd=c(1,1,2,2))

}

\keyword{spatial}%-- one or more ...

