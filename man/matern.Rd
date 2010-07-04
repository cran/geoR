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

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }


\seealso{\code{\link{cov.spatial}} for the correlation functions
  implemented in \pkg{geoR}, and \code{\link{besselK}} for computation
  of the Bessel functions.}

\examples{
#
# Models with fixed range and varying smoothness parameter
#
curve(matern(x, phi= 0.25, kappa = 0.5),from = 0, to = 1.5,
      xlab = "distance", ylab = expression(rho(h)), lty = 2,
      main=expression(paste("varying  ", kappa, "  and fixed  ", phi)))
curve(matern(x, phi= 0.25, kappa = 1),from = 0, to = 1.5, add = TRUE)
curve(matern(x, phi= 0.25, kappa = 2),from = 0, to = 1.5, add = TRUE,
      lwd = 2, lty=2)
curve(matern(x, phi= 0.25, kappa = 3),from = 0, to = 1.5, add = TRUE,
      lwd = 2)
legend("topright", expression(kappa==0.5, kappa==1.5, kappa==2, kappa==3),
    lty=c(2,1,2,1), lwd=c(1,1,2,2))

#
# Correlations with equivalent "practical range"
# and varying smoothness parameter
#
curve(matern(x, phi = 0.25, kappa = 0.5),from = 0, to = 1,
      xlab = "distance", ylab = expression(gamma(h)), lty = 2,
      main = "models with equivalent \"practical\" range")
curve(matern(x, phi = 0.188, kappa = 1),from = 0, to = 1, add = TRUE)      
curve(matern(x, phi = 0.14, kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2, lty=2)      
curve(matern(x, phi = 0.117, kappa = 2), from = 0, to = 1,
      add = TRUE, lwd=2)      
legend("topright", expression(list(kappa == 0.5, phi == 0.250),
       list(kappa == 1, phi == 0.188), list(kappa == 2, phi == 0.140),
       list(kappa == 3, phi == 0.117)), lty=c(2,1,2,1), lwd=c(1,1,2,2))
}

\keyword{spatial}%-- one or more ...

