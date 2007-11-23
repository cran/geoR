\name{coords.aniso}

\alias{coords.aniso}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Geometric Anisotropy Correction}

\description{
  Transforms or back-transforms a set of coordinates according to
  the geometric anisotropy parameters.
}

\usage{
coords.aniso(coords, aniso.pars, reverse = FALSE)
}

\arguments{
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix with the coordinates
    to be transformed.  }
  \item{aniso.pars}{a vector with two elements, \eqn{\psi_A}{psiA} and
    \eqn{\psi_R}{psiR},  the \emph{anisotropy
      angle} and the \emph{anisotropy ratio}, respectively. Notice that the
    parameters must be provided in this order. 
    See section \code{DETAILS} below for more information on anisotropy parameters.  }
  \item{reverse}{logical. Defaults to \code{FALSE}. If \code{TRUE} the reverse
    transformation is performed.  }
}

\details{
  Geometric anisotropy is defined by two parameters:
  \describe{
    \item{Anisotropy angle}{defined here as the azimuth angle of the
    direction with greater spatial continuity, i.e. the angle between the
    \emph{y-axis} and the direction with the maximum range.  }
    
    \item{Anisotropy ratio}{defined here as the ratio between the ranges
    of the directions with greater and smaller continuity, i.e. the ratio
    between maximum and minimum ranges. Therefore, its value is always
    greater or equal to one.  }
  }
    
  If \code{reverse = FALSE} (the default) the
  coordinates are transformed  from the \emph{anisotropic space} to the \emph{isotropic
    space}.
  The transformation consists in multiplying the original
  coordinates by a rotation matrix \eqn{R}
  and a
  shrinking matrix \eqn{T}, as follows:
  \deqn{X_m = X  R T , }{X_m = X \%*\% R \%*\% T , }
  where \eqn{X_m} is a matrix with the modified coordinates (isotropic
  space) , \eqn{X} is a matrix with original coordinates (anisotropic
  space), \eqn{R} rotates coordinates according to the anisotropy angle
  \eqn{\psi_A}{psiA} and \eqn{T} shrinks the coordinates according to
  the anisotropy ratio \eqn{\psi_R}{psiR}.
  
  If \code{reverse = TRUE}, the back-transformation is performed, i.e.
  transforming the coordinates from the \emph{isotropic space} to the
  \emph{anisotropic  space} by computing:
  \deqn{X = X_m (R T)^{-1}}{X = X_m \%*\% solve( R \%*\% T ).}

}


\value{
  An \eqn{n \times 2}{n x 2} matrix with the transformed coordinates.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br} \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

 
%\seealso{ ~~objects to SEE ALSO as \code{\link{~~fun~~}}, ~~~ }

\examples{
op <- par(no.readonly = TRUE)
par(mfrow=c(3,2))
par(mar=c(2.5,0,0,0))
par(mgp=c(2,.5,0))
par(pty="s")
## Defining a set of coordinates
coords <- expand.grid(seq(-1, 1, l=3), seq(-1, 1, l=5))
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type="n")
text(coords[,1], coords[,2], 1:nrow(coords))
## Transforming coordinates according to some anisotropy parameters
coordsA <- coords.aniso(coords, aniso.pars=c(0, 2))
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type="n")
text(coordsA[,1], coordsA[,2], 1:nrow(coords))
##
coordsB <- coords.aniso(coords, aniso.pars=c(pi/2, 2))
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type="n")
text(coordsB[,1], coordsB[,2], 1:nrow(coords))
##
coordsC <- coords.aniso(coords, aniso.pars=c(pi/4, 2))
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type="n")
text(coordsC[,1], coordsC[,2], 1:nrow(coords))
##
coordsD <- coords.aniso(coords, aniso.pars=c(3*pi/4, 2))
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type="n")
text(coordsD[,1], coordsD[,2], 1:nrow(coords))
##
coordsE <- coords.aniso(coords, aniso.pars=c(0, 5))
plot(c(-1.5, 1.5), c(-1.5, 1.5), xlab="", ylab="", type="n")
text(coordsE[,1], coordsE[,2], 1:nrow(coords))
##
par(op)
}

\keyword{spatial}
