\name{cov.spatial}

\alias{cov.spatial}
\alias{cor.number}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Computes Value of the Covariance Function}

\description{
  Computes the covariances for pairs variables, given the separation
  distance of their locations.
  Options for different correlation functions are available.
  The results can be seen as a change of metric,
  from the \emph{Euclidean distances} to \emph{covariances}.
}

\usage{
cov.spatial(obj, cov.model=c("exponential", "matern", "gaussian",
            "spherical", "circular", "cubic", "wave",
            "powered.exponential", "cauchy", 
            "gneiting", "gneiting.matern", "pure.nugget"),
            cov.pars=stop("no cov.pars argument provided"), kappa)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{a numeric object (vector or matrix), typically with values
    of distances between pairs of spatial locations.  }
  \item{cov.model}{a string indicating the type of the correlation
    function.
    See section \code{DETAILS} for available options and expressions of the correlation
    functions.  }
  \item{cov.pars}{a vector with 2 elements or an \eqn{ns \times 2}{ns x 2} matrix with
    the covariance parameters. The first element (if a vector) or first
    column (if a matrix) corresponds to the variance parameter \eqn{\sigma^2}{sigma^2}. The
    second element or column corresponds to the range parameter
    \eqn{\phi}{phi} of the correlation function.
    If a matrix is provided, each row corresponds to
    the parameters of one \emph{spatial structure}.
    Models with several structures
    are also called \emph{nested models}
    in the geostatistical literature.  }
  \item{kappa}{numerical value for the additional smoothness parameter of the
    correlation function.
    Only required by the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"},
    \code{"gneiting"} and \code{"gneiting.matern"}.  }
}


\details{
  
  Covariance functions return the value of the covariance
  \eqn{C(h)} between a pair variables located at points separated by the
  distance \eqn{h}.
  The covariance function can be written as a product of a variance
  parameter \eqn{\sigma^2}{sigma^2} times a positive definite
  \emph{correlation function} \eqn{\rho(h)}{rho(h)}:
  \deqn{C(h) = \sigma^2 \rho(h).}{C(h) = sigma^2 * rho(h).}
  The expressions of the covariance functions available in \bold{geoR}
  are given below. We recommend the \emph{LaTeX} (and/or the corresponding
  \emph{.dvi}, \emph{.pdf} or \emph{.ps}) version of this document for
  better visualization of the formulas.
  
  Denote \eqn{\phi}{phi} the basic parameter of the correlation
  function and name it the \emph{range parameter}.
  Some of the correlation functions will have an extra parameter
  \eqn{\kappa}{kappa}, the \emph{smoothness parameter}.
  \eqn{K_\kappa(x)} denotes the modified Bessel
  function of the third kind of order \eqn{\kappa}{kappa}. See
  documentation of the function \code{\link{besselK}} for further details.
  In the equations below the functions are valid for  \eqn{\phi>0}{phi >
    0} and \eqn{\kappa>0}{kappa > 0}, unless stated otherwise.
  
  \bold{exponential}\cr
  \deqn{\rho(h) = \exp(-\frac{h}{\phi})}{exp(-h/phi)}

  
  \bold{wave}\cr
  \deqn{\rho(h) = \frac{\phi}{h}\sin(\frac{h}{\phi})}{(phi/h) *
    sin(h/phi)} 

  
  \bold{matern}\cr
  \deqn{\rho(h) =
    \frac{1}{2^{\kappa-1}\Gamma(\kappa)}(\frac{h}{\phi})^\kappa
    K_{\kappa}(\frac{h}{\phi})}{rho(h) =
    (1/(2^(kappa-1) * \Gamma(kappa))) * ((h/phi)^kappa) * K_{kappa}(h/phi)} 

  
  \bold{gaussian}\cr
  \deqn{\rho(h) = \exp[(-\frac{h}{\phi})^2]}{exp(-(h/phi)^2)} 

  
  \bold{spherical}\cr
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    1 - 1.5\frac{h}{\phi} + 0.5(\frac{h}{\phi})^3
    \mbox{ , if $h$ < $\phi$} \cr
    0    \mbox{ ,  otherwise}
    \end{array} \right.}{rho(h) = 
    1 - 1.5 * (h/phi) + 0.5(h/phi)^3   if h < $phi$ , 
    0   otherwise}
    
  
  \bold{circular}\cr
  Let \eqn{\theta = \min(\frac{h}{\phi},1)}{theta = min(h/phi,1)} and
  \deqn{g(h)= 2\frac{(\theta\sqrt{1-\theta^2}+
      \sin^{-1}\sqrt{\theta})}{\pi}.}{gamma(h)= 2 * ((\theta\sqrt{1-\theta^2}+
      \sin^{-1}\sqrt{theta}))/pi.}
  Then, the circular model is given by:
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    1 - g(h)  \mbox{ , if $h < \phi$}\cr
    0    \mbox{ , otherwise}
    \end{array} \right.}{rho(h) =
    1 - gamma(h)  if h < phi , 0 otherwise}

  
  \bold{cubic}\cr
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    7(\frac{h}{\phi})^2-8.75(\frac{h}{\phi})^3+
    3.5(\frac{h}{\phi})^5-0.75(\frac{h}{\phi})^7 \mbox{ , if $h<\phi$} \cr
    0  \mbox{ , otherwise}
    \end{array} \right.}{rho(h) = 
    7 * ((h/phi)^2) - 8.75 * ((h/phi)^3) +
    3.5 * ((h/phi)^5) - 0.75 * ((h/phi)^7) if h < phi
    , 0  otherwise}

  
  \bold{powered.exponential}\cr
  \deqn{\rho(h) = \exp[-(\frac{h}{\phi})^\kappa] \mbox{ , } 0 < \kappa
    \leq 2}{rho(h) = exp[-(h/phi)^kappa] if  0 < \kappa <= 2}

  
  \bold{cauchy} \cr
  \deqn{\rho(h) = [1+(\frac{h}{\phi})^2]^{-\kappa}}{rho(h) = [1+(h/phi)^2]^(-kappa)}

  
  \bold{gneiting}\cr
  Let \eqn{\theta = \min(\frac{h}{\phi},1)}{theta = min((h/phi),1)}.
  The Gneiting model is given by:
  \deqn{\rho(h) = (1+8\theta+25\theta^2+32\theta^3)(1-\theta)^8}{rho(h)
    = (1 + 8 * theta + 25 * (theta^2) + 32 * (theta^3)) * ((1 - theta)^8)}

  
  \bold{gneiting.matern}\cr
  Let \eqn{\alpha=\phi\kappa_2}{alpha=phi * kappa_2}, \eqn{\rho_m(\cdot)}{rho_m(.)} denotes the \eqn{\mbox{Mat\'{e}rn}}{Matern} model
  and \eqn{\rho_g(\cdot)}{rho_g(.)} the Gneiting model. Then the
  \eqn{\mbox{Gneiting-Mat\'{e}rn}}{Gneiting-Matern} is given by
  \deqn{\rho(h) = \rho_g(h|\phi=\alpha) \,
    \rho_m(h|\phi=\phi,\kappa=\kappa_1)}{rho(h) = rho_g(h | phi=alpha) *
    rho_m(h | phi=phi, kappa = kappa_1)}

}


\value{
  The function returns values of the covariances corresponding to  the
  given distances. 
  The type of output is the same as the type of the object provided in the
  argument \code{obj}, typically a vector, matrix or array.
}

\references{

  For a review on correlation functions:\cr
  Schlather, M. (1999) \emph{An introduction to positive definite functions and to unconditional
    simulation of random fields}. Technical report ST 99-10, Dept. of Maths and Statistics,
    Lancaster University. 

%    For \bold{geoR} implementation: \cr   
%    Ribeiro, P.J. Jr. and Diggle, P.J. (1999) \emph{geoR: a geostatistical package/library for R/S-PLUS}. Tech. Report ST-99-09, Dept
%  Maths and Stats, Lancaster University.\cr
%  Available on-line at: 
%  \url{http://www.maths.lancs.ac.uk/~ribeiro/publications.html}.
  

  Further information about \bold{geoR} can be found at:\cr
  \url{http://www.maths.lancs.ac.uk/~ribeiro/geoR.html}.
}

\author{
  Paulo J. Ribeiro Jr. \email{Paulo.Ribeiro@est.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\seealso{\code{\link{matern}} for computation of the
  \eqn{\mbox{Mat\'{e}rn}}{Matern} model, \code{\link{besselK}} for
  computation of the Bessel function  and
  \code{\link{varcov.spatial}} for computations related to the covariance matrix. }

\examples{
#
# Variogram models with the same "practical" range:
#
xval <- seq(0,1,l=101)
vexp <- cov.spatial(xval, cov.pars=c(1, .2))
vsph <- cov.spatial(xval, cov.pars=c(1, .60), cov.model="sph")
vgau <- cov.spatial(xval, cov.pars=c(1, .60/sqrt(3)),
                    cov.model="gau")
plot(0:1, 0:1, type="n", xlab="distance",
     ylab=expression(gamma(h)),
     main="variograms with equivalent \"practical range\"")
lines(xval, (1-vexp))
lines(xval, (1-vsph), lty=2)
lines(xval, (1-vgau), lwd=2)
legend(0.5,.3, c("exponential", "spherical", "gaussian"),
       lty=c(1,2,1), lwd=c(1,1,2))
#
# Matern models with equivalent "practical range"
# and varying smoothness parameter
#
dval <- seq(0,1,l=101)
mat.5 <- cov.spatial(dval, cov.pars = c(1, 0.25), kappa = 0.5)
mat1 <-  cov.spatial(dval, cov.pars = c(1, 0.188), kappa = 1,
                     cov.model="mat")
mat2 <-  cov.spatial(dval, cov.pars = c(1, 0.14), kappa = 2,
                     cov.model="mat")
mat3 <-  cov.spatial(dval, cov.pars = c(1, 0.117), kappa = 3,
                     cov.model="mat")
plot(0:1, 0:1, type="n", xlab="distance",
     ylab=expression(gamma(h)),
     main="models with equivalent \"practical\" range")
lines(dval, 1-mat.5, lty=2)
lines(dval, 1-mat1)
lines(dval, 1-mat2, lwd=2, lty=2)
lines(dval, 1-mat3, lwd=2)
legend(0.5,.3, c(expression(paste(kappa == 0.5, "  and  ",
       phi == 0.250)),
expression(paste(kappa == 1, "  and  ", phi == 0.188)),
expression(paste(kappa == 2, "  and  ", phi == 0.140)),
expression(paste(kappa == 3, "  and  ", phi == 0.117))),
lty=c(2,1,2,1), lwd=c(1,1,2,2))

}


\keyword{spatial}%-- one or more ...





