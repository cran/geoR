\name{cov.spatial}

\alias{cov.spatial}
\alias{.cor.number}
\alias{.check.cov.model}
\alias{geoRCovModels}

\encoding{latin1}

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
cov.spatial(obj, cov.model= "matern",
            cov.pars=stop("no cov.pars argument provided"),
            kappa = 0.5)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj}{a numeric object (vector or matrix), typically with values
    of distances between pairs of spatial locations.  }
  \item{cov.model}{string indicating the type of the correlation
    function. Available choices are: "matern", "exponential", "gaussian",
            "spherical", "circular", "cubic", "wave",
            "power", "powered.exponential", "cauchy", "gencauchy", 
            "gneiting", "gneiting.matern", "pure.nugget".
    See section \code{DETAILS} for available options and expressions of the correlation
    functions.  }
  \item{cov.pars}{a vector with 2 elements or an \eqn{ns \times 2}{ns x 2} matrix with
    the covariance parameters. The first element (if a vector) or first
    column (if a matrix) corresponds to the variance parameter \eqn{\sigma^2}{sigma^2}. The second element or column corresponds to the range parameter
    \eqn{\phi}{phi} of the correlation function.
    If a matrix is provided, each row corresponds to
    the parameters of one \emph{spatial structure} (see DETAILS below).
    }
  \item{kappa}{numerical value for the additional smoothness parameter of the
    correlation function.
    Only required by the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"},
    \code{"cauchy"}, \code{"gencauchy"} and \code{"gneiting.matern"}.  }
}


\details{
  Covariance functions return the value of the covariance
  \eqn{C(h)} between a pair variables located at points separated by the
  distance \eqn{h}.
  The covariance function can be written as a product of a variance
  parameter \eqn{\sigma^2}{sigma^2} times a positive definite
  \emph{correlation function} \eqn{\rho(h)}{rho(h)}:
  \deqn{C(h) = \sigma^2 \rho(h).}{C(h) = sigma^2 * rho(h).}
  The expressions of the covariance functions available in \pkg{geoR}
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
  
  \bold{cauchy} \cr
  \deqn{\rho(h) = [1+(\frac{h}{\phi})^2]^{-\kappa}}{rho(h) = [1+(h/phi)^2]^(-kappa)}

  \bold{gencauchy (generalised Cauchy)} \cr
  \deqn{\rho(h) = [1+(\frac{h}{\phi})^{\kappa_{2}}]^{-{\kappa_1}/{\kappa_2}},
    \kappa_1 > 0, 0 < \kappa_2 \leq 2 }{rho(h) =
    [1+(h/phi)^kappa2]^(-kappa1/kappa2), kappa1 > 0, 0 < kappa2 <=0}
    
  \bold{circular}\cr
  Let \eqn{\theta = \min(\frac{h}{\phi},1)}{theta = min(h/phi,1)} and
  \deqn{g(h)= 2\frac{(\theta\sqrt{1-\theta^2}+
      \sin^{-1}\sqrt{\theta})}{\pi}.}{gamma(h)= 2 * ((\theta *
    \sqrt{1-\theta^2} +
      \sin^{-1} \sqrt{theta}))/pi.}
  Then, the circular model is given by:
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    1 - g(h)  \mbox{ , if $h < \phi$}\cr
    0    \mbox{ , otherwise}
    \end{array} \right.}{rho(h) =
    1 - gamma(h)  if h < phi , 0 otherwise}
  
  \bold{cubic}\cr
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    1 - [7(\frac{h}{\phi})^2 - 8.75(\frac{h}{\phi})^3 +
    3.5(\frac{h}{\phi})^5-0.75(\frac{h}{\phi})^7] \mbox{ , if $h<\phi$} \cr
    0  \mbox{ , otherwise.}
    \end{array} \right.}{rho(h) = 
    1 - (7 * ((h/phi)^2) - 8.75 * ((h/phi)^3) +
    3.5 * ((h/phi)^5) - 0.75 * ((h/phi)^7)) if h < phi
    , 0  otherwise.}
  
  \bold{gaussian}\cr
  \deqn{\rho(h) = \exp[-(\frac{h}{\phi})^2]}{rho(h) = exp(-(h/phi)^2)} 

  \bold{exponential}\cr
  \deqn{\rho(h) = \exp(-\frac{h}{\phi})}{rho(h) = exp(-h/phi)}

  \bold{matern}\cr
  \deqn{\rho(h) =
    \frac{1}{2^{\kappa-1}\Gamma(\kappa)}(\frac{h}{\phi})^\kappa
    K_{\kappa}(\frac{h}{\phi})}{\rho(h) =
    (1/(2^(\kappa-1) * \Gamma(\kappa))) * ((h/\phi)^\kappa) * K_{\kappa}(h/\phi)} 

  \bold{spherical}\cr
  \deqn{\rho(h) = \left\{ \begin{array}{ll}
    1 - 1.5\frac{h}{\phi} + 0.5(\frac{h}{\phi})^3
    \mbox{ , if $h$ < $\phi$} \cr
    0    \mbox{ ,  otherwise}
    \end{array} \right.}{rho(h) = 
    1 - 1.5 * (h/phi) + 0.5(h/phi)^3   if h < phi , 
    0   otherwise}
    
  \bold{power (and linear)}\cr
  The parameters of the this model
  \eqn{\sigma^2}{sigma^2} and \eqn{\phi}{phi} can not be  
  interpreted as \emph{partial sill} and \emph{range}
  as for the other models.
  This model implies an unlimited dispersion and,
  therefore, has no sill and corresponds to a process which is only
  intrinsically stationary.
  The variogram function is given by:
  \deqn{\gamma(h) = \sigma^2 {h}^{\phi} \mbox{ , } 0 < \phi < 2,
    \sigma^2 > 0}{sigma^2
    * h^phi , 0 < phi < 2, sigma^2 > 0} 

  Since the corresponding process is not second order stationary the
  covariance and correlation functions are not defined.
  For internal calculations the \emph{geoR}
  functions uses the fact the this model possesses locally
  stationary representations with covariance functions of the form:
  \deqn{C_(h) = \sigma^2 (A - h^\phi),}{C_(h) = sigma^2 (A - h^phi),}
  where \eqn{A}{A} is a suitable constant as given in
  Chilès \& Delfiner (pag. 511, eq. 7.35).
  
  The \emph{linear} model corresponds a particular case with
  \eqn{\phi = 1}{phi = 1}.
  
  \bold{powered.exponential (or stable)}\cr
  \deqn{\rho(h) = \exp[-(\frac{h}{\phi})^\kappa] \mbox{ , } 0 < \kappa
    \leq 2}{rho(h) = exp[-(h/phi)^kappa] if  0 < kappa <= 2}

  
  \bold{gneiting}\cr
  \deqn{C(h)=\left(1 + 8 sh + 25 (sh)^2 + 32
    (sh)^3\right)(1-sh)^8 1_{[0,1]}(sh)}{C(x)=
    (1 + 8 s h + 25 s^2 h^2 + 32
    s^3 h^3)*(1-s h)^8   if 0<=s h<=1, 0 otherwise}
  where
  \eqn{s=0.301187465825}{}.
  For further details see documentation of the function
  \code{\link[RandomFields]{CovarianceFct}} in the package
  \code{RandomFields} from where we extract the following :\cr
  \emph{It is an alternative to the \code{gaussian} model since
  its graph is visually hardly distinguishable from the graph of
  the Gaussian model, but possesses neither the mathematical and nor the
  numerical disadvantages of the Gaussian model.}\cr


  \bold{gneiting.matern}\cr
  Let \eqn{\alpha=\phi\kappa_2}{alpha=phi * kappa_2}, \eqn{\rho_m(\cdot)}{rho_m(.)} denotes the \eqn{\mbox{Mat\'{e}rn}}{Matern} model
  and \eqn{\rho_g(\cdot)}{rho_g(.)} the Gneiting model. Then the
  \eqn{\mbox{Gneiting-Mat\'{e}rn}}{Gneiting-Matern} is given by
  \deqn{\rho(h) = \rho_g(h|\phi=\alpha) \,
    \rho_m(h|\phi=\phi,\kappa=\kappa_1)}{rho(h) = rho_g(h | phi=alpha) *
    rho_m(h | phi=phi, kappa = kappa_1)}\cr

  \bold{wave}\cr
  \deqn{\rho(h) = \frac{\phi}{h}\sin(\frac{h}{\phi})}{(phi/h) *
    sin(h/phi)} \cr

  \bold{pure.nugget} \cr
  \deqn{\rho(h) = k}{rho(h) = k}\cr
  where k is a constant value. This model corresponds to
  no spatial correlation.

  \bold{Nested models}
  Models with several structures
  usually called \emph{nested models}
  in the geostatistical literature are also allowed.
  In this case the argument \code{cov.pars} takes a matrix and
  \code{cov.model} and \code{lambda} can either have length equal to
  the number of rows of this matrix or length 1.
  For the latter  cov.model and/or lambda are recycled, i.e. the same
  value is used for all structures.

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

  Chilès, J.P. and Delfiner, P. (1999)
  \bold{Geostatistics: Modelling Spatial Uncertainty}, Wiley.

%    For \pkg{geoR} implementation: \cr   
%    Ribeiro, P.J. Jr. and Diggle, P.J. (1999) \emph{geoR: a geostatistical package/library for R/S-PLUS}. Tech. Report ST-99-09, Dept
%  Maths and Stats, Lancaster University.\cr
%  Available on-line at: 
%  \url{http://www.maths.lancs.ac.uk/~ribeiro/publications.html}.
  

  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
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
v.f <- function(x, ...){1-cov.spatial(x, ...)}
#
curve(v.f(x, cov.pars=c(1, .2)), from = 0, to = 1,
      xlab = "distance", ylab = expression(gamma(h)),
      main = "variograms with equivalent \"practical range\"")
curve(v.f(x, cov.pars = c(1, .6), cov.model = "sph"), 0, 1,
      add = TRUE, lty = 2)
curve(v.f(x, cov.pars = c(1, .6/sqrt(3)), cov.model = "gau"),
      0, 1, add = TRUE, lwd = 2)
legend("topleft", c("exponential", "spherical", "gaussian"),
       lty=c(1,2,1), lwd=c(1,1,2))
#
# Matern models with equivalent "practical range"
# and varying smoothness parameter
#
curve(v.f(x, cov.pars = c(1, 0.25), kappa = 0.5),from = 0, to = 1,
      xlab = "distance", ylab = expression(gamma(h)), lty = 2,
      main = "models with equivalent \"practical\" range")
curve(v.f(x, cov.pars = c(1, 0.188), kappa = 1),from = 0, to = 1,
      add = TRUE)      
curve(v.f(x, cov.pars = c(1, 0.14), kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2, lty=2)      
curve(v.f(x, cov.pars = c(1, 0.117), kappa = 2),from = 0, to = 1,
      add = TRUE, lwd=2)      
legend("bottomright",
       expression(list(kappa == 0.5, phi == 0.250), 
         list(kappa == 1, phi == 0.188), list(kappa == 2, phi == 0.140),
         list(kappa == 3, phi == 0.117)), lty=c(2,1,2,1), lwd=c(1,1,2,2))
# plotting a nested variogram model
curve(v.f(x, cov.pars = rbind(c(.4, .2), c(.6,.3)),
          cov.model = c("sph","exp")), 0, 1, ylab='nested model')
}

\keyword{spatial}%-- one or more ...
\keyword{models}%-- one or more ...





