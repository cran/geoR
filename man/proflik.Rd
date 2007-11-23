\name{proflik}

\concept{profile likelihood}

\alias{proflik}

\alias{.proflik.aux0}
\alias{.proflik.aux1}
\alias{.proflik.aux10}
\alias{.proflik.aux11}
\alias{.proflik.aux1.1}
\alias{.proflik.aux12}
\alias{.proflik.aux13}
\alias{.proflik.aux14}
\alias{.proflik.aux15}
\alias{.proflik.aux16}
\alias{.proflik.aux17}
\alias{.proflik.aux18}
\alias{.proflik.aux19}
\alias{.proflik.aux2}
\alias{.proflik.aux20}
\alias{.proflik.aux21}
\alias{.proflik.aux21.1}
\alias{.proflik.aux22}
\alias{.proflik.aux23}
\alias{.proflik.aux24}
\alias{.proflik.aux27}
\alias{.proflik.aux28}
\alias{.proflik.aux30}
\alias{.proflik.aux3}
\alias{.proflik.aux31}
\alias{.proflik.aux32}
\alias{.proflik.aux33}
\alias{.proflik.aux4}
\alias{.proflik.aux5}
\alias{.proflik.aux6}
\alias{.proflik.aux7}
\alias{.proflik.aux8}
\alias{.proflik.aux9}
\alias{.proflik.cov}
\alias{.proflik.lambda}
\alias{.proflik.main}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Computes Profile Likelihoods}

\description{
  Computes profile likelihoods for model parameters
  previously estimated using the function
  \code{\link{likfit}}. 
}

\usage{
proflik(obj.likfit, geodata, coords = geodata$coords,
        data = geodata$data, sill.values, range.values,
        nugget.values, nugget.rel.values, lambda.values, 
        sillrange.values = TRUE, sillnugget.values = TRUE,
        rangenugget.values = TRUE, sillnugget.rel.values = FALSE,
        rangenugget.rel.values = FALSE, silllambda.values = FALSE,
        rangelambda.values = TRUE,  nuggetlambda.values = FALSE,
        nugget.rellambda.values = FALSE,
        uni.only = TRUE, bi.only = FALSE, messages, \dots)
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{obj.likfit}{an object of the class \code{likfit}, typically an
    output of the function \code{\link{likfit}}.}
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} described next. Typically an object of the class
    \code{"geodata"} - a \pkg{geoR} data-set. If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{ an \eqn{n \times 2}{n x 2} matrix containing in each row Euclidean
    coordinates of the \eqn{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
  \item{sill.values}{set of values of the partial sill parameter
    \eqn{\sigma^2}{sigma^2} for which the profile likelihood will be
    computed.  }
  \item{range.values}{set of values of the range parameter
    \eqn{\phi}{phi} for which the profile likelihood will be computed.  }
  \item{nugget.values}{set of values of the nugget parameter
    \eqn{\tau^2}{tau^2} for which the profile likelihood will be
    computed. Only used if the model was fitted using the function
    \code{\link{likfit}} with the option \code{fix.nugget = FALSE}.  }
  \item{nugget.rel.values}{set of values of the relative nugget parameter
    \eqn{\tau_{R}^{2}}{tauR^2} for which the profile likelihood will be
    computed. Only used if the model was fitted using the function
    \code{\link{likfit}} with the option \code{fix.nugget = FALSE}.  }
  \item{lambda.values}{set of values of the Box-Cox transformation parameter
    \eqn{\lambda}{lambda} for which the profile likelihood will be
    computed. Only to be used if the model was fitted using the function
    \code{\link{likfit}} with the option \code{fix.lambda = FALSE}.  }
  \item{sillrange.values}{logical indicating
    whether or not the 2-D profile likelihood should be computed.
    Only valid if \code{uni.only = FALSE}.  }
  \item{sillnugget.values}{as above.  }
  \item{rangenugget.values}{as above.  }
  \item{sillnugget.rel.values}{as above.  }
  \item{rangenugget.rel.values}{as above.  }
  \item{silllambda.values}{as above.  }
  \item{rangelambda.values}{as above.  }
  \item{nuggetlambda.values}{as above.  }
  \item{nugget.rellambda.values}{as above.  }
  \item{uni.only}{as above.  }
  \item{bi.only}{as above.  }
  \item{messages}{logical. Indicates whether
    status messages should be printed on the screen (i.e. current output device)
    while the function is running.  }
%  \item{minimisation.function}{minimization function to be
%    used. Defaults to \code{\link{optim}}.  }
  \item{\dots}{additional parameters to be passed to the minimization
    function.  }
}

\details{
  The functions \code{.proflik.*} are auxiliary functions used to
  compute the profile likelihoods. These functions are
  internally called by the
  minimization functions when estimating the model parameters.
}

\value{
  An object of the class \code{"proflik"} which is 
  a list. Each element contains values of a parameter (or a pair of
  parameters for 2-D profiles) and the
  corresponding value of the profile likelihood.
  The components of the output will vary according to the
  input options.
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

\note{
  \enumerate{
    \item Profile likelihoods for Gaussian Random Fields are usually
    uni-modal. 
    Unusual or jagged shapes can
    be due to the lack of the convergence in the numerical minimization
    for particular values of the parameter(s).  
    If this is the case it might be necessary to pass \code{control} arguments 
    to the minimization functions using the argument \dots.
    It's also advisable to try the different options for the
    \code{minimisation.function} argument.
    See documentation of the functions \code{\link{optim}} and/or
    \code{\link{nlm}} for further details.
    
    \item 2-D profiles can be computed by setting the argument
    \code{uni.only = FALSE}. However, before computing 2-D profiles be
    sure they are really necessary.
    Their computation can be time demanding since it 
    is performed on a grid determined by the
    cross-product of the values defining the 1-D profiles.
        
    \item There is no "default strategy" to find reasonable values for the
    x-axis.
    They must be found in a "try-and-error" exercise. It's recommended
    to use short sequences in the initial attempts.
    The \code{EXAMPLE} section below illustrates this.
  }

}

\seealso{\code{\link{plot.proflik}} for graphical output,
  \code{\link{likfit}} for the parameter estimation,
  \code{\link{optim}} and \code{\link{nlm}} for further details about
  the minimization functions.  }

\examples{
op <- par(no.readonly=TRUE)
ml <- likfit(s100, ini=c(.5, .5), fix.nug=TRUE)
## a first atempt to find reasonable values for the x-axis:
prof <- proflik(ml, s100, sill.values=seq(0.5, 1.5, l=4),
                range.val=seq(0.1, .5, l=4))
par(mfrow=c(1,2))
plot(prof)
## a nicer setting 
\dontrun{
prof <- proflik(ml, s100, sill.values=seq(0.45, 2, l=11),
                range.val=seq(0.1, .55, l=11))
plot(prof)
## to include 2-D profiles use:
## (commented because this is time demanding)
#prof <- proflik(ml, s100, sill.values=seq(0.45, 2, l=11),
#                range.val=seq(0.1, .55, l=11), uni.only=FALSE)
#par(mfrow=c(2,2))
#plot(prof, nlevels=16)
}
par(op)
}

\keyword{spatial}
