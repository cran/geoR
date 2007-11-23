\name{variog.mc.env}

\alias{variog.mc.env}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Envelops for Empirical Variograms Based on Permutation}

\description{
  Computes envelops for empirical variograms by
  permutation of the data values on the spatial locations.
}

\usage{
variog.mc.env(geodata, coords = geodata$coords, data = geodata$data,
              obj.variog, nsim = 99, save.sim = FALSE, messages) 
}

%- maybe also `usage' for other objects documented here.

\arguments{
  \item{geodata}{a list containing elements \code{coords} and
    \code{data} as described next. 
    Typically an object of the class \code{"geodata"} - a \pkg{geoR}
    data-set.
    If not provided the arguments
    \code{coords} and \code{data} must be provided instead.  }
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix, each row containing Euclidean
    coordinates of the \emph{n} data locations. By default it takes the
    element \code{coords} of the argument \code{geodata}.  }
  \item{data}{a vector with the data values. By default it takes the
    element \code{data} of the argument \code{geodata}.  }
  \item{obj.variog}{an object of the class \code{"variogram"}, typically
    an output of the function \code{\link{variog}}.  }
  \item{nsim}{number of simulations used to compute the envelope. Defaults to
    99.  }
  \item{save.sim}{logical. Indicates whether or not the simulated data
    are included in the output. Defaults to \code{FALSE}.  }
  \item{messages}{logical. If \code{TRUE}, the default, status messages
    are printed while the function is running.  }
}

\details{
  The envelops are obtained by permutation.
  For each simulations data values are randomly allocated to the spatial locations.
  The empirical variogram is computed for each simulation
  using the same lags as for the variogram originally computed for the data.
  The envelops are computed by taking, at each lag, the maximum and minimum
  values of the variograms for the simulated data.
}

\value{
  An object of the \code{\link{class}} \code{"variogram.envelope"} which is a
  list with the following components:
  \item{u}{a vector with distances. }
  \item{v.lower}{a vector with the minimum variogram values at each
    distance in \code{u}.  }
  \item{v.upper}{a vector with the maximum variogram values at each
    distance in \code{u}.  }
  \item{simulations}{a matrix with simulated data.
    Only returned if \code{save.sim = TRUE}.  }
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}


\author{Paulo Justiniano Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{variog.model.env}} for envelops computed by
  from a model specification, 
   \code{\link{variog}} for variogram calculations,
  \code{\link{plot.variogram}} and \code{\link{variog.mc.env}} for
  graphical output.
}

\examples{
s100.vario <- variog(s100, max.dist=1)
s100.env <- variog.mc.env(s100, obj.var = s100.vario)
plot(s100.vario, envelope = s100.env)
}

\keyword{spatial}
\keyword{nonparametric}


