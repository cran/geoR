\name{varcov.spatial}

\alias{varcov.spatial}

%- Also NEED an `\alias' for EACH other topic documented here.

\title{Computes Covariance Matrix and Related Results}

\description{
  This function builds the covariance matrix for a set of spatial
  locations, given the covariance parameters.
  According to the input options  other results related to the covariance matrix
  (such as decompositions, determinants, inverse. etc)
  can also be returned.
}

\usage{
varcov.spatial(coords = NULL, dists.lowertri = NULL,
               cov.model = "matern", kappa = 0.5, nugget = 0,
               cov.pars = stop("no cov.pars argument"), 
               inv = FALSE, det = FALSE,
               func.inv = c("cholesky", "eigen", "svd", "solve"),
               scaled = FALSE,  only.decomposition = FALSE,
               sqrt.inv = FALSE, try.another.decomposition = TRUE,
               only.inv.lower.diag = FALSE, \dots)

}
	     
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{coords}{an \eqn{n \times 2}{n x 2} matrix with the coordinates
    of the data locations. If not provided
    the argument \code{dists.lowertri} should be provided instead.}
  \item{dists.lowertri}{a vector with the lower triangle of the matrix
    of distances between pairs of data points. If not provided
    the argument \code{coords} should be provided instead.}
  \item{cov.model}{a string indicating the type of the correlation
    function. More details in the
    documentation for \code{\link{cov.spatial}}.
    Defaults are equivalent to the \emph{exponential} model.  }
  \item{kappa}{values of the additional smoothness parameter, only required by
    the following correlation
    functions: \code{"matern"}, \code{"powered.exponential"}, \code{"cauchy"} and
    \code{"gneiting.matern"}.  }
  \item{nugget}{the value of the nugget parameter \eqn{\tau^2}{tau^2}.  }
  \item{cov.pars}{a vector with 2 elements or an \eqn{ns \times 2} matrix with
    the covariance parameters. The first element (if a vector) or first
    column (if a matrix) corresponds to the variance parameter \eqn{\sigma^2}{sigma^2}.
    second element or column corresponds to the correlation function parameter
    \eqn{\phi}{phi}. If a matrix is provided each row corresponds to
    the parameters of one \emph{spatial structure}. Models with several structures
    are also called \emph{nested models} in the geostatistical
    literature.  }
  \item{inv}{if \code{TRUE} the inverse of covariance
    matrix is returned. Defaults to \code{FALSE}.}
  \item{det}{if \code{TRUE} the  logarithmic of the square root of the
    determinant of the covariance
    matrix is returned. Defaults to \code{FALSE}.}
  \item{func.inv}{algorithm used for the decomposition and inversion of the covariance
    matrix. Options are \code{"chol"} for Cholesky decomposition,
    \code{"svd"} for singular value decomposition and \code{"eigen"} for
    eigenvalues/eigenvectors decomposition. Defaults to \code{"chol"}.}
  \item{scaled}{logical indicating whether the covariance matrix should
    be scaled. If \code{TRUE} the partial sill
    parameter \eqn{\sigma^2}{sigma^2} is set to 1. Defaults to \code{FALSE}.}
  \item{only.decomposition}{logical. If \code{TRUE} only the square root
    of the covariance matrix is
    returned. Defaults to \code{FALSE}.  }
  \item{sqrt.inv}{if \code{TRUE} the square root of the inverse of  covariance
    matrix is returned. Defaults to \code{FALSE}.}
  \item{try.another.decomposition}{logical. If \code{TRUE} and the argument
    \code{func.inv} is one of \code{"cholesky"}, \code{"svd"} or
    \code{"solve"}, the matrix decomposition or inversion is tested and,
    if it fails, the argument \code{func.inv} is re-set to \code{"eigen"}.  }
  \item{only.inv.lower.diag}{logical. If \code{TRUE} only the lower triangle and
    the diagonal of the inverse of the covariance matrix are
    returned. Defaults to \code{FALSE}.}
  \item{\dots}{for naw, only for internal usage.}
}


\details{
  The elements of the covariance matrix are computed by the function
  \code{\link{cov.spatial}}. Typically  this is an auxiliary function  called by other
  functions in the \pkg{geoR} package.
}


\value{
  The result is always list. The components will vary according to the input
  options. The possible components are:
  
  \item{varcov}{the covariance matrix.  }
  \item{sqrt.varcov}{a square root of the covariance matrix.  }
  \item{lower.inverse}{the lower triangle of the inverse of covariance
    matrix.  }
  \item{diag.inverse}{the diagonal of the inverse of covariance matrix.
  }
  \item{inverse}{the inverse of covariance matrix.  }
  \item{sqrt.inverse}{a square root of the inverse of covariance
    matrix.  } 
  \item{log.det.to.half}{the logarithmic of the square root of the
    determinant of the covariance matrix.   }
}

\references{
  Further information on the package \pkg{geoR} can be found at:\cr
  \url{http://www.leg.ufpr.br/geoR}.
}

\author{
  Paulo J. Ribeiro Jr. \email{paulojus@leg.ufpr.br}, \cr
  Peter J. Diggle \email{p.diggle@lancaster.ac.uk}.
}

%\note{ ~~further notes~~ }

\seealso{\code{\link{cov.spatial}} for more information on the
  correlation functions; \code{\link{chol}}, \code{\link{solve}},
  \code{\link{svd}} and \code{\link{eigen}} for matrix inversion and/or decomposition.
}

%\examples{
%}

\keyword{spatial}%-- one or more ...

