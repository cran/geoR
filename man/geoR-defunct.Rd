\name{geoR-defunct}

\title{Defunct Functions in the Package geoR}

\alias{geoRdefunct}
\alias{olsfit}
\alias{wlsfit}
\alias{likfit.old}
\alias{likfit.nospatial}
\alias{loglik.spatial}
\alias{.proflik.nug}
\alias{.proflik.phi}
\alias{.proflik.ftau}
\alias{distdiag}

\usage{
geoRdefunct()
}

\description{
  The functions listed here are no longer part of the package \pkg{geoR} as they are not needed (any more).
}

\details{
  
  The following functions are now defunct:
  
  \enumerate{
    
    \item{olsfit}{functionality incorporated by \code{\link{variofit}}.
      From \code{geoR_1.0-6}.  }
    
    \item{wlsfit}{functionality incorporated by \code{\link{variofit}}.
      From \code{geoR_1.0-6}.  }
    
    \item{likfit.old}{functionality incorporated by \code{\link{likfit}}.
      From \code{geoR_1.0-6}. The associated functions were also made
      defunct: \cr
      \code{likfit.nospatial}, \code{loglik.spatial}, \code{proflik.nug}, \code{proflik.phi}, \code{proflik.ftau}.  }

    \item{distdiag}{functionally is redundant with
      \code{\link{dist}}.  }
  }


}

\seealso{\code{\link{variofit}}}

\keyword{spatial}
\keyword{utilities}

