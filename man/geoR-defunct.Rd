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
  The functions listed here are no longer part of the package \pkg{geoR}
  as they are no longer needed.
}

\details{
  
  The following functions are now defunct:
  
  \describe{
    
    \item{olsfit}{functionality incorporated by \code{\link{variofit}}
      starting from package version \samp{1.0-6}.  }
    
    \item{wlsfit}{functionality incorporated by \code{\link{variofit}}
      starting from package version \samp{1.0-6}.  }
    
    \item{likfit.old}{functionality incorporated by \code{\link{likfit}}
      starting from package version \samp{1.0-6}.  
      The related functions were also made defunct: \cr
      \code{likfit.nospatial}, \code{loglik.spatial},
      \code{proflik.nug}, \code{proflik.phi}, \code{proflik.ftau}.  }
    
    \item{distdiag}{functionally is redundant with \code{\link{dist}}.  }
  }
  
  
}

\seealso{\code{\link{variofit}}}

\keyword{spatial}
\keyword{utilities}

